"""Experiment driver: trains fold models, simulates difficulty-model training data, evaluates target systems."""
from __future__ import annotations
import os, json, time
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .protocol import Data, FoldModels, max_similarity, simulate_calibration_rows, FOLD_NAMES
from .system_eval import evaluate_system, SIM_BINS
from ..models.base_models import RetentionOrderModel, SMRTRegressor
from ..models.difficulty import DifficultyModel, GROUPS, encode_system
from ..utils import get_logger, Timer

ABLATIONS = {"P": ["compound", "ion", "novelty", "geometry", "system"],
             "P-novelty": ["compound", "ion", "geometry", "system"],
             "P-ion": ["compound", "novelty", "geometry", "system"],
             "P-system": ["compound", "ion", "novelty", "geometry"],
             "P-geometry": ["compound", "ion", "novelty", "system"],
             "P-compound": ["novelty", "geometry", "system"]}
STRATEGY_CODE = {"random": 0, "span": 1, "common": 2}


def _compound_frame(data: Data, rows, maxsim, seen):
    C = data.fs.physchem(rows)
    C["maxsim"] = maxsim; C["seen"] = seen.astype(np.float32)
    return C


def _strata_labels(data: Data, ik, y, maxsim, seen, C, top_superclasses):
    lab = {}
    lab["seen"] = np.where(seen, "seen", "unseen")
    sb = np.empty(len(ik), dtype=object)
    for lo, hi, name in SIM_BINS:
        sb[(maxsim >= lo) & (maxsim < hi)] = name
    lab["similarity"] = sb
    ion = (C["n_acidic_groups"].values + C["n_basic_groups"].values) > 0
    lab["ionisable"] = np.where(ion, "ionisable", "neutral")
    span = y.max() - y.min(); rel = (y - y.min()) / max(span, 1e-9)
    lab["rt_tertile"] = np.where(rel < 1 / 3, "early", np.where(rel < 2 / 3, "middle", "late"))
    sc = data.superclass.reindex(ik).fillna("Unknown").astype(str).str.replace(r" \(CHEMONTID.*", "", regex=True).values
    lab["superclass"] = np.where(np.isin(sc, top_superclasses), sc, "Other")
    return lab


class Experiment:
    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.out = cfg["out_dir"]; os.makedirs(self.out, exist_ok=True)
        self.log = get_logger(os.path.join(self.out, "experiment_logs", f"{cfg['name']}.log"), name=cfg["name"])
        self.data = Data(cfg["processed_dir"], method_type=cfg.get("method_type", "RP"))
        self.fm = FoldModels(self.data, os.path.join(cfg["cache_dir"], "base_models"), cfg["base_model"], cfg["smrt_n"],
                             cfg["n_jobs"], self.log)
        self.formula_groups = self.data.fs.comp.groupby("formula").indices
        sc = self.data.superclass.astype(str).str.replace(r" \(CHEMONTID.*", "", regex=True)
        tgt_ik = self.data.rt[self.data.rt.system.isin(self.data.targets(cfg["min_target_compounds"]).index)].ik14.unique()
        self.top_superclasses = sc.reindex(tgt_ik).value_counts().head(8).index.tolist()
        self.timings = {}

    # ------------------------------------------------------------------
    def scores_and_novelty(self, model_pair, systems, extra_rows=None):
        model, train_ik = model_pair
        tab = self.data.rt[self.data.rt.system.isin(systems)]
        rows = np.unique(np.r_[self.data.fs.rows(tab.ik14.unique()), extra_rows if extra_rows is not None else []]).astype(int)
        sc = model.predict(self.data.fs, rows)
        ref = self.data.fs.rows(train_ik)
        t = Timer(); ms = max_similarity(self.data, rows, ref, n_jobs=self.cfg["n_jobs"]); self.timings.setdefault("maxsim_s", []).append(t())
        seen = np.isin(self.data.fs.comp.ik14.values[rows], train_ik)
        ms = np.where(seen, 1.0, ms)
        return pd.DataFrame({"row": rows, "score": sc, "maxsim": ms, "seen": seen}).set_index("row")

    def difficulty_training_rows(self, f: int, base_kind: str):
        cfg = self.cfg; data = self.data
        path = os.path.join(cfg["cache_dir"], f"sigma_rows_fold{f}_{base_kind}_t{cfg['base_model']['n_estimators']}.parquet")
        if os.path.exists(path):
            return pd.read_parquet(path)
        parts = []
        for g in sorted(set(data.meta.fold.unique()) - {f}):
            systems = [s for s in data.meta[data.meta.fold == g].index if s != data.smrt_id]
            if base_kind == "order":
                mp = self.fm.get(frozenset({f, g}))
            else:
                mp = self.smrt_model
            nov = self.scores_and_novelty(mp, systems)
            scores, novelty = {}, {}
            for sid in systems:
                tab = data.rt[data.rt.system == sid]; r = data.fs.rows(tab.ik14)
                scores[sid] = nov.loc[r, "score"].values
                novelty[sid] = (nov.loc[r, "maxsim"].values, nov.loc[r, "seen"].values.astype(np.float32))
            parts.append(simulate_calibration_rows(data, systems, scores, novelty, n_list=tuple(cfg["sim_n_list"]),
                                                   draws=cfg["sim_draws"], seed=cfg.get("seed", 0)))
            self.log.info(f"fold {f}: simulated difficulty rows from fold {g}: {len(parts[-1])}")
        df = pd.concat(parts, ignore_index=True); df.to_parquet(path)
        return df

    def train_sigma_models(self, f: int, base_kind: str, which=None):
        df = self.difficulty_training_rows(f, base_kind)
        models = {}
        for name, groups in ABLATIONS.items():
            if which is not None and name not in which:
                continue
            t = Timer()
            models[name] = DifficultyModel(groups=groups, n_jobs=self.cfg["n_jobs"], **self.cfg["difficulty_model"]).fit(df, df["target"].values)
            self.timings.setdefault("sigma_train_s", []).append(t())
            models[name].model.set_params(n_jobs=1)
        return models, len(df)

    def ensemble_models(self, f: int):
        cfg = self.cfg; ens = []
        for e in range(cfg["n_ensemble"]):
            path = os.path.join(cfg["cache_dir"], "base_models", f"ens_fold{f}_rep{e}_t{cfg['ens_trees']}.joblib")
            import joblib
            if os.path.exists(path):
                ens.append(joblib.load(path)); continue
            tr = self.data.training_table({f}, cfg["smrt_n"])
            sys_ = tr.system.unique(); rng = np.random.default_rng([f, e, 7])
            pick = rng.choice(sys_, len(sys_), replace=True)
            parts = [tr[tr.system == s].assign(system=f"{s}_{k}") for k, s in enumerate(pick)]
            kw = dict(self.cfg["base_model"]); kw["n_estimators"] = cfg["ens_trees"]; kw["seed"] = e + 1
            m = RetentionOrderModel(n_jobs=cfg["n_jobs"], **kw).fit(self.data.fs, pd.concat(parts))
            joblib.dump(m, path); ens.append(m)
        return ens

    def decoys_for(self, ik, rng_seed: int):
        fs = self.data.fs; comp = fs.comp
        parents, drows = [], []
        rng = np.random.default_rng(rng_seed)
        for j, k in enumerate(ik):
            frm = self.data.formula.get(k)
            if frm is None or frm not in self.formula_groups:
                continue
            cand = self.formula_groups[frm]
            cand = cand[comp.ik14.values[cand] != k]
            if len(cand) == 0:
                continue
            if len(cand) > self.cfg["max_decoys"]:
                cand = rng.choice(cand, self.cfg["max_decoys"], replace=False)
            parents.extend([j] * len(cand)); drows.extend(cand.tolist())
        return np.array(parents, int), np.array(drows, int)

    # ------------------------------------------------------------------
    def run(self, base_kind: str = "order", strategy: str = "random", perturb=None, folds=None, systems_limit=None,
            sigma_which=None, methods=None, tag="main"):
        cfg = self.cfg; data = self.data
        targets = data.targets(cfg["min_target_compounds"])
        if base_kind == "smrt" and not hasattr(self, "smrt_model"):
            sm = data.rt[data.rt.system == data.smrt_id]
            t = Timer(); m = SMRTRegressor(n_jobs=cfg["n_jobs"], **cfg["smrt_model"]).fit(data.fs, sm)
            self.timings["smrt_train_s"] = t(); self.smrt_model = (m, np.unique(sm.ik14.values))
        results = {"rows": [], "strata": [], "decoys": [], "errors": [], "timing": []}
        for f in (folds if folds is not None else sorted(targets.fold.unique())):
            tsys = list(targets[targets.fold == f].index)
            if systems_limit:
                tsys = tsys[:systems_limit]
            if not tsys:
                continue
            t = Timer()
            if base_kind == "order":
                outer = self.fm.get(frozenset({f}))
            else:
                outer = self.smrt_model
            sig_models, n_sig_rows = self.train_sigma_models(f, base_kind, which=sigma_which)
            ens = self.ensemble_models(f) if (base_kind == "order" and "B6_jk+_ensemble" in (methods or cfg["methods"])) else None
            # decoys
            dec_info = {}
            extra = []
            for sid in tsys:
                tab = data.rt[data.rt.system == sid]
                par, drows = self.decoys_for(tab.ik14.values, int(sid))
                dec_info[sid] = (par, drows); extra.append(drows)
            extra = np.unique(np.concatenate(extra)) if extra else np.array([], int)
            nov = self.scores_and_novelty(outer, tsys, extra_rows=extra if cfg["decoy_draws"] > 0 else None)
            ens_scores = None
            if ens is not None:
                allrows = nov.index.values
                ens_scores = pd.DataFrame({e: m.predict(data.fs, allrows) for e, m in enumerate(ens)}, index=allrows)
            self.log.info(f"[{tag}] fold {f} ({FOLD_NAMES[f]}): {len(tsys)} targets, sigma rows {n_sig_rows}, setup {t():.0f}s")
            jobs = []
            for sid in tsys:
                tab = data.rt[data.rt.system == sid].reset_index(drop=True)
                rows = data.fs.rows(tab.ik14)
                ms = nov.loc[rows, "maxsim"].values; seen = nov.loc[rows, "seen"].values
                C = _compound_frame(data, rows, ms, seen)
                job = {"sid": sid, "y": tab.rt.values.astype(float), "s": nov.loc[rows, "score"].values, "C": C, "ik": tab.ik14.values,
                       "sysf": encode_system(data.meta.loc[sid].to_dict()), "sig_models": sig_models,
                       "ens": ens_scores.loc[rows].values.T if ens_scores is not None else None,
                       "strata": _strata_labels(data, tab.ik14.values, tab.rt.values, ms, seen, C, self.top_superclasses),
                       "n_systems_global": data.n_systems.reindex(tab.ik14.values).fillna(1).values,
                       "cfg": {"alphas": cfg["alphas"], "n_list": cfg["n_list"], "draws": cfg["draws"], "primary_n": cfg["primary_n"],
                               "min_test": cfg["min_test"], "strategy": strategy, "strategy_code": STRATEGY_CODE[strategy],
                               "perturb": perturb, "methods": methods or cfg["methods"], "expensive_n": cfg["expensive_n"],
                               "decoy_draws": cfg["decoy_draws"] if strategy == "random" and perturb is None else 0,
                               "strata": True, "seed": cfg.get("seed", 0)}}
                par, drows = dec_info[sid]
                if len(drows) and job["cfg"]["decoy_draws"] > 0:
                    job["decoys"] = {"parent": par, "s": nov.loc[drows, "score"].values,
                                     "C": _compound_frame(data, drows, nov.loc[drows, "maxsim"].values, nov.loc[drows, "seen"].values)}
                jobs.append(job)
            t = Timer()
            outs = Parallel(n_jobs=cfg["n_jobs"], backend="loky")(delayed(evaluate_system)(j) for j in jobs)
            self.log.info(f"[{tag}] fold {f}: evaluated {len(jobs)} systems in {t():.0f}s")
            for o in outs:
                for k in results:
                    for rr in o[k]:
                        rr["fold"] = f
                    results[k].extend(o[k])
        res = {k: pd.DataFrame(v) for k, v in results.items()}
        for k, df in res.items():
            if len(df):
                df["base"] = base_kind; df["strategy"] = strategy; df["perturb"] = perturb or "none"; df["experiment"] = tag
        return res


def save_results(res: dict, out_dir: str, tag: str):
    os.makedirs(out_dir, exist_ok=True)
    for k, df in res.items():
        if len(df):
            df.to_parquet(os.path.join(out_dir, f"{tag}_{k}.parquet"))

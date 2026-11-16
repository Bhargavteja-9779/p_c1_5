"""Addendum D: selection-aware system-level calibration (SASC) — learn c on other systems, apply to targets and NORMAN.
Usage: python experiments/run_sasc.py [config]"""
import sys, os, json, yaml, time, zlib
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
from joblib import Parallel, delayed
from rdkit import Chem
from conformrt.evaluation.driver import Experiment, STRATEGY_CODE, _compound_frame
from conformrt.evaluation.system_eval import select_calibrants, _sigma_matrix, _sigma_cal_loo
from conformrt.evaluation.sasc import loo_quantities, windows, coverage_curve, choose_c, GRID
from conformrt.evaluation.intervals import winkler
from conformrt.models.difficulty import DifficultyModel, encode_system
from conformrt.preprocessing.features import mol_from_smiles

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/full.yaml"))
quick = cfg["name"] == "quick"; N_SIM_DRAWS = 2 if quick else 5; N_EVAL_DRAWS = 2 if quick else 10; n = cfg["primary_n"]; ALPHA = 0.1
ex = Experiment(dict(cfg, name=cfg["name"] + "_sasc")); data = ex.data; fs = data.fs
RULES = ["random", "span", "common"]
# NORMAN calibration chemicals (connectivity InChIKey blocks)
nd = pd.read_csv("data/raw/norman_ilc/extracted/SI/data/NORMAN_interlab_RT_data.csv")
ncal = nd[nd.sample_type == "cal"].drop_duplicates("compound")
NORMAN_CAL_IK = {Chem.MolToInchiKey(m)[:14] for m in (mol_from_smiles(s) for s in ncal.SMILES) if m is not None}


def sigma_fn(model, C, sysf):
    def f(s_te, s_cal, y_cal, te_idx=None, cal_idx=None):
        return (_sigma_matrix(model, C.iloc[te_idx].reset_index(drop=True), s_te, s_cal, y_cal, sysf, loo=True),
                _sigma_cal_loo(model, C.iloc[cal_idx].reset_index(drop=True), s_cal, y_cal, sysf))
    return f


def sim_system(job):
    """Coverage curves for one training system under each rule (P and J+)."""
    y, s, C, sysf, sig, sid = job["y"], job["s"], job["C"], job["sysf"], job["sig"], job["sid"]
    N = len(y); out = {}
    rules = list(RULES) + (["norman"] if job["norman_pool"] is not None else [])
    for rule in rules:
        curves_P, curves_J = [], []
        draws = 1 if rule == "norman" else N_SIM_DRAWS
        for r in range(draws):
            if rule == "norman":
                cal = job["norman_pool"]
            else:
                if N < n + 20: break
                rng = np.random.default_rng([int(sid), n, r, STRATEGY_CODE[rule], 99]); cal = select_calibrants(rng, y, n, rule, job["nsys"])
            te = np.setdiff1d(np.arange(N), cal)
            if len(te) < 10: continue
            mu, R, sl = loo_quantities(s[cal], y[cal], s[te], lambda st, sc, yc: sigma_fn(sig, C, sysf)(st, sc, yc, te, cal))
            curves_P.append(coverage_curve(mu, R, sl, y[te], ALPHA))
            muJ, RJ, _ = loo_quantities(s[cal], y[cal], s[te], None)
            curves_J.append(coverage_curve(muJ, RJ, None, y[te], ALPHA))
        if curves_P:
            out[rule] = (np.mean(curves_P, 0), np.mean(curves_J, 0))
    return sid, out


def make_job(sid, scores_df, sig, extra=None):
    tab = data.rt[data.rt.system == sid].reset_index(drop=True); rows = fs.rows(tab.ik14)
    C = _compound_frame(data, rows, scores_df.loc[rows, "maxsim"].values, scores_df.loc[rows, "seen"].values)
    pool = np.flatnonzero(np.isin(tab.ik14.values, list(NORMAN_CAL_IK)))
    return {"sid": sid, "y": tab.rt.values.astype(float), "s": scores_df.loc[rows, "score"].values, "C": C,
            "sysf": encode_system(data.meta.loc[sid].to_dict()), "sig": sig, "nsys": data.n_systems.reindex(tab.ik14.values).fillna(1).values,
            "norman_pool": pool if len(pool) >= 10 and len(tab) - len(pool) >= 10 else None}


results_c, eval_rows = [], []
norman_curves_P, norman_curves_J = [], []
targets = data.targets(cfg["min_target_compounds"])
for f in sorted(targets.fold.unique()):
    t0 = time.time()
    rows_f = ex.difficulty_training_rows(f, "order")
    rows_f["sfold"] = data.meta.loc[rows_f.system.astype(str), "fold"].values
    jobs = []
    for g in sorted(set(data.meta.fold.unique()) - {f}):
        sig_fg = DifficultyModel(n_jobs=cfg["n_jobs"], **cfg["difficulty_model"]).fit(rows_f[rows_f.sfold != g], rows_f[rows_f.sfold != g]["target"].values)
        sig_fg.model.set_params(n_jobs=1)
        systems = [s for s in data.meta[(data.meta.fold == g) & (data.meta.n_compounds >= 40)].index if s != data.smrt_id]
        if quick: systems = systems[:4]
        sc = ex.scores_and_novelty(ex.fm.get(frozenset({f, g})), systems)
        jobs += [make_job(sid, sc, sig_fg) for sid in systems]
    outs = Parallel(n_jobs=cfg["n_jobs"])(delayed(sim_system)(j) for j in jobs)
    for rule in RULES:
        cp = np.array([o[rule][0] for _, o in outs if rule in o]); cj = np.array([o[rule][1] for _, o in outs if rule in o])
        results_c.append({"fold": int(f), "rule": rule, "n_train_systems": len(cp), "c_P": choose_c(cp, ALPHA), "c_J": choose_c(cj, ALPHA)})
    print(f"fold {f}: simulation on {len(jobs)} training systems in {time.time()-t0:.0f}s", results_c[-3:], flush=True)
    # ---- apply to target systems of fold f (outer models; same draws as robustness runs) ----
    sig_f = ex.train_sigma_models(f, "order", which=["P"])[0]["P"]
    tsys = list(targets[targets.fold == f].index)[: (2 if quick else None)]
    sc_t = ex.scores_and_novelty(ex.fm.get(frozenset({f})), tsys)
    cP = {r["rule"]: r["c_P"] for r in results_c if r["fold"] == f}; cJ = {r["rule"]: r["c_J"] for r in results_c if r["fold"] == f}

    def eval_target(job):
        y, s, C, sysf, sid = job["y"], job["s"], job["C"], job["sysf"], job["sid"]; N = len(y); rows_ = []
        for rule in RULES:
            for r in range(N_EVAL_DRAWS):
                rng = np.random.default_rng([int(sid), n, r, STRATEGY_CODE[rule], cfg.get("seed", 0)])
                cal = select_calibrants(rng, y, n, rule, job["nsys"]); te = np.setdiff1d(np.arange(N), cal)
                mu, R, sl = loo_quantities(s[cal], y[cal], s[te], lambda st, scal, yc: sigma_fn(sig_f, C, sysf)(st, scal, yc, te, cal))
                muJ, RJ, _ = loo_quantities(s[cal], y[cal], s[te], None)
                for meth, (m_, R_, sl_, c_) in {"P": (mu, R, sl, 1.0), "P+SASC": (mu, R, sl, cP[rule]), "J+": (muJ, RJ, None, 1.0), "J+SASC": (muJ, RJ, None, cJ[rule])}.items():
                    lo, hi = windows(m_, R_, sl_, ALPHA, c_); cov = (lo <= y[te]) & (y[te] <= hi)
                    rows_.append({"system": sid, "fold": int(f), "rule": rule, "draw": r, "method": meth, "c": c_, "coverage": cov.mean(),
                                  "width_mean": (hi - lo).mean(), "width_median": np.median(hi - lo), "winkler": winkler(lo, hi, y[te], ALPHA).mean()})
        return rows_
    tj = [make_job(sid, sc_t, sig_f) for sid in tsys]
    for rr in Parallel(n_jobs=cfg["n_jobs"])(delayed(eval_target)(j) for j in tj):
        eval_rows += rr
    # NORMAN-mix curves on this fold's systems scored by the outer (out-of-source) model and σ_f
    nsys_f = [s for s in data.meta[(data.meta.fold == f) & (data.meta.n_compounds >= 40)].index if s != data.smrt_id]
    sc_n = ex.scores_and_novelty(ex.fm.get(frozenset({f})), nsys_f)
    nj = [j for j in (make_job(sid, sc_n, sig_f) for sid in nsys_f) if j["norman_pool"] is not None]
    for sid, o in Parallel(n_jobs=cfg["n_jobs"])(delayed(sim_system)(j) for j in nj):
        if "norman" in o:
            norman_curves_P.append(o["norman"][0]); norman_curves_J.append(o["norman"][1])
    print(f"fold {f}: targets evaluated; NORMAN-mix training systems so far: {len(norman_curves_P)}", flush=True)

cN = {"n_systems": len(norman_curves_P), "c_P": choose_c(np.array(norman_curves_P), ALPHA) if norman_curves_P else None,
      "c_J": choose_c(np.array(norman_curves_J), ALPHA) if norman_curves_J else None}
raw = os.path.join(cfg["out_dir"], "raw"); os.makedirs(raw, exist_ok=True)
pd.DataFrame(results_c).to_csv(os.path.join(raw, "sasc_c_values.csv"), index=False)
pd.DataFrame(eval_rows).to_parquet(os.path.join(raw, "sasc_rows.parquet"))
json.dump(cN, open(os.path.join(raw, "sasc_norman_c.json"), "w"), indent=1)
print("done", cN)

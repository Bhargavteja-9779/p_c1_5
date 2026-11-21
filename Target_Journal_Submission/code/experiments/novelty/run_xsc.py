"""Addendum G: cross-system conformal calibration (XSC) with laboratory-level (ε, δ) retention guarantees.

For every leave-source-out fold f:
  1. replay the calibrant protocol (rule ρ, number n) on every training system s of the other source groups g, with base
     scores from a model that excluded {f, g} and a difficulty model trained without group g (no system calibrates itself);
  2. choose the multiplier c per (ρ, n, window family) for an average retention of 1 − α, or for a laboratory-level
     guarantee P{own retention ≥ 1 − ε} ≥ 1 − δ;
  3. evaluate on the held-out target systems of fold f with the calibrant draws of the main experiment.
Usage: python experiments/novelty/run_xsc.py [config]
"""
import sys, os, json, time, math
import numpy as np, pandas as pd, yaml
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
os.chdir(os.path.join(HERE, "..", ".."))
from joblib import Parallel, delayed
from scipy import stats
from conformrt.evaluation.driver import Experiment, STRATEGY_CODE, _compound_frame
from conformrt.evaluation.system_eval import select_calibrants, _sigma_matrix, _sigma_cal_loo
from conformrt.evaluation.intervals import jackknife_plus, winkler
from conformrt.evaluation import sasc as S
from conformrt.evaluation.crosssys import (loo_residuals, scale_pivot, pivot_summary, choose_avg, choose_lab,
                                            choose_lab_exact, C_GRID, EPS_LIST)
from conformrt.models.difficulty import DifficultyModel, encode_system
from conformrt.models.mapping import MonotoneMap

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/full.yaml"))
quick = cfg["name"] == "quick"
N_LIST = [3, 5, 10, 20, 40]
JP_N = [10, 20, 40]                      # jackknife+-based families need n >= 1/α − 1
RULES = ["random", "span", "common", "polar"]
D_TRAIN = 2 if quick else 8
D_EVAL = {"random": 2 if quick else 20, "span": 2 if quick else 10, "common": 2 if quick else 10, "polar": 2 if quick else 10}
ALPHAS = [0.05, 0.10, 0.20]
DELTA = 0.10
SHRINK = 5.0  # pseudo-calibrants for the shrunk-scale family 'sigs' (fixed a priori)
OUT = os.path.join(cfg["out_dir"], "raw", "xsc"); os.makedirs(OUT, exist_ok=True)
ex = Experiment(dict(cfg, name=cfg["name"] + "_xsc")); data = ex.data; fs = data.fs
log = ex.log


def make_job(sid, scores_df, sig, excl):
    tab = data.rt[data.rt.system == sid].reset_index(drop=True); rows = fs.rows(tab.ik14)
    C = _compound_frame(data, rows, scores_df.loc[rows, "maxsim"].values, scores_df.loc[rows, "seen"].values)
    return {"sid": sid, "y": tab.rt.values.astype(float), "s": scores_df.loc[rows, "score"].values, "C": C,
            "sysf": encode_system(data.meta.loc[sid].to_dict()), "sig": sig,
            "nsys": data.fold_frequency(excl).reindex(tab.ik14.values).fillna(0).values, "logp": C["MolLogP"].values}


def draw_quantities(job, cal, te, want_P=False):
    y, s, C, sysf, sig = job["y"], job["s"], job["C"], job["sysf"], job["sig"]
    s_cal, y_cal, s_te, y_te = s[cal], y[cal], s[te], y[te]
    Cc, Ct = C.iloc[cal].reset_index(drop=True), C.iloc[te].reset_index(drop=True)
    loo_res, loo = loo_residuals(s_cal, y_cal)
    q = {"loo_res": loo_res, "loo": loo}
    q["flat"] = scale_pivot(s_cal, y_cal, s_te, y_te, None, None, loo_res)
    sig_te = _sigma_matrix(sig, Ct, s_te, s_cal, y_cal, sysf, loo=False)
    sig_cal_loo = _sigma_cal_loo(sig, Cc, s_cal, y_cal, sysf)
    q["sig"] = scale_pivot(s_cal, y_cal, s_te, y_te, sig_te, sig_cal_loo, loo_res)
    q["sigs"] = scale_pivot(s_cal, y_cal, s_te, y_te, sig_te, sig_cal_loo, loo_res, shrink=SHRINK)
    q["sig_te"] = sig_te
    if len(cal) in JP_N:
        q["mu_loo_te"] = np.vstack([m.predict(s_te) for m in loo])
        if want_P:
            q["P_sl_te"] = _sigma_matrix(sig, Ct, s_te, s_cal, y_cal, sysf, loo=True)
            q["P_R"] = np.abs(loo_res) / sig_cal_loo
    return q


def sim_system(job):
    """Replay the protocol on one training system; returns per (rule, n, family) curves and required multipliers."""
    y, sid = job["y"], job["sid"]; N = len(y); out = {}
    for rule in RULES:
        for n in N_LIST:
            if N < n + 20:
                continue
            acc = {}
            for r in range(D_TRAIN):
                rng = np.random.default_rng([int(sid), n, r, STRATEGY_CODE[rule], 99])
                cal = select_calibrants(rng, y, n, rule, job["nsys"], None, job["logp"])
                te = np.setdiff1d(np.arange(N), cal)
                q = draw_quantities(job, cal, te)
                for fam in ("flat", "sig", "sigs"):
                    curve, req = pivot_summary(q[fam][2])
                    acc.setdefault(fam, []).append((curve, req))
                if n in JP_N:  # jackknife+ with a multiplier on the scores (flat); coverage curve on the SASC grid
                    cur = S.coverage_curve(q["mu_loo_te"], np.abs(q["loo_res"]), None, y[te], 0.10)
                    req = {e: (float(S.GRID[np.flatnonzero(cur >= 1 - e)[0]]) if (cur >= 1 - e).any() else np.inf) for e in EPS_LIST}
                    acc.setdefault("jp", []).append((cur.astype(np.float32), req))
            for fam, lst in acc.items():
                out[(rule, n, fam)] = (np.mean([c for c, _ in lst], 0), {e: np.array([q_[e] for _, q_ in lst]) for e in EPS_LIST})
    return sid, out


def calibrate(sim_outs):
    """Multipliers per (rule, n, family): avg(α) and lab(ε, δ) (draw-averaged) and lab_exact (first draw)."""
    keys = set(k for _, o in sim_outs for k in o)
    table = []
    for (rule, n, fam) in sorted(keys):
        per = [o[(rule, n, fam)] for _, o in sim_outs if (rule, n, fam) in o]
        curves = np.array([p[0] for p in per]); grid = S.GRID if fam == "jp" else C_GRID
        row = {"rule": rule, "n": n, "family": fam, "S": len(per)}
        for a in ALPHAS:
            row[f"c_avg_{a}"] = choose_avg(curves, a, grid)
        for e in EPS_LIST:
            row[f"c_lab_{e}"] = choose_lab([p[1][e] for p in per], e, DELTA, grid)
            row[f"c_labx_{e}"] = choose_lab_exact(np.array([p[1][e][0] for p in per]), DELTA)
        table.append(row)
    return pd.DataFrame(table)


def eval_target(job, ctab):
    y, s, sid = job["y"], job["s"], job["sid"]; N = len(y); rows_ = []
    cdict = {(r["rule"], r["n"], r["family"]): r for r in ctab.to_dict("records")}
    for rule in RULES:
        for n in N_LIST:
            if N < n + cfg["min_test"]:
                continue
            for r in range(D_EVAL[rule]):
                rng = np.random.default_rng([int(sid), n, r, STRATEGY_CODE[rule], cfg.get("seed", 0)])
                cal = select_calibrants(rng, y, n, rule, job["nsys"], None, job["logp"])
                te = np.setdiff1d(np.arange(N), cal); y_te = y[te]
                q = draw_quantities(job, cal, te, want_P=True)
                levels = ALPHAS if rule == "random" else [0.10]
                cands = {}
                for a in levels:
                    mu_f, tau_f, _ = q["flat"]; mu_s, tau_s, _ = q["sig"]
                    z = stats.norm.ppf(1 - a / 2); rmse = np.sqrt(np.mean(q["loo_res"] ** 2))
                    cands[("B1_gauss_rmse", a)] = (mu_f - z * rmse, mu_f + z * rmse)
                    if n in JP_N:
                        cands[("B5_jackknife_plus", a)] = jackknife_plus(q["mu_loo_te"], np.abs(q["loo_res"]), None, a)
                        cands[("P", a)] = jackknife_plus(q["mu_loo_te"], q["P_R"], q["P_sl_te"], a)
                    for fam, (mu, tau, _) in (("flat", q["flat"]), ("sig", q["sig"]), ("sigs", q["sigs"])):
                        c = cdict.get((rule, n, fam))
                        if c is None:
                            continue
                        sc = tau * (1.0 if fam == "flat" else q["sig_te"])
                        ca = c[f"c_avg_{a}"]; cl = c[f"c_lab_{a}"]
                        cands[(f"XSC_{fam}_avg", a)] = (mu - ca * sc, mu + ca * sc)
                        cands[(f"XSC_{fam}_lab", a)] = (mu - cl * sc, mu + cl * sc)
                    c = cdict.get((rule, n, "jp"))
                    if c is not None and n in JP_N:
                        cl = c[f"c_lab_{a}"]
                        cands[("JP_lab", a)] = S.windows(q["mu_loo_te"], np.abs(q["loo_res"]), None, 0.10, cl)
                for (meth, a), (lo, hi) in cands.items():
                    cov = (lo <= y_te) & (y_te <= hi); w = hi - lo; fin = bool(np.isfinite(w).all())
                    rows_.append({"system": sid, "rule": rule, "n": n, "draw": r, "method": meth, "level": a,
                                  "coverage": float(cov.mean()), "width_median": float(np.median(w)) if fin else np.inf,
                                  "width_mean": float(w.mean()) if fin else np.inf,
                                  "winkler": float(winkler(lo, hi, y_te, a).mean()) if fin else np.inf,
                                  "span": float(y.max() - y.min()), "m_test": len(te)})
    return rows_


if __name__ == "__main__":
    targets = data.targets(cfg["min_target_compounds"])
    folds = sorted(targets.fold.unique())
    for f in folds:
        done = os.path.join(OUT, f"fold{f}_DONE")
        if os.path.exists(done):
            log.info(f"[xsc] fold {f}: done (checkpoint)"); continue
        t0 = time.time()
        rows_f = ex.difficulty_training_rows(f, "order")
        rows_f["sfold"] = data.meta.loc[rows_f.system.astype(str), "fold"].values
        sim_path = os.path.join(OUT, f"fold{f}_sim.pkl")
        if os.path.exists(sim_path):
            sim_outs = pd.read_pickle(sim_path)
        else:
            jobs = []
            for g in sorted(set(data.meta.fold.unique()) - {f}):
                sub = rows_f[rows_f.sfold != g]
                sig_fg = DifficultyModel(n_jobs=cfg["n_jobs"], **cfg["difficulty_model"]).fit(sub, sub["target"].values)
                sig_fg.model.set_params(n_jobs=1)
                systems = [s_ for s_ in data.meta[(data.meta.fold == g) & (data.meta.n_compounds >= 40)].index if s_ != data.smrt_id]
                if quick:
                    systems = systems[:4]
                sc = ex.scores_and_novelty(ex.fm.get(frozenset({f, g})), systems)
                jobs += [make_job(sid, sc, sig_fg, {f, g}) for sid in systems]
            sim_outs = Parallel(n_jobs=cfg["n_jobs"], backend="loky")(delayed(sim_system)(j) for j in jobs)
            pd.to_pickle(sim_outs, sim_path)
        ctab = calibrate(sim_outs); ctab["fold"] = f
        ctab.to_csv(os.path.join(OUT, f"fold{f}_multipliers.csv"), index=False)
        log.info(f"[xsc] fold {f}: protocol replayed on {len(sim_outs)} training systems in {time.time()-t0:.0f}s")
        sig_f = ex.train_sigma_models(f, "order", which=["P"])[0]["P"]
        tsys = list(targets[targets.fold == f].index)[: (3 if quick else None)]
        sc_t = ex.scores_and_novelty(ex.fm.get(frozenset({f})), tsys)
        tj = [make_job(sid, sc_t, sig_f, {f}) for sid in tsys]
        res = Parallel(n_jobs=cfg["n_jobs"], backend="loky")(delayed(eval_target)(j, ctab) for j in tj)
        df = pd.DataFrame([r for rr in res for r in rr]); df["fold"] = f
        df.to_parquet(os.path.join(OUT, f"fold{f}_eval.parquet"))
        open(done, "w").write("ok")
        log.info(f"[xsc] fold {f}: evaluated {len(tsys)} targets; total {time.time()-t0:.0f}s")
    print("XSC_DONE", flush=True)

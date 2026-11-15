"""Addendum B (post hoc, added after reading the ACA journal benchmark, before inspecting main results):
fixed-tolerance RT windows used in screening practice, evaluated on the IDENTICAL calibrant draws as the main run
(same seeds; same cached fold base models), paired by (system, n, draw).

Windows (centred on the monotone calibrant map g(s(x))): ±0.25, ±0.5, ±1.0, ±2.0 min; ±5 %, ±10 % of predicted RT;
±2·RMSE_LOO and ±3·RMSE_LOO. Metrics: coverage, width, Winkler score at α = 0.10 (for comparability with Table 2), and
isomeric-decoy elimination at n = 20 (draws 0–4).
Usage: python experiments/run_extra_baselines.py [config]"""
import sys, os, yaml, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
from joblib import Parallel, delayed
from conformrt.evaluation.driver import Experiment, STRATEGY_CODE
from conformrt.evaluation.system_eval import select_calibrants
from conformrt.evaluation.intervals import winkler
from conformrt.models.mapping import MonotoneMap

WINDOWS = [("FIX_0.25min", "abs", 0.25), ("FIX_0.5min", "abs", 0.5), ("FIX_1min", "abs", 1.0), ("FIX_2min", "abs", 2.0),
           ("FIX_5pct", "rel", 0.05), ("FIX_10pct", "rel", 0.10), ("RMSE_x2", "rmse", 2.0), ("RMSE_x3", "rmse", 3.0)]


def eval_sys(job):
    sid, y, s, cfg, dec = job["sid"], job["y"], job["s"], job["cfg"], job.get("dec")
    N = len(y); rows, drows = [], []
    for n in cfg["n_list"]:
        if N < n + cfg["min_test"]:
            continue
        for r in range(cfg["draws"]):
            rng = np.random.default_rng([int(sid), n, r, STRATEGY_CODE["random"], cfg.get("seed", 0)])
            cal = select_calibrants(rng, y, n, "random")
            te = np.setdiff1d(np.arange(N), cal)
            full = MonotoneMap().fit(s[cal], y[cal]); mu = full.predict(s[te])
            loo_res = np.array([y[cal][i] - MonotoneMap().fit(np.delete(s[cal], i), np.delete(y[cal], i)).predict(s[cal][i:i + 1])[0] for i in range(n)])
            rmse = np.sqrt(np.mean(loo_res ** 2))
            def hw(kind, v, m):
                return np.full_like(m, v) if kind == "abs" else (np.abs(m) * v if kind == "rel" else np.full_like(m, v * rmse))
            for name, kind, v in WINDOWS:
                h = hw(kind, v, mu); lo, hi = mu - h, mu + h; cov = (lo <= y[te]) & (y[te] <= hi)
                rows.append({"system": sid, "n": n, "draw": r, "method": name, "alpha": 0.1, "m_test": len(te), "coverage": cov.mean(),
                             "width_mean": (hi - lo).mean(), "width_median": np.median(hi - lo),
                             "width_rel_mean": (hi - lo).mean() / (y.max() - y.min()), "winkler": winkler(lo, hi, y[te], 0.1).mean()})
            if n == cfg["primary_n"] and r < cfg["decoy_draws"] and dec is not None:
                par, sd = dec; keep = np.isin(par, te)
                if keep.sum():
                    mu_d = full.predict(sd[keep]); yp = y[par[keep]]
                    for name, kind, v in WINDOWS:
                        h = hw(kind, v, mu_d); elim = ~((mu_d - h <= yp) & (yp <= mu_d + h))
                        e = pd.DataFrame({"p": par[keep], "e": elim}).groupby("p").e.mean()
                        hm = hw(kind, v, mu); cov_t = pd.Series((mu - hm <= y[te]) & (y[te] <= mu + hm), index=te)
                        drows.append({"system": sid, "draw": r, "method": name, "alpha": 0.1, "n_parents": len(e),
                                      "decoy_elim_rate": float(e.mean()), "true_retained": float(cov_t.loc[e.index].mean())})
    return rows, drows


if __name__ == "__main__":
    cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/full.yaml"))
    ex = Experiment(cfg); data = ex.data; targets = data.targets(cfg["min_target_compounds"])
    allr, alld = [], []
    for f in sorted(targets.fold.unique()):
        model, train_ik = ex.fm.get(frozenset({f}))
        tsys = list(targets[targets.fold == f].index)[: cfg.get("systems_limit") or None]
        jobs = []
        for sid in tsys:
            tab = data.rt[data.rt.system == sid].reset_index(drop=True); rows = data.fs.rows(tab.ik14)
            par, drows = ex.decoys_for(tab.ik14.values, int(sid))
            jobs.append({"sid": sid, "y": tab.rt.values.astype(float), "s": model.predict(data.fs, rows), "cfg": cfg,
                         "dec": (par, model.predict(data.fs, drows)) if len(drows) else None})
        out = Parallel(n_jobs=cfg["n_jobs"])(delayed(eval_sys)(j) for j in jobs)
        for r_, d_ in out:
            allr += [dict(x, fold=f) for x in r_]; alld += [dict(x, fold=f) for x in d_]
        print(f"fold {f}: {len(tsys)} systems", flush=True)
    raw = os.path.join(cfg["out_dir"], "raw"); os.makedirs(raw, exist_ok=True)
    pd.DataFrame(allr).to_parquet(os.path.join(raw, "fixedwin_rows.parquet")); pd.DataFrame(alld).to_parquet(os.path.join(raw, "fixedwin_decoys.parquet"))
    print("done", len(allr), len(alld))

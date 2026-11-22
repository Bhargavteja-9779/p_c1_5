"""Addendum G4: per-draw pivot summaries on target systems, so that any multiplier (ε/δ sweeps, metadata-weighted XSC,
exact one-draw calibration) can be evaluated offline. Same calibrant draws, models and σ as run_xsc.py.
Stores per (system, rule, n, draw, family): τ̂, median σ(x) of the test compounds, and 201 quantiles of the pivot u(x).
Usage: python experiments/novelty/run_xsc_pivots.py [config]"""
import sys, os, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_xsc as X  # noqa: E402
from joblib import Parallel, delayed
from conformrt.evaluation.system_eval import select_calibrants

cfg, data, ex = X.cfg, X.data, X.ex
QL = np.linspace(0, 1, 201)


def pivots(job):
    y, sid = job["y"], job["sid"]; N = len(y); rows = []
    for rule in X.RULES:
        for n in X.N_LIST:
            if N < n + cfg["min_test"]:
                continue
            for r in range(X.D_EVAL[rule]):
                rng = np.random.default_rng([int(sid), n, r, X.STRATEGY_CODE[rule], cfg.get("seed", 0)])
                cal = select_calibrants(rng, y, n, rule, job["nsys"], None, job["logp"])
                te = np.setdiff1d(np.arange(N), cal)
                q = X.draw_quantities(job, cal, te)
                for fam in ("flat", "sigs"):
                    mu, tau, u = q[fam]
                    sm = 1.0 if fam == "flat" else float(np.median(q["sig_te"]))
                    rows.append({"system": sid, "rule": rule, "n": n, "draw": r, "family": fam, "tau": tau, "sig_med": sm,
                                 "m_test": len(te), "span": float(y.max() - y.min()), "uq": np.quantile(u, QL).astype(np.float32)})
    return rows


if __name__ == "__main__":
    targets = data.targets(cfg["min_target_compounds"])
    for f in sorted(targets.fold.unique()):
        pth = os.path.join(X.OUT, f"fold{f}_pivots.parquet")
        if os.path.exists(pth):
            continue
        t0 = time.time()
        sig_f = ex.train_sigma_models(f, "order", which=["P"])[0]["P"]
        tsys = list(targets[targets.fold == f].index)[: (3 if X.quick else None)]
        sc_t = ex.scores_and_novelty(ex.fm.get(frozenset({f})), tsys)
        res = Parallel(n_jobs=cfg["n_jobs"], backend="loky")(delayed(pivots)(X.make_job(s, sc_t, sig_f, {f})) for s in tsys)
        df = pd.DataFrame([r for rr in res for r in rr]); df["fold"] = f
        df["uq"] = df["uq"].apply(list)
        df.to_parquet(pth)
        X.log.info(f"[pivots] fold {f}: {len(tsys)} systems in {time.time()-t0:.0f}s")
    print("PIVOTS_DONE", flush=True)

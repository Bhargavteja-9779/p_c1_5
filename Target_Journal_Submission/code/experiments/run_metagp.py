"""Addendum E (post hoc; prompted by the citation audit, 2026-11-16): comparator B7 = Gaussian-process projection with
META-LEARNED hyper-parameters, in the spirit of the meta-learned GP projections of Garcia et al. (J. Cheminform. 2022).

For each leave-source-out fold f, the kernel and noise hyper-parameters of a GP that maps the (calibrant-standardised)
base score to the (calibrant-standardised) retention time are meta-learned on the TRAINING systems only: each training
source group g is scored by the base model trained without f and g (the inner cross-fitting used for the difficulty model),
random calibrant sets of size 10, 20 and 40 are drawn (3 draws, same seeds as the difficulty-model simulation), and the
hyper-parameters maximise the mean Gaussian predictive log-likelihood of the held-out compounds (each system weighted
equally). The fixed hyper-parameters are then applied to every target system of fold f, scored by the outer model, on
exactly the calibrant draws of the main study (same seeds); nothing is tuned on target data. As a check that the draws are
identical, jackknife+ coverage is recomputed for n = 20, draw 0 and compared with the main run.

Outputs: <out_dir>/raw/metagp_rows.parquet (same schema as main_rows), metagp_theta.csv, metagp_check.csv
Usage: python experiments/run_metagp.py [config]"""
import os, sys, time, json
import numpy as np, pandas as pd, yaml
from scipy import optimize, stats
from scipy.linalg import cho_solve, solve_triangular
from joblib import Parallel, delayed
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
from conformrt.evaluation.driver import Experiment, STRATEGY_CODE  # noqa: E402
from conformrt.evaluation.system_eval import select_calibrants  # noqa: E402
from conformrt.evaluation.intervals import jackknife_plus, winkler  # noqa: E402
from conformrt.models.mapping import MonotoneMap  # noqa: E402

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/full.yaml"))
SEED = cfg.get("seed", 0); NJ = cfg["n_jobs"]
MAX_TEST_META = 400  # held-out compounds per simulated draw used in the meta-objective (random subset, fixed seed)
NAMES = ["vc", "vl", "vr", "ell", "vn"]  # constant, linear, RBF variances; RBF length scale; noise variance (standardised units)
BOUNDS = [(-9, 3), (-9, 3), (-9, 3), (np.log(0.03), np.log(30.0)), (-9, 2)]


def unpack(x):
    return dict(zip(NAMES, np.exp(np.asarray(x, float))))


def kern(a, b, th):
    d = a[:, None] - b[None, :]
    return th["vc"] + th["vl"] * a[:, None] * b[None, :] + th["vr"] * np.exp(-0.5 * d * d / th["ell"] ** 2)


def gp_predict(s_cal, y_cal, s_q, th):
    """GP posterior predictive (mean, sd) in minutes; inputs and outputs standardised by the calibrants."""
    ms, ss = s_cal.mean(), s_cal.std() + 1e-9; my, sy = y_cal.mean(), y_cal.std() + 1e-9
    a = (s_cal - ms) / ss; b = (s_q - ms) / ss; yc = (y_cal - my) / sy
    K = kern(a, a, th) + (th["vn"] + 1e-8) * np.eye(len(a))
    L = np.linalg.cholesky(K)
    w = cho_solve((L, True), yc)
    Ks = kern(a, b, th)
    mu = Ks.T @ w
    v = solve_triangular(L, Ks, lower=True)
    var = np.maximum(th["vc"] + th["vl"] * b * b + th["vr"] - np.sum(v * v, axis=0), 1e-12) + th["vn"]
    return mu * sy + my, np.sqrt(var) * sy


def meta_tasks(train_sets):
    """Simulated calibrant draws on training systems (same scheme and seeds as simulate_calibration_rows)."""
    tasks = []
    for sid, s, y in train_sets:
        N = len(y)
        if N < 30:
            continue
        for n in cfg["sim_n_list"]:
            if N < n + 10:
                continue
            for r in range(cfg["sim_draws"]):
                rng = np.random.default_rng([int(sid), n, r, SEED])
                cal = rng.choice(N, n, replace=False); te = np.setdiff1d(np.arange(N), cal)
                if len(te) > MAX_TEST_META:
                    te = np.random.default_rng([int(sid), n, r, SEED, 1]).choice(te, MAX_TEST_META, replace=False)
                tasks.append((sid, s[cal], y[cal].astype(float), s[te], y[te].astype(float)))
    return tasks


def neg_mean_loglik(x, tasks):
    th = unpack(x); per_sys = {}
    for sid, sc, yc, st, yt in tasks:
        try:
            mu, sd = gp_predict(sc, yc, st, th)
        except np.linalg.LinAlgError:
            return 1e6
        per_sys.setdefault(sid, []).append(np.mean(stats.norm.logpdf(yt, mu, sd)))
    return -float(np.mean([np.mean(v) for v in per_sys.values()]))


def meta_learn(f, train_sets):
    tasks = meta_tasks(train_sets); t = time.time()
    x0 = np.log([0.05, 1.0, 0.1, 1.0, 0.05])
    best = None
    for start in [x0, np.log([0.5, 0.5, 0.5, 0.3, 0.02])]:  # two starts, keep the better optimum
        r = optimize.minimize(neg_mean_loglik, start, args=(tasks,), method="L-BFGS-B", bounds=BOUNDS, options={"maxiter": 200})
        if best is None or r.fun < best.fun:
            best = r
    th = unpack(best.x)
    return {"fold": int(f), "n_train_systems": len({t_[0] for t_ in tasks}), "n_tasks": len(tasks), "neg_mean_loglik": float(best.fun),
            "converged": bool(best.success), "seconds": round(time.time() - t, 1), **{k: float(v) for k, v in th.items()}}


def evaluate_target(sid, f, s, y, th):
    out, chk = [], []
    N = len(y); span = float(y.max() - y.min())
    for n in cfg["n_list"]:
        if N < n + cfg["min_test"]:
            continue
        for r in range(cfg["draws"]):
            rng = np.random.default_rng([int(sid), n, r, STRATEGY_CODE["random"], SEED])  # identical to the main run
            cal = select_calibrants(rng, y, n, "random")
            te = np.setdiff1d(np.arange(N), cal)
            mu, sd = gp_predict(s[cal], y[cal].astype(float), s[te], th)
            for a in cfg["alphas"]:
                z = stats.norm.ppf(1 - a / 2); lo, hi = mu - z * sd, mu + z * sd; w = hi - lo
                out.append({"system": sid, "n": n, "draw": r, "method": "B7_gp_meta", "alpha": a, "m_test": len(te),
                            "coverage": float(np.mean((lo <= y[te]) & (y[te] <= hi))), "width_mean": float(w.mean()),
                            "width_median": float(np.median(w)), "width_rel_mean": float(w.mean() / span),
                            "winkler": float(winkler(lo, hi, y[te], a).mean()), "fold": int(f)})
            if n == cfg["primary_n"] and r == 0:  # consistency check against the main run (same draws -> same jackknife+ coverage)
                sc, yc = s[cal], y[cal].astype(float)
                loo = [MonotoneMap().fit(np.delete(sc, i), np.delete(yc, i)) for i in range(n)]
                res = np.abs(yc - np.array([loo[i].predict(sc[i:i + 1])[0] for i in range(n)]))
                lo, hi = jackknife_plus(np.vstack([m.predict(s[te]) for m in loo]), res, None, 0.1)
                chk.append({"system": sid, "jplus_cov_recomputed": float(np.mean((lo <= y[te]) & (y[te] <= hi)))})
    return out, chk


if __name__ == "__main__":
    t0 = time.time()
    ex = Experiment(dict(cfg, name=cfg["name"] + "_metagp")); data = ex.data
    targets = data.targets(cfg["min_target_compounds"])
    meta_in, tgt_in = {}, []
    for f in sorted(targets.fold.unique()):
        train_sets = []
        for g in sorted(set(data.meta.fold.unique()) - {f}):
            model, _ = ex.fm.get(frozenset({f, g}))
            for sid in [s_ for s_ in data.meta[data.meta.fold == g].index if s_ != data.smrt_id]:
                tab = data.rt[data.rt.system == sid]
                if len(tab) < 30:
                    continue
                train_sets.append((sid, model.predict(data.fs, data.fs.rows(tab.ik14)), tab.rt.values))
        meta_in[f] = train_sets
        outer, _ = ex.fm.get(frozenset({f}))
        for sid in targets[targets.fold == f].index:
            tab = data.rt[data.rt.system == sid].reset_index(drop=True)
            tgt_in.append((sid, f, outer.predict(data.fs, data.fs.rows(tab.ik14)), tab.rt.values.astype(float)))
        print(f"fold {f}: {len(train_sets)} training systems scored", flush=True)
    thetas = Parallel(n_jobs=NJ)(delayed(meta_learn)(f, ts) for f, ts in meta_in.items())
    theta = pd.DataFrame(thetas); print(theta.to_string(), flush=True)
    TH = {int(r.fold): {k: r[k] for k in NAMES} for _, r in theta.iterrows()}
    outs = Parallel(n_jobs=NJ)(delayed(evaluate_target)(sid, f, s, y, TH[int(f)]) for sid, f, s, y in tgt_in)
    rows = pd.DataFrame([r for o, _ in outs for r in o]); chk = pd.DataFrame([x for _, cl in outs for x in cl])
    rows["base"] = "order"; rows["strategy"] = "random"; rows["perturb"] = "none"; rows["experiment"] = "metagp"
    raw = os.path.join(cfg["out_dir"], "raw"); os.makedirs(raw, exist_ok=True)
    main = pd.read_parquet(os.path.join(raw, "main_rows.parquet"))
    m20 = main[(main.method == "B5_jackknife_plus") & (main.n == cfg["primary_n"]) & (main.draw == 0) & np.isclose(main.alpha, 0.1)].set_index("system").coverage
    chk["jplus_cov_main"] = chk.system.map(m20); chk["abs_diff"] = (chk.jplus_cov_recomputed - chk.jplus_cov_main).abs()
    rows.to_parquet(os.path.join(raw, "metagp_rows.parquet")); theta.to_csv(os.path.join(raw, "metagp_theta.csv"), index=False)
    chk.to_csv(os.path.join(raw, "metagp_check.csv"), index=False)
    print(f"draw check: max |coverage difference| vs main run = {chk.abs_diff.max():.2e} over {len(chk)} systems", flush=True)
    s20 = rows[(rows.n == 20) & np.isclose(rows.alpha, 0.1)].groupby("system").coverage.mean()
    print(f"B7 n=20 a=0.10: mean coverage {s20.mean():.4f} over {len(s20)} systems; done in {time.time() - t0:.0f} s", flush=True)

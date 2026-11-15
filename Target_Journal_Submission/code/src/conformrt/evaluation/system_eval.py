"""Evaluation of all interval methods on one target chromatographic system (run in a worker process)."""
from __future__ import annotations
import time
import warnings
import numpy as np
import pandas as pd

from ..models.mapping import MonotoneMap, PSplineMap, GPMap, LinearMap
from ..models.difficulty import geometry, geometry_loo
from .intervals import jackknife_plus, split_conformal, gaussian_rmse, winkler, k_upper

warnings.filterwarnings("ignore")
SIM_BINS = [(-1, 0.4, "sim<0.4"), (0.4, 0.7, "sim0.4-0.7"), (0.7, 1.01, "sim>=0.7")]


def _sigma_matrix(model, Cq: pd.DataFrame, s_q, s_cal, y_cal, sysf, loo: bool):
    """σ predictions. loo=True → (n, m) array of σ_{-i}(x_j); else (m,) with full-calibrant geometry."""
    cols = model.cols
    m = len(s_q)
    if loo:
        G = geometry_loo(s_q, s_cal, y_cal); n = len(s_cal)
        X = np.empty((n * m, len(cols)), dtype=np.float32)
        for c, name in enumerate(cols):
            if name in ("pos", "extrap", "gap"):
                X[:, c] = G[name].ravel()
            elif name == "cal_span":
                X[:, c] = np.repeat(G["cal_span"], m)
            elif name == "n_cal":
                X[:, c] = G["n_cal"]
            elif name in sysf:
                X[:, c] = sysf[name]
            else:
                X[:, c] = np.tile(Cq[name].values, n)
        return model.sigma(X).reshape(n, m)
    G = geometry(s_q, s_cal, y_cal)
    X = np.empty((m, len(cols)), dtype=np.float32)
    for c, name in enumerate(cols):
        if name in ("pos", "extrap", "gap"):
            X[:, c] = G[name]
        elif name in ("cal_span", "n_cal"):
            X[:, c] = G[name]
        elif name in sysf:
            X[:, c] = sysf[name]
        else:
            X[:, c] = Cq[name].values
    return model.sigma(X)


def _sigma_cal_loo(model, Cc: pd.DataFrame, s_cal, y_cal, sysf):
    """σ_{-i}(x_i) for each calibrant i (geometry of the LOO calibrant set evaluated at calibrant i)."""
    S = _sigma_matrix(model, Cc, s_cal, s_cal, y_cal, sysf, loo=True)
    return np.diag(S).copy()


def select_calibrants(rng, y, n, strategy, n_systems_global=None):
    N = len(y)
    if strategy == "random":
        return rng.choice(N, n, replace=False)
    if strategy == "span":
        order = np.argsort(y)
        return np.array([rng.choice(b) for b in np.array_split(order, n)])
    if strategy == "common":
        key = np.lexsort((rng.random(N), -np.asarray(n_systems_global, float)))
        return key[:n]
    raise ValueError(strategy)


def evaluate_system(job: dict) -> dict:
    sid = job["sid"]; y = job["y"]; s = job["s"]; C = job["C"]; sysf = job["sysf"]; cfg = job["cfg"]
    sig_models = job["sig_models"]; ens = job.get("ens"); strata_lab = job["strata"]
    N = len(y); span_sys = float(y.max() - y.min())
    out_rows, strata_rows, decoy_rows, err_rows, timing = [], [], [], [], []
    alphas = cfg["alphas"]; primary_n = cfg["primary_n"]
    for n in cfg["n_list"]:
        if N < n + cfg["min_test"]:
            continue
        for r in range(cfg["draws"]):
            rng = np.random.default_rng([int(sid), n, r, cfg["strategy_code"], cfg.get("seed", 0)])
            cal = select_calibrants(rng, y, n, cfg["strategy"], job.get("n_systems_global"))
            te = np.setdiff1d(np.arange(N), cal)
            s_cal, s_te, y_te = s[cal], s[te], y[te]
            y_cal = y[cal].astype(float).copy()
            if cfg.get("perturb") == "noise":
                y_cal = y_cal + rng.normal(0, 0.02 * span_sys, size=n)
            elif cfg.get("perturb") == "outlier":
                y_cal[rng.integers(n)] += 0.10 * span_sys
            Cc, Ct = C.iloc[cal].reset_index(drop=True), C.iloc[te].reset_index(drop=True)
            is_primary = (n == primary_n)
            t0 = time.perf_counter()
            # ---- monotone LOO fits (shared by B1, B5, B6, P*) ----
            loo = [MonotoneMap().fit(np.delete(s_cal, i), np.delete(y_cal, i)) for i in range(n)]
            full = MonotoneMap().fit(s_cal, y_cal)
            mu_loo_cal = np.array([loo[i].predict(s_cal[i:i + 1])[0] for i in range(n)])
            mu_loo_te = np.vstack([m.predict(s_te) for m in loo])
            mu_full_te = full.predict(s_te)
            loo_res = y_cal - mu_loo_cal
            t_map = time.perf_counter() - t0
            intervals = {}  # method -> function(alpha) -> (lo, hi)
            methods = cfg["methods"]
            if "B1_gauss_rmse" in methods:
                intervals["B1_gauss_rmse"] = lambda a: gaussian_rmse(mu_full_te, loo_res, a)
            if "B2_pspline_PI" in methods and (n in cfg.get("expensive_n", cfg["n_list"])):
                ps = PSplineMap().fit(s_cal, y_cal)
                intervals["B2_pspline_PI"] = lambda a, ps=ps: ps.predict_interval(s_te, a)
            if "B3_gp_PI" in methods and (n in cfg.get("expensive_n", cfg["n_list"])):
                try:
                    gp = GPMap().fit(s_cal, y_cal)
                    intervals["B3_gp_PI"] = lambda a, gp=gp: gp.predict_interval(s_te, a)
                except Exception:
                    pass
            if "B4_split_conformal" in methods:
                perm = rng.permutation(n); A, B = perm[: (n + 1) // 2], perm[(n + 1) // 2:]
                mA = MonotoneMap().fit(s_cal[A], y_cal[A]); rB = y_cal[B] - mA.predict(s_cal[B]); muA = mA.predict(s_te)
                intervals["B4_split_conformal"] = lambda a, muA=muA, rB=rB: split_conformal(muA, rB, a)
            if "B5_jackknife_plus" in methods:
                intervals["B5_jackknife_plus"] = lambda a: jackknife_plus(mu_loo_te, np.abs(loo_res), None, a)
            if "B6_jk+_ensemble" in methods and ens is not None and (n in cfg.get("expensive_n", cfg["n_list"])):
                # LOO-consistent ensemble spread: replicate maps refitted without calibrant i (sigma_{-i})
                E_ = ens.shape[0]; sd_te = np.empty((n, len(te))); sd_cal = np.empty(n)
                for i in range(n):
                    keep_i = np.arange(n) != i
                    pt, pc = [], []
                    for e in range(E_):
                        me = MonotoneMap().fit(ens[e, cal][keep_i], y_cal[keep_i])
                        pt.append(me.predict(ens[e, te])); pc.append(me.predict(ens[e, cal[i:i + 1]])[0])
                    sd_te[i] = np.std(pt, axis=0) + 0.01; sd_cal[i] = np.std(pc) + 0.01
                Rens = np.abs(loo_res) / sd_cal
                intervals["B6_jk+_ensemble"] = lambda a, Rens=Rens, sd_te=sd_te: jackknife_plus(mu_loo_te, Rens, sd_te, a)
            t1 = time.perf_counter()
            sig_cache = {}
            for name, model in sig_models.items():
                if name != "P" and not is_primary:
                    continue
                sl_te = _sigma_matrix(model, Ct, s_te, s_cal, y_cal, sysf, loo=True)
                sl_cal = _sigma_cal_loo(model, Cc, s_cal, y_cal, sysf)
                sig_cache[name] = (sl_te, sl_cal)
                R = np.abs(loo_res) / sl_cal
                intervals[name] = (lambda a, sl_te=sl_te, R=R: jackknife_plus(mu_loo_te, R, sl_te, a))
            t_sig = time.perf_counter() - t1
            if is_primary and "P" in sig_models:
                model = sig_models["P"]
                sf_te = _sigma_matrix(model, Ct, s_te, s_cal, y_cal, sysf, loo=False)
                sf_cal = _sigma_matrix(model, Cc, s_cal, s_cal, y_cal, sysf, loo=False)
                Rn = np.abs(y_cal - full.predict(s_cal)) / sf_cal
                def naive(a, Rn=Rn, sf_te=sf_te):
                    k = k_upper(n, a); q = np.inf if k > n else np.sort(Rn)[k - 1]
                    return mu_full_te - q * sf_te, mu_full_te + q * sf_te
                intervals["P_naive_noLOO"] = naive
                lin = [LinearMap().fit(np.delete(s_cal, i), np.delete(y_cal, i)) for i in range(n)]
                mu_lin_cal = np.array([lin[i].predict(s_cal[i:i + 1])[0] for i in range(n)])
                mu_lin_te = np.vstack([m.predict(s_te) for m in lin])
                sl_te, sl_cal = sig_cache["P"]
                Rl = np.abs(y_cal - mu_lin_cal) / sl_cal
                intervals["P_linear_map"] = lambda a, Rl=Rl, mu_lin_te=mu_lin_te, sl_te=sl_te: jackknife_plus(mu_lin_te, Rl, sl_te, a)
            timing.append({"system": sid, "n": n, "draw": r, "t_map_loo_s": t_map, "t_sigma_s": t_sig, "m_test": len(te)})
            # ---- metrics ----
            for meth, f in intervals.items():
                for a in alphas:
                    lo, hi = f(a)
                    cov = (lo <= y_te) & (y_te <= hi)
                    w = hi - lo
                    finite = np.isfinite(w).all()
                    out_rows.append({"system": sid, "n": n, "draw": r, "method": meth, "alpha": a, "m_test": len(te),
                                     "coverage": cov.mean(), "width_mean": w.mean() if finite else np.inf,
                                     "width_median": np.median(w) if finite else np.inf,
                                     "width_rel_mean": (w.mean() / span_sys) if finite else np.inf,
                                     "winkler": winkler(lo, hi, y_te, a).mean() if finite else np.inf})
                    if is_primary and cfg.get("strata", True):
                        dyn = {"region": np.where(geometry(s_te, s_cal, y_cal)["extrap"] > 0, "outside_cal_range", "inside_cal_range")}
                        for lab_name, lab in list(strata_lab.items()) + list(dyn.items()):
                            lt = lab[te] if lab_name != "region" else lab
                            for g in np.unique(lt):
                                msk = lt == g
                                strata_rows.append({"system": sid, "draw": r, "method": meth, "alpha": a, "stratum": lab_name,
                                                    "group": g, "covered": int(cov[msk].sum()), "total": int(msk.sum()),
                                                    "width_sum": float(w[msk].sum()) if finite else np.inf})
                    if is_primary and r == 0 and a == 0.10 and meth in ("P", "B5_jackknife_plus", "B1_gauss_rmse"):
                        for j, tj in enumerate(te):
                            err_rows.append({"system": sid, "method": meth, "ik14": job["ik"][tj], "y": y_te[j],
                                             "lo": lo[j], "hi": hi[j], "covered": bool(cov[j])})
            # ---- decoy filtering (primary n, first decoy_draws draws) ----
            if is_primary and r < cfg.get("decoy_draws", 0) and job.get("decoys") is not None:
                dec = job["decoys"]  # dict: parent (index into system), s, C, ens
                par = dec["parent"]; keep = np.isin(par, te)
                if keep.sum() > 0:
                    par_k = par[keep]; s_d = dec["s"][keep]; Cd = dec["C"].iloc[np.flatnonzero(keep)].reset_index(drop=True)
                    y_par = y[par_k]
                    mu_loo_d = np.vstack([m.predict(s_d) for m in loo]); mu_full_d = full.predict(s_d)
                    dints = {}
                    dints["B1_gauss_rmse"] = lambda a: gaussian_rmse(mu_full_d, loo_res, a)
                    dints["B5_jackknife_plus"] = lambda a: jackknife_plus(mu_loo_d, np.abs(loo_res), None, a)
                    if "B2_pspline_PI" in intervals:
                        dints["B2_pspline_PI"] = lambda a: ps.predict_interval(s_d, a)
                    if "B3_gp_PI" in intervals:
                        dints["B3_gp_PI"] = lambda a: gp.predict_interval(s_d, a)
                    if "P" in sig_cache:
                        sd_te = _sigma_matrix(sig_models["P"], Cd, s_d, s_cal, y_cal, sysf, loo=True)
                        Rp = np.abs(loo_res) / sig_cache["P"][1]
                        dints["P"] = lambda a, sd_te=sd_te, Rp=Rp: jackknife_plus(mu_loo_d, Rp, sd_te, a)
                    for meth, f in dints.items():
                        for a in alphas:
                            lo_d, hi_d = f(a)
                            elim = ~((lo_d <= y_par) & (y_par <= hi_d))
                            dfp = pd.DataFrame({"p": par_k, "e": elim}).groupby("p").e.mean()
                            lo_t, hi_t = intervals[meth](a)
                            cov_t = pd.Series((lo_t <= y_te) & (y_te <= hi_t), index=te)
                            decoy_rows.append({"system": sid, "draw": r, "method": meth, "alpha": a,
                                               "n_parents": int(len(dfp)), "n_decoys": int(len(par_k)),
                                               "decoy_elim_rate": float(dfp.mean()),
                                               "true_retained": float(cov_t.loc[dfp.index].mean())})
    return {"rows": out_rows, "strata": strata_rows, "decoys": decoy_rows, "errors": err_rows, "timing": timing}

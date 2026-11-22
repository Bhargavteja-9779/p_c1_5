"""Offline analyses of XSC from replay summaries (training systems) and per-draw pivot summaries (targets):
  (1) laboratory-level calibration curve: observed share of laboratories with retention < 1−ε vs nominal δ;
  (2) exact (one replay per system) vs draw-averaged multipliers;
  (3) localized XSC: multiplier calibrated on the K training systems most similar to the target in routinely reported
      chromatographic metadata (column class, column length, particle size, temperature, flow, gradient duration, pH class,
      organic modifier, RT range), K = 60 fixed a priori; also the gradient-duration-only kernel variant.
Usage: python experiments/novelty/analyze_xsc_offline.py"""
import os, sys, glob, json, math
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
os.chdir(os.path.join(HERE, "..", ".."))
from conformrt.evaluation.crosssys import C_GRID
from conformrt.visualization.figures import plt, INK2, MUTED, GRID, W1, W2, _save as _fsave

RAW = "../results/raw/xsc"; TAB = "../tables"; FIG = "../figures"
QL = np.linspace(0, 1, 201)
K_LOCAL = 60
meta = pd.read_csv("data/processed/systems.csv", dtype={"system": str}).set_index("system")


def descriptors(ids):
    m = meta.loc[ids]
    num = pd.DataFrame({"lgrad": np.log(m.grad_dur.fillna(meta.grad_dur.median()).clip(lower=1)),
                        "lrt": np.log(m.rt_max.clip(lower=1)), "lcol": np.log(m.col_len.fillna(100).clip(lower=10)),
                        "part": m.particle.fillna(meta.particle.median()), "temp": m.temp.replace(0, np.nan).fillna(35),
                        "flow": m.flow.replace(0, np.nan).fillna(0.3)}, index=ids)
    num = (num - num_mu) / num_sd
    cat = pd.get_dummies(m[["ph_class", "modifier", "usp"]].astype(str)).reindex(columns=cat_cols, fill_value=0).astype(float)
    return np.hstack([num.values, 0.7 * cat.values])


def cov_at(uq, c):
    """retention at multiplier c from 201 pivot quantiles (linear interpolation of the empirical CDF)."""
    return float(np.interp(c, uq, QL, left=0.0, right=1.0))


def lab_mult(V, delta, w=None):
    """laboratory-level multiplier: exact ⌈(1−δ)(S+1)⌉-th order statistic (w=None) or weighted conformal quantile."""
    V = np.asarray(V, float)
    if w is None:
        S = len(V); k = math.ceil((1 - delta) * (S + 1))
        return np.inf if k > S else float(np.sort(V)[k - 1])
    o = np.argsort(V); V, w = V[o], np.asarray(w, float)[o]
    p = np.r_[w, 1.0] / (w.sum() + 1.0)  # test point weight 1 (kernel at distance 0), mass at +∞
    cdf = np.cumsum(p[:-1]); i = np.searchsorted(cdf, 1 - delta)
    return np.inf if i >= len(V) else float(V[i])


def lab_mult_avg(Vdraws, delta):
    S = len(Vdraws)
    F = np.mean([np.searchsorted(np.sort(r), C_GRID, side="right") / len(r) for r in Vdraws], 0)
    ok = np.flatnonzero(F >= (1 - delta) * (1 + 1 / S))
    return float(C_GRID[ok[0]]) if len(ok) else np.inf


if __name__ == "__main__":
    allmeta = meta[meta.method_type == "RP"]
    _num = pd.DataFrame({"lgrad": np.log(allmeta.grad_dur.fillna(allmeta.grad_dur.median()).clip(lower=1)), "lrt": np.log(allmeta.rt_max.clip(lower=1)),
                         "lcol": np.log(allmeta.col_len.fillna(100).clip(lower=10)), "part": allmeta.particle.fillna(allmeta.particle.median()),
                         "temp": allmeta.temp.replace(0, np.nan).fillna(35), "flow": allmeta.flow.replace(0, np.nan).fillna(0.3)})
    num_mu, num_sd = _num.mean(), _num.std().replace(0, 1)
    cat_cols = pd.get_dummies(allmeta[["ph_class", "modifier", "usp"]].astype(str)).columns
    res, curves = [], []
    for f in range(5):
        pv = f"{RAW}/fold{f}_pivots.parquet"
        if not os.path.exists(pv):
            continue
        piv = pd.read_parquet(pv); piv["uq"] = piv.uq.apply(np.asarray)
        sim = pd.read_pickle(f"{RAW}/fold{f}_sim.pkl")
        for (rule, n), P in piv[piv.family == "sigs"].groupby(["rule", "n"]):
            key = (rule, n, "sigs")
            tr = [(sid, o[key]) for sid, o in sim if key in o]
            if not tr:
                continue
            tids = [str(s) for s, _ in tr]; Dtr = descriptors(tids)
            for eps in (0.05, 0.10, 0.20):
                Vd = [o[1][eps] for _, o in tr]; V1 = np.array([v[0] for v in Vd])
                for delta in (0.02, 0.05, 0.10, 0.15, 0.20, 0.30):
                    c_avg = lab_mult_avg(Vd, delta); c_ex = lab_mult(V1, delta)
                    cov = [cov_at(u, c_avg) for u in P.uq]; covx = [cov_at(u, c_ex) for u in P.uq]
                    curves.append({"fold": f, "rule": rule, "n": n, "eps": eps, "delta": delta, "c_avg": c_avg, "c_exact": c_ex,
                                   "fail_avg": float(np.mean(np.array(cov) < 1 - eps)), "fail_exact": float(np.mean(np.array(covx) < 1 - eps)),
                                   "n_labs": len(P)})
            # localized XSC at eps = delta = 0.10 (per target system: K nearest training systems in metadata)
            eps, delta = 0.10, 0.10
            Vmed = np.array([np.median(o[1][eps]) for _, o in tr])  # per-system median requirement (robust)
            Vfirst = np.array([o[1][eps][0] for _, o in tr])
            Vall = [o[1][eps] for _, o in tr]
            c_glob = lab_mult_avg(Vall, delta)
            for sid, Ps in P.groupby("system"):
                dz = np.sqrt(((Dtr - descriptors([str(sid)])) ** 2).sum(1))
                nn = np.argsort(dz)[:K_LOCAL]
                c_loc = lab_mult_avg([Vall[i] for i in nn], delta)
                h = np.median(dz); w = np.exp(-(dz / h) ** 2)
                c_w = lab_mult(Vfirst, delta, w)
                for _, r in Ps.iterrows():
                    for meth, c in (("XSC_lab", c_glob), ("XSC_local_K60", c_loc), ("XSC_kernel", c_w)):
                        res.append({"fold": f, "rule": rule, "n": n, "system": sid, "draw": r.draw, "method": meth, "c": c,
                                    "coverage": cov_at(r.uq, c), "width_median": 2 * c * r.tau * r.sig_med, "span": r.span})
    C = pd.DataFrame(curves); R = pd.DataFrame(res)
    C.to_csv(f"{TAB}/xsc_lab_calibration_curves.csv", index=False); R.to_parquet(f"{RAW}/xsc_offline_local.parquet")
    M = {}
    g = C.groupby(["rule", "eps", "delta"]).apply(lambda d: pd.Series({"fail_avg": np.average(d.fail_avg, weights=d.n_labs), "fail_exact": np.average(d.fail_exact, weights=d.n_labs)}))
    g = g.reset_index(); g.to_csv(f"{TAB}/xsc_lab_calibration_summary.csv", index=False)
    M["lab_calibration"] = {f"{r.rule}|{r.eps}|{r.delta}": [r.fail_avg, r.fail_exact] for r in g.itertuples()}
    L = R.groupby(["rule", "n", "method"]).apply(lambda d: pd.Series({"fail": (d.coverage < 0.9).mean(), "coverage_mean": d.coverage.mean(),
                                                                       "width_median": d.width_median.median()})).reset_index()
    L.to_csv(f"{TAB}/xsc_local_summary.csv", index=False)
    M["local"] = {f"{r.rule}|{r.n}|{r.method}": {"fail": r.fail, "coverage_mean": r.coverage_mean, "width_median": r.width_median} for r in L.itertuples()}
    Lf = R[(R.rule == "random") & (R.n == 20)].groupby(["fold", "method"]).apply(lambda d: pd.Series({"fail": (d.coverage < 0.9).mean(), "width_median": d.width_median.median()})).reset_index()
    Lf.to_csv(f"{TAB}/xsc_local_per_source.csv", index=False)
    M["local_per_source"] = {f"{r.fold}|{r.method}": {"fail": r.fail, "width_median": r.width_median} for r in Lf.itertuples()}
    json.dump(M, open("../results/xsc_offline_metrics.json", "w"), indent=1, default=float)
    print(L[L.rule == "random"].to_string()); print(Lf.to_string())
    # calibration figure: observed vs nominal δ, random rule, all n pooled
    fig, ax = plt.subplots(1, 1, figsize=(W1, 70 / 25.4))
    cols = {0.05: "#eb6834", 0.10: "#2a78d6", 0.20: "#1baf7a"}; mk = {0.05: "s", 0.10: "o", 0.20: "D"}
    for eps in (0.05, 0.10, 0.20):
        d = g[(g.rule == "random") & (g.eps == eps)]
        ax.plot(d.delta, d.fail_avg, marker=mk[eps], ms=4, color=cols[eps], lw=1.2, label=f"retention ≥ {1-eps:.2f}")
    ax.plot([0, 0.32], [0, 0.32], color=MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.set_xlabel("Nominal share of laboratories allowed to fail (δ)"); ax.set_ylabel("Observed share of laboratories failing")
    ax.set_xlim(0, 0.32); ax.set_ylim(0, 0.32); ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    fig.tight_layout(); _fsave(fig, FIG, "figX3_lab_calibration"); plt.close(fig)

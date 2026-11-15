"""Statistical analysis. Unit of analysis = target chromatographic system (draws are averaged within system)."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats

METHOD_LABELS = {
    "B1_gauss_rmse": "Gaussian ±z·RMSE",
    "B2_pspline_PI": "P-spline PI",
    "B3_gp_PI": "GP projection PI",
    "B4_split_conformal": "Split conformal",
    "B5_jackknife_plus": "Jackknife+",
    "B6_jk+_ensemble": "Jackknife+ (ensemble σ)",
    "P": "Proposed (meta-learned σ)",
    "P-novelty": "– novelty features",
    "P-ion": "– ionisation features",
    "P-system": "– system features",
    "P-geometry": "– calibrant-geometry features",
    "P-compound": "– physicochemical features",
    "P_naive_noLOO": "no leave-one-out refit",
    "P_linear_map": "linear instead of monotone map",
}
BASELINES = ["B1_gauss_rmse", "B2_pspline_PI", "B3_gp_PI", "B4_split_conformal", "B5_jackknife_plus", "B6_jk+_ensemble"]
ABLATIONS = ["P-novelty", "P-ion", "P-system", "P-geometry", "P-compound", "P_naive_noLOO", "P_linear_map"]


def per_system(rows: pd.DataFrame) -> pd.DataFrame:
    """Average over draws within system; keeps fold."""
    g = rows.groupby(["system", "fold", "n", "method", "alpha"], as_index=False)
    out = g.agg(coverage=("coverage", "mean"), width=("width_mean", "mean"), width_rel=("width_rel_mean", "mean"),
                winkler=("winkler", "mean"), m_test=("m_test", "mean"), draws=("draw", "nunique"))
    out["gap"] = (1 - out["alpha"]) - out["coverage"]
    return out


def boot_ci(x, stat=np.mean, B=10000, seed=0, clusters=None):
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    if clusters is None:
        idx = rng.integers(0, len(x), size=(B, len(x)))
        bs = np.array([stat(x[i]) for i in idx])
    else:
        clusters = np.asarray(clusters); uc = np.unique(clusters); groups = [x[clusters == c] for c in uc]
        bs = []
        for _ in range(B):
            pick = rng.integers(0, len(uc), len(uc))
            bs.append(stat(np.concatenate([groups[p] for p in pick])))
        bs = np.array(bs)
    return float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def holm(pvals):
    p = np.asarray(pvals, float); order = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
    for k, i in enumerate(order):
        run = max(run, (m - k) * p[i]); adj[i] = min(1.0, run)
    return adj


def rank_biserial(d):
    d = np.asarray(d, float); d = d[d != 0]
    if len(d) == 0:
        return 0.0
    r = stats.rankdata(np.abs(d))
    wp, wm = r[d > 0].sum(), r[d < 0].sum()
    return float((wp - wm) / (wp + wm))


def summary_table(ps: pd.DataFrame, n: int, alpha: float, methods=None, B=10000) -> pd.DataFrame:
    sub = ps[(ps.n == n) & (np.isclose(ps.alpha, alpha))]
    rows = []
    for m in (methods or sorted(sub.method.unique())):
        s = sub[sub.method == m]
        if s.empty:
            continue
        fin = np.isfinite(s.width)
        lo, hi = boot_ci(s.coverage.values, B=B)
        rows.append({"method": m, "label": METHOD_LABELS.get(m, m), "n_systems": len(s),
                     "coverage_mean": s.coverage.mean(), "coverage_ci_lo": lo, "coverage_ci_hi": hi,
                     "coverage_median": s.coverage.median(), "coverage_p10": s.coverage.quantile(0.10),
                     "share_severe_under": float((s.coverage < (1 - alpha) - 0.05).mean()),
                     "share_below_nominal": float((s.coverage < (1 - alpha)).mean()),
                     "width_mean_min": s.width[fin].mean() if fin.any() else np.inf,
                     "width_median_min": s.width[fin].median() if fin.any() else np.inf,
                     "width_rel_median": s.width_rel[fin].median() if fin.any() else np.inf,
                     "winkler_mean": s.winkler[fin].mean() if fin.any() else np.inf,
                     "winkler_median": s.winkler[fin].median() if fin.any() else np.inf,
                     "share_infinite": float((~fin).mean())})
    return pd.DataFrame(rows)


def paired_tests(ps: pd.DataFrame, n: int, alpha: float, ref: str, comparators, metric: str, B=10000) -> pd.DataFrame:
    """Two-sided Wilcoxon signed-rank tests of ref vs each comparator on per-system `metric` (paired by system);
    Holm-adjusted within the family; rank-biserial effect size; bootstrap CI of the median paired difference."""
    sub = ps[(ps.n == n) & (np.isclose(ps.alpha, alpha))]
    a = sub[sub.method == ref].set_index("system")[metric]
    out = []
    for c in comparators:
        b = sub[sub.method == c].set_index("system")[metric]
        j = a.index.intersection(b.index)
        x, y = a.loc[j].values, b.loc[j].values
        ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]
        if len(x) < 5:
            continue
        d = x - y
        try:
            p = stats.wilcoxon(x, y, zero_method="wilcox", alternative="two-sided").pvalue
        except ValueError:
            p = 1.0
        lo, hi = boot_ci(d, stat=np.median, B=B)
        out.append({"ref": ref, "comparator": c, "metric": metric, "n_systems": len(x), "median_ref": np.median(x),
                    "median_comp": np.median(y), "median_diff": np.median(d), "diff_ci_lo": lo, "diff_ci_hi": hi,
                    "share_ref_better": float(np.mean(d < 0)) if metric != "coverage" else float(np.mean(d > 0)),
                    "rank_biserial": rank_biserial(d), "p_wilcoxon": p})
    df = pd.DataFrame(out)
    if len(df):
        df["p_holm"] = holm(df.p_wilcoxon.values)
    return df


def strata_coverage(strata: pd.DataFrame, alpha: float) -> pd.DataFrame:
    s = strata[np.isclose(strata.alpha, alpha)]
    # system-level coverage within stratum, then average over systems (each system equally weighted)
    g = s.groupby(["stratum", "group", "method", "system"], as_index=False)[["covered", "total"]].sum()
    g["cov"] = g.covered / g.total
    agg = g.groupby(["stratum", "group", "method"]).agg(coverage_sysavg=("cov", "mean"), n_systems=("system", "nunique"),
                                                        covered=("covered", "sum"), total=("total", "sum")).reset_index()
    agg["coverage_pooled"] = agg.covered / agg.total
    return agg

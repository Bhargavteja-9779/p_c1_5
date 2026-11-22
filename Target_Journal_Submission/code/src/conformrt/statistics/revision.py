"""Additional analyses: source-level (random-effects) inference, calibration-conditional
(per-draw) coverage, stratum estimators with cluster-bootstrap CIs and a continuous source-adjusted gradient-length model."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats


def source_meta(lr, src) -> dict:
    """Random-effects meta-analysis over sources of per-source mean log-ratios (DerSimonian-Laird tau^2;
    Hartung-Knapp-Sidik-Jonkman t-based CI with k-1 df; prediction interval for a new source with k-2 df)."""
    df = pd.DataFrame({"lr": np.asarray(lr, float), "src": np.asarray(src)}).replace([np.inf, -np.inf], np.nan).dropna()
    g = df.groupby("src").lr
    m = g.mean(); v = (g.var(ddof=1) / g.size()).clip(lower=1e-12); k = len(m)
    w = 1 / v; mu_fe = float((w * m).sum() / w.sum()); Q = float((w * (m - mu_fe) ** 2).sum())
    C = float(w.sum() - (w ** 2).sum() / w.sum()); tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    ws = 1 / (v + tau2); mu = float((ws * m).sum() / ws.sum())
    q = float((ws * (m - mu) ** 2).sum() / (k - 1)); se = float(np.sqrt(max(q, 1e-12) / ws.sum()))
    t = stats.t.ppf(0.975, k - 1)
    tpi = stats.t.ppf(0.975, k - 2) if k > 2 else np.nan
    half_pi = tpi * np.sqrt(tau2 + se ** 2)
    out = {"k_sources": int(k), "n_systems": int(len(df)), "pooled_log_ratio": mu, "ratio": float(np.exp(mu)),
           "ratio_ci_lo": float(np.exp(mu - t * se)), "ratio_ci_hi": float(np.exp(mu + t * se)),
           "ratio_pi_lo": float(np.exp(mu - half_pi)), "ratio_pi_hi": float(np.exp(mu + half_pi)),
           "p_hksj": float(2 * stats.t.sf(abs(mu / se), k - 1)), "tau2": float(tau2),
           "I2": float(max(0.0, (Q - (k - 1)) / Q)) if Q > 0 else 0.0,
           "per_source_ratio": {str(s_): float(np.exp(v_)) for s_, v_ in m.items()},
           "p_wilcoxon_log_ratio": float(stats.wilcoxon(df.lr).pvalue) if (df.lr != 0).any() else 1.0}
    try:  # linear mixed model with a random intercept per source (REML)
        import statsmodels.formula.api as smf
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = smf.mixedlm("lr ~ 1", df, groups=df["src"]).fit(reml=True)
        lo, hi = fit.conf_int().loc["Intercept"]
        out.update({"lmm_ratio": float(np.exp(fit.params["Intercept"])), "lmm_ci_lo": float(np.exp(lo)), "lmm_ci_hi": float(np.exp(hi)),
                    "lmm_p": float(fit.pvalues["Intercept"])})
    except Exception as e:  # pragma: no cover
        out["lmm_error"] = str(e)[:200]
    return out


def draw_coverage(rows: pd.DataFrame, methods) -> pd.DataFrame:
    """Calibration-conditional view: distribution over (system, draw) of the coverage obtained with ONE calibrant set."""
    out = []
    for (n, al, m), g in rows[rows.method.isin(methods)].groupby(["n", "alpha", "method"]):
        c = g.coverage.values; nom = 1 - al
        out.append({"n": int(n), "alpha": float(al), "method": m, "draws": int(len(c)), "q05": float(np.quantile(c, 0.05)),
                    "q10": float(np.quantile(c, 0.10)), "median": float(np.median(c)), "share_below_nom_minus_0.05": float((c < nom - 0.05).mean()),
                    "share_below_nom_minus_0.10": float((c < nom - 0.10).mean())})
    return pd.DataFrame(out)


def split_conformal_reference(n: int, alpha: float) -> dict:
    """Exact calibration-conditional coverage of split conformal with n_cal = n // 2 exchangeable calibration points:
    coverage | calibration set ~ Beta(k, n_cal + 1 - k), k = ceil((1 - alpha)(n_cal + 1)) (infinite test set)."""
    n_cal = n // 2; k = int(np.ceil((1 - alpha) * (n_cal + 1)))
    if k > n_cal:
        return {"n_cal": n_cal, "k": k, "finite": False}
    b = stats.beta(k, n_cal + 1 - k); nom = 1 - alpha
    return {"n_cal": n_cal, "k": k, "finite": True, "q05": float(b.ppf(0.05)), "q10": float(b.ppf(0.10)), "median": float(b.median()),
            "share_below_nom_minus_0.05": float(b.cdf(nom - 0.05)), "share_below_nom_minus_0.10": float(b.cdf(nom - 0.10))}


def strata_estimators(strata: pd.DataFrame, alpha: float, methods, B: int = 2000, seed: int = 0) -> pd.DataFrame:
    """System-average and pooled coverage per stratum with bootstrap CIs over systems (and the source-cluster range)."""
    s = strata[np.isclose(strata.alpha, alpha) & strata.method.isin(methods)]
    g = s.groupby(["stratum", "group", "method", "system", "fold"], as_index=False).agg(covered=("covered", "sum"), total=("total", "sum"),
                                                                                         draws=("draw", "nunique"))
    rng = np.random.default_rng(seed); out = []
    for (st, gr, m), d in g.groupby(["stratum", "group", "method"]):
        d = d[d.total > 0]
        if d.empty:
            continue
        cov_sys = (d.covered / d.total).values; cv, tt = d.covered.values, d.total.values
        idx = rng.integers(0, len(d), size=(B, len(d)))
        bs_sys = cov_sys[idx].mean(axis=1); bs_pool = cv[idx].sum(axis=1) / tt[idx].sum(axis=1)
        per_src = d.groupby("fold").apply(lambda x: x.covered.sum() / x.total.sum(), include_groups=False)
        big = d[(d.total / d.draws) >= 10]
        out.append({"stratum": st, "group": gr, "method": m, "systems": int(len(d)), "median_per_draw_count": float(np.median(d.total / d.draws)),
                    "sysavg": float(cov_sys.mean()), "sysavg_ci_lo": float(np.percentile(bs_sys, 2.5)), "sysavg_ci_hi": float(np.percentile(bs_sys, 97.5)),
                    "pooled": float(cv.sum() / tt.sum()), "pooled_ci_lo": float(np.percentile(bs_pool, 2.5)), "pooled_ci_hi": float(np.percentile(bs_pool, 97.5)),
                    "pooled_source_min": float(per_src.min()), "pooled_source_max": float(per_src.max()),
                    "systems_ge10_per_draw": int(len(big)), "sysavg_ge10": float((big.covered / big.total).mean()) if len(big) else np.nan})
    return pd.DataFrame(out)


def span_model(lr, log_span, src) -> dict:
    """Linear mixed model lr ~ log(span) + (1 | source) on RepoRT (exploratory, defined after the results were seen)."""
    import statsmodels.formula.api as smf
    import warnings
    df = pd.DataFrame({"lr": lr, "x": log_span, "src": src}).replace([np.inf, -np.inf], np.nan).dropna()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = smf.mixedlm("lr ~ x", df, groups=df["src"]).fit(reml=True)
    ci = fit.conf_int()
    return {"n_systems": int(len(df)), "intercept": float(fit.params["Intercept"]), "slope": float(fit.params["x"]),
            "slope_ci_lo": float(ci.loc["x", 0]), "slope_ci_hi": float(ci.loc["x", 1]), "slope_p": float(fit.pvalues["x"]),
            "x_mean": float(df.x.mean())}


def external_span_check(lr, log_span, model: dict) -> dict:
    """Compare an external set with the RepoRT span model: observed vs predicted mean log-ratio, and the external slope (OLS)."""
    lr = np.asarray(lr, float); x = np.asarray(log_span, float); ok = np.isfinite(lr) & np.isfinite(x)
    lr, x = lr[ok], x[ok]
    pred = model["intercept"] + model["slope"] * x
    res = stats.linregress(x, lr) if len(lr) > 3 else None
    se = lr.std(ddof=1) / np.sqrt(len(lr))
    long_ = x > np.median(x)
    return {"n_systems": int(len(lr)), "observed_mean_lr": float(lr.mean()), "observed_ci_lo": float(lr.mean() - 1.96 * se),
            "observed_ci_hi": float(lr.mean() + 1.96 * se), "predicted_mean_lr": float(pred.mean()),
            "slope": float(res.slope) if res else np.nan, "slope_p": float(res.pvalue) if res else np.nan,
            "ratio_long": float(np.exp(np.median(lr[long_]))), "ratio_short": float(np.exp(np.median(lr[~long_])))}

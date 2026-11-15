"""Derive every number quoted in the manuscript from machine-generated result files.
build_facts(results_dir, tables_dir) -> (facts: dict[str, value], tables: dict[str, markdown_table])"""
import os, json
import numpy as np, pandas as pd

LAB = {"B1_gauss_rmse": "Gaussian ±z·RMSE", "B2_pspline_PI": "P-spline PI", "B3_gp_PI": "GP projection PI",
       "B4_split_conformal": "Split conformal", "B5_jackknife_plus": "Jackknife+", "B6_jk+_ensemble": "Jackknife+ (ensemble σ)",
       "P": "**Proposed**", "FIX_0.25min": "Fixed ±0.25 min", "FIX_0.5min": "Fixed ±0.5 min", "FIX_1min": "Fixed ±1 min",
       "FIX_2min": "Fixed ±2 min", "FIX_5pct": "Fixed ±5 %", "FIX_10pct": "Fixed ±10 %", "RMSE_x2": "±2·RMSE", "RMSE_x3": "±3·RMSE",
       "P-novelty": "– novelty features", "P-ion": "– ionisation features", "P-system": "– system features",
       "P-geometry": "– calibrant-geometry features", "P-compound": "– physicochemical features",
       "P_naive_noLOO": "no leave-one-out refit", "P_linear_map": "linear map"}


def _f(x, d=3):
    if x is None or (isinstance(x, float) and (np.isnan(x))):
        return "–"
    if isinstance(x, float) and np.isinf(x):
        return "∞"
    return f"{x:.{d}f}"


def _p(x):
    return "< 0.001" if x < 0.001 else f"{x:.3f}"


def build_facts(RES, TAB):
    F, T = {}, {}
    m = json.load(open(os.path.join(RES, "metrics.json"))); st = json.load(open(os.path.join(RES, "statistics.json")))
    d = m["dataset"]; F["n_target_systems"] = d["n_target_systems"]; F["n_rt_values_targets"] = f"{d['n_rt_values_targets']:,}"
    F["median_compounds"] = f"{d['median_compounds_per_system']:.0f}"
    F["fold_list"] = "; ".join(f"{r['fold_name']}: {r['systems']} systems" for r in d["per_fold"])
    for r in m["main_n20_a010"]:
        for k, v in r.items():
            F[f"main.{r['method']}.{k}"] = v
    for r in st.get("paired_tests_n20_a010", []):
        for k, v in r.items():
            F[f"tests.{r['metric']}.{r['comparator']}.{k}"] = v
    for r in st.get("ablation", []):
        for k, v in r.items():
            F[f"abl.{r['method']}.{k}"] = v
    for r in st.get("strata", []):
        F[f"strata.{r['stratum']}.{r['group']}.{r['method']}"] = r["coverage_sysavg"]
        F[f"strata_n.{r['stratum']}.{r['group']}"] = r["total"]
    for r in m.get("candidate_filtering", []):
        F[f"dec.{r['method']}.{r['alpha']:.2f}.elim"] = r["decoy_elim_mean"]; F[f"dec.{r['method']}.{r['alpha']:.2f}.ret"] = r["true_retained_mean"]
    for r in st.get("coverage_validity", []):
        F[f"val.{r['method']}.n{int(r['n'])}.a{r['alpha']:.2f}"] = r["mean_cov"]
        F[f"valshare.{r['method']}.n{int(r['n'])}.a{r['alpha']:.2f}"] = r["share_systems_below_nominal_minus_0.05"]
    # tables (markdown) ------------------------------------------------------------
    t2 = pd.read_csv(os.path.join(TAB, "table2_main_n20_a010.csv"))
    rows = ["| Method | Mean coverage [95 % CI] | Systems with coverage < 0.85 | Median width (min) | Median width (% of RT span) | Mean interval score (min) |", "|---|---|---|---|---|---|"]
    for _, r in t2.iterrows():
        rows.append(f"| {LAB.get(r.method, r.method)} | {_f(r.coverage_mean)} [{_f(r.coverage_ci_lo)}, {_f(r.coverage_ci_hi)}] | "
                    f"{100*r.share_severe_under:.0f} % | {_f(r.width_median_min, 2)} | {_f(100*r.width_rel_median, 1)} | {_f(r.winkler_mean, 2)} |")
    T["table2"] = "\n".join(rows)
    return F, T

"""Generate all tables, figures and statistics from raw experiment outputs (no hand-typed numbers).
Usage: python experiments/analyze.py --raw ../results/raw --out ../ [--B 10000]"""
import argparse, os, sys, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from conformrt.statistics.analysis import (one_sided_coverage_test, per_system, summary_table, paired_tests, strata_coverage, boot_ci, rank_biserial,
                                           METHOD_LABELS, BASELINES, ABLATIONS)
from conformrt.visualization import figures as F

ap = argparse.ArgumentParser(); ap.add_argument("--raw", default="../results/raw"); ap.add_argument("--out", default="..")
ap.add_argument("--processed", default="data/processed"); ap.add_argument("--B", type=int, default=10000)
a = ap.parse_args()
RES, TAB, FIG = (os.path.join(a.out, d) for d in ("results", "tables", "figures"))
for d in (RES, TAB, FIG):
    os.makedirs(d, exist_ok=True)
PRIMARY_N, PRIMARY_A = 20, 0.1
ORDER = ["P", "B1_gauss_rmse", "B5_jackknife_plus", "B4_split_conformal", "B3_gp_PI", "B7_gp_meta", "B2_pspline_PI", "B6_jk+_ensemble"]
L = METHOD_LABELS


def load(tag, kind):
    p = os.path.join(a.raw, f"{tag}_{kind}.parquet")
    return pd.read_parquet(p) if os.path.exists(p) else pd.DataFrame()


def md(df, path, floatfmt=3):
    df.to_csv(path + ".csv", index=False)
    d = df.copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].map(lambda v: "∞" if np.isinf(v) else (f"{v:.{floatfmt}g}" if abs(v) >= 1000 else f"{v:.{floatfmt}f}"))
    with open(path + ".md", "w") as f:
        f.write("| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n")
        for _, r in d.iterrows():
            f.write("| " + " | ".join(map(str, r.values)) + " |\n")


metrics, stats_out = {}, {}
main_rows = load("main", "rows")
_mg = load("metagp", "rows")  # Addendum E: comparator B7 evaluated post hoc on the identical calibrant draws
if not _mg.empty:
    main_rows = pd.concat([main_rows, _mg[main_rows.columns.intersection(_mg.columns)]], ignore_index=True)
if main_rows.empty:
    sys.exit("no main results yet")
ps = per_system(main_rows)
allps = [ps.assign(experiment="main", base="order", strategy="random", perturb="none")]

# ---------------- Table 1: dataset ----------------
meta = pd.read_csv(os.path.join(a.processed, "systems.csv"), dtype={"system": str})
from conformrt.evaluation.protocol import fold_of_source, FOLD_NAMES
tsys = ps.system.unique()
mt = meta[meta.system.isin(tsys)].copy(); mt["fold"] = mt.source.map(fold_of_source); mt["fold_name"] = mt.fold.map(FOLD_NAMES)
mt["rt_span"] = mt.rt_max - mt.rt_min  # RT span = last minus first eluting compound
t1 = mt.groupby("fold_name").agg(systems=("system", "size"), compounds_median=("n_compounds", "median"),
                                  compounds_total=("n_compounds", "sum"), rt_span_median_min=("rt_span", "median"),
                                  acidic=("ph_class", lambda x: int((x == "acidic").sum())),
                                  non_acidic=("ph_class", lambda x: int((x != "acidic").sum())),
                                  acn=("modifier", lambda x: int(x.str.startswith("acn").sum())),
                                  meoh=("modifier", lambda x: int(x.str.startswith("meoh").sum()))).reset_index()
md(t1, os.path.join(TAB, "table1_dataset"))
F.fig_dataset(mt, FIG)
metrics["dataset"] = {"n_target_systems": int(len(mt)), "n_rt_values_targets": int(mt.n_compounds.sum()),
                      "median_compounds_per_system": float(mt.n_compounds.median()), "per_fold": t1.to_dict("records")}

# ---------------- Table 2: main results ----------------
methods_main = [m for m in ORDER if m in ps.method.unique()]
t2 = summary_table(ps, PRIMARY_N, PRIMARY_A, methods_main, B=a.B)
md(t2, os.path.join(TAB, "table2_main_n20_a010"))
metrics["main_n20_a010"] = t2.to_dict("records")
F.fig_main(ps, FIG, PRIMARY_N, PRIMARY_A, methods_main, L)

# all n x alpha summary (supplement)
full = []
for n in sorted(ps.n.unique()):
    for al in sorted(ps.alpha.unique()):
        s = summary_table(ps, n, al, [m for m in ORDER + ABLATIONS if m in ps.method.unique()], B=min(a.B, 2000))
        s.insert(0, "alpha", al); s.insert(0, "n", n); full.append(s)
full = pd.concat(full, ignore_index=True); md(full, os.path.join(TAB, "tableS1_all_n_alpha"))
F.fig_calibrants(full, FIG, [m for m in ORDER if m in ps.method.unique()], L, alpha=PRIMARY_A)
F.fig_alpha_calibration(full, FIG, [m for m in ["P", "B5_jackknife_plus", "B4_split_conformal", "B3_gp_PI", "B7_gp_meta", "B2_pspline_PI", "B1_gauss_rmse"] if m in ps.method.unique()], L, n=PRIMARY_N)

# coverage validity (system-bootstrap and source-cluster bootstrap) for every n/alpha
val = []
for (n, al, m), s in ps.groupby(["n", "alpha", "method"]):
    lo, hi = boot_ci(s.coverage.values, B=min(a.B, 4000))
    clo, chi = boot_ci(s.coverage.values, B=min(a.B, 4000), clusters=s.fold.values)
    ub, p1 = one_sided_coverage_test(s.coverage.values, 1 - al, B=a.B)  # pre-specified: H1 mean coverage < 1 - alpha
    cub, cp1 = one_sided_coverage_test(s.coverage.values, 1 - al, B=min(a.B, 4000), clusters=s.fold.values)
    val.append({"n": n, "alpha": al, "method": m, "mean_cov": s.coverage.mean(), "ci_lo": lo, "ci_hi": hi, "cluster_ci_lo": clo,
                "cluster_ci_hi": chi, "upper95_one_sided": ub, "p_under_one_sided": p1, "cluster_upper95_one_sided": cub,
                "cluster_p_under_one_sided": cp1, "nominal": 1 - al, "mean_minus_nominal": s.coverage.mean() - (1 - al),
                "share_systems_below_nominal_minus_0.05": float((s.coverage < 1 - al - 0.05).mean())})
val = pd.DataFrame(val); md(val, os.path.join(TAB, "tableS2_coverage_validity"))
stats_out["coverage_validity"] = val.to_dict("records")

# ---------------- Table 3: paired tests ----------------
comps = [m for m in BASELINES if m in ps.method.unique()]
t3 = pd.concat([paired_tests(ps, PRIMARY_N, PRIMARY_A, "P", comps, met, B=a.B) for met in ["winkler", "width", "gap"]], ignore_index=True)
t3["comparator_label"] = t3.comparator.map(L)
md(t3, os.path.join(TAB, "table3_paired_tests_n20_a010"))
stats_out["paired_tests_n20_a010"] = t3.to_dict("records")
# all n / alpha paired tests on winkler (supplement)
t3s = []
for n in sorted(ps.n.unique()):
    for al in sorted(ps.alpha.unique()):
        cm = [m for m in BASELINES if m in ps[(ps.n == n)].method.unique()]
        r = paired_tests(ps, n, al, "P", cm, "winkler", B=1000)
        if len(r):
            r.insert(0, "alpha", al); r.insert(0, "n", n); t3s.append(r)
if t3s:
    t3s = pd.concat(t3s, ignore_index=True); md(t3s, os.path.join(TAB, "tableS3_paired_winkler_all_n_alpha"))
    stats_out["paired_winkler_all"] = t3s.to_dict("records")
# per-fold consistency
pf = []
for f, s in ps[(ps.n == PRIMARY_N) & np.isclose(ps.alpha, PRIMARY_A)].groupby("fold"):
    for m in methods_main:
        q = s[s.method == m]
        pf.append({"fold": int(f), "fold_name": FOLD_NAMES[int(f)], "method": m, "systems": len(q), "coverage": q.coverage.mean(),
                   "winkler": q.winkler.mean(), "width_median": q.width.median()})
pf = pd.DataFrame(pf); md(pf, os.path.join(TAB, "tableS4_per_fold_n20_a010")); stats_out["per_fold"] = pf.to_dict("records")

# ---------------- Table 4: ablation ----------------
abl_m = [m for m in ABLATIONS if m in ps.method.unique()]
if abl_m:
    t4a = paired_tests(ps, PRIMARY_N, PRIMARY_A, "P", abl_m, "winkler", B=a.B)
    # paired_tests returns ref - comparator; for the ablation figure we want comparator - ref (positive = ablation worse)
    abl = t4a.rename(columns={"comparator": "method"}).copy()
    abl["median_diff"], abl["diff_ci_lo"], abl["diff_ci_hi"] = -abl.median_diff, -abl.diff_ci_hi, -abl.diff_ci_lo
    abl["rank_biserial"] = -abl.rank_biserial
    cov = ps[(ps.n == PRIMARY_N) & np.isclose(ps.alpha, PRIMARY_A)].groupby("method").coverage.mean()
    abl["coverage_mean"] = abl.method.map(cov); abl["label"] = abl.method.map(L)
    full_row = pd.DataFrame([{"method": "P", "label": L["P"], "median_diff": 0.0, "diff_ci_lo": 0.0, "diff_ci_hi": 0.0, "coverage_mean": cov.get("P")}])
    abl_plot = pd.concat([full_row, abl], ignore_index=True)
    md(abl[["method", "label", "n_systems", "median_diff", "diff_ci_lo", "diff_ci_hi", "rank_biserial", "p_wilcoxon", "p_holm", "coverage_mean"]],
       os.path.join(TAB, "table4_ablation_n20_a010"))
    stats_out["ablation"] = abl.to_dict("records")
    F.fig_ablation(abl_plot, FIG, L)

# ---------------- strata / conditional coverage ----------------
st = load("main", "strata")
if not st.empty:
    sc = strata_coverage(st, PRIMARY_A); md(sc, os.path.join(TAB, "tableS5_strata_n20_a010"))
    stats_out["strata"] = sc.to_dict("records")
    order = [("seen", ["seen", "unseen"]), ("similarity", ["sim>=0.7", "sim0.4-0.7", "sim<0.4"]), ("ionisable", ["neutral", "ionisable"]),
             ("region", ["inside_cal_range", "outside_cal_range"]), ("rt_tertile", ["early", "middle", "late"])]
    F.fig_strata(sc, FIG, [m for m in ["P", "B1_gauss_rmse", "B5_jackknife_plus"] if m in sc.method.unique()], L, PRIMARY_A, order)
    sup = sc[sc.stratum == "superclass"].pivot_table(index="group", columns="method", values="coverage_sysavg")
    sup.reset_index().to_csv(os.path.join(TAB, "tableS6_superclass_coverage.csv"), index=False)
    # strata-specific worst-case coverage per method (max deviation below nominal)
    worst = sc[sc.stratum != "superclass"].groupby("method").coverage_sysavg.min().rename("worst_stratum_coverage").reset_index()
    md(worst, os.path.join(TAB, "tableS7_worst_stratum")); stats_out["worst_stratum"] = worst.to_dict("records")

# ---------------- candidate filtering ----------------
dec = load("main", "decoys")
if not dec.empty:
    dsys = dec.groupby(["system", "method", "alpha"]).agg(decoy_elim_rate=("decoy_elim_rate", "mean"), true_retained=("true_retained", "mean"),
                                                         n_parents=("n_parents", "mean"), n_decoys=("n_decoys", "mean")).reset_index()
    t7 = dsys.groupby(["method", "alpha"]).agg(systems=("system", "nunique"), decoy_elim_mean=("decoy_elim_rate", "mean"),
                                               decoy_elim_median=("decoy_elim_rate", "median"), true_retained_mean=("true_retained", "mean"),
                                               parents_per_system=("n_parents", "mean"), decoys_per_system=("n_decoys", "mean")).reset_index()
    t7["label"] = t7.method.map(L); md(t7, os.path.join(TAB, "table7_candidate_filtering"))
    metrics["candidate_filtering"] = t7.to_dict("records")
    F.fig_decoys(t7.rename(columns={"decoy_elim_mean": "decoy_elim_rate", "true_retained_mean": "true_retained"}), FIG,
                 [m for m in ["P", "B1_gauss_rmse", "B5_jackknife_plus", "B3_gp_PI", "B2_pspline_PI"] if m in t7.method.unique()], L)
    # paired test: decoy elimination P vs jackknife+ at matched alpha (system level)
    dt = []
    for al in sorted(dsys.alpha.unique()):
        A = dsys[(dsys.method == "P") & np.isclose(dsys.alpha, al)].set_index("system")
        for c in [m for m in ["B5_jackknife_plus", "B1_gauss_rmse", "B2_pspline_PI", "B3_gp_PI"] if m in dsys.method.unique()]:
            Bm = dsys[(dsys.method == c) & np.isclose(dsys.alpha, al)].set_index("system"); j = A.index.intersection(Bm.index)
            if len(j) < 5:
                continue
            from scipy import stats as _st
            d = A.loc[j, "decoy_elim_rate"].values - Bm.loc[j, "decoy_elim_rate"].values
            dt.append({"alpha": al, "comparator": c, "n_systems": len(j), "median_diff_elim": float(np.median(d)),
                       "rank_biserial": rank_biserial(d), "p_wilcoxon": float(_st.wilcoxon(d).pvalue) if np.any(d != 0) else 1.0})
    if dt:
        dt = pd.DataFrame(dt); md(dt, os.path.join(TAB, "tableS8_decoy_tests")); stats_out["decoy_tests"] = dt.to_dict("records")

# ---------------- errors (per-compound) ----------------
err = load("main", "errors")
if not err.empty:
    comp = pd.read_parquet(os.path.join(a.processed, "compounds.parquet"))[["ik14", "superclass", "cls", "n_systems"]]
    e = err.merge(comp, on="ik14", how="left")
    e["superclass"] = e.superclass.astype(str).str.replace(r" \(CHEMONTID.*", "", regex=True)
    e["cls"] = e.cls.astype(str).str.replace(r" \(CHEMONTID.*", "", regex=True)
    miss = e[(e.method == "P") & (~e.covered)]
    top = miss.groupby("cls").size().sort_values(ascending=False).head(15).rename("misses").reset_index()
    tot = e[e.method == "P"].groupby("cls").size().rename("total")
    top["miss_rate"] = top.cls.map(lambda c: top.set_index("cls").misses[c] / tot.get(c, np.nan))
    md(top, os.path.join(TAB, "tableS9_error_classes"))
    side = miss.assign(side=np.where(miss.y < miss.lo, "eluted_earlier_than_window", "eluted_later_than_window")).side.value_counts()
    metrics["error_analysis"] = {"P_misses_draw0": int(len(miss)), "P_evaluated_draw0": int((e.method == "P").sum()), "side": side.to_dict()}

# ---------------- timing ----------------
tm = load("main", "timing")
if not tm.empty:
    t8 = tm.groupby("n").agg(draws=("draw", "size"), t_map_loo_ms=("t_map_loo_s", lambda x: 1000 * x.median()),
                             t_sigma_ms=("t_sigma_s", lambda x: 1000 * x.median()), m_test_median=("m_test", "median")).reset_index()
    md(t8, os.path.join(TAB, "table8_compute"))
    metrics["compute_per_draw"] = t8.to_dict("records")
cb = os.path.join(a.out, "results", "experiment_logs", "compute_benchmark.json")
if os.path.exists(cb):
    metrics["compute_benchmark"] = json.load(open(cb))
man = sorted(glob.glob(os.path.join(a.out, "results", "experiment_logs", "manifest_main_*.json")))
if man:
    metrics["run_manifest"] = json.load(open(man[-1]))

# ---------------- robustness ----------------
rob = []
for tag in sorted(glob.glob(os.path.join(a.raw, "robust_*_rows.parquet"))):
    t = os.path.basename(tag).replace("_rows.parquet", ""); r = per_system(pd.read_parquet(tag)); allps.append(r.assign(experiment=t))
    s = summary_table(r, PRIMARY_N, PRIMARY_A, [m for m in ORDER if m in r.method.unique()], B=1000); s.insert(0, "scenario", t); rob.append(s)
    pt = paired_tests(r, PRIMARY_N, PRIMARY_A, "P", [m for m in BASELINES if m in r.method.unique()], "winkler", B=1000)
    stats_out.setdefault("robustness_tests", {})[t] = pt.to_dict("records")
if rob:
    base_s = t2.copy(); base_s.insert(0, "scenario", "robust_random_none (main)")
    rob = pd.concat([base_s] + rob, ignore_index=True)
    names = {"robust_random_none (main)": "random calibrants (main)", "robust_span_none": "RT-spanning calibrants",
             "robust_common_none": "commonly measured calibrants", "robust_random_noise": "calibrant RT noise (2% of span)",
             "robust_random_outlier": "one outlying calibrant (+10% span)"}
    rob["scenario"] = rob.scenario.map(lambda x: names.get(x, x))
    md(rob, os.path.join(TAB, "table5_robustness")); metrics["robustness"] = rob.to_dict("records")
    F.fig_robustness(rob, FIG, [m for m in ["P", "B1_gauss_rmse", "B5_jackknife_plus"] if m in rob.method.unique()], L)

# robustness mechanism: share of test compounds outside the calibrant score range and coverage there, per scenario
mech = []
for tag in ["main"] + [os.path.basename(x).replace("_strata.parquet", "") for x in sorted(glob.glob(os.path.join(a.raw, "robust_*_strata.parquet")))]:
    pth = os.path.join(a.raw, f"{tag}_strata.parquet")
    if not os.path.exists(pth):
        continue
    s_ = pd.read_parquet(pth, filters=[("stratum", "in", ["region", "rt_tertile"])])
    s_ = s_[np.isclose(s_.alpha, PRIMARY_A)]
    reg = s_[(s_.stratum == "region") & (s_.method == "P")].groupby("group").total.sum()
    cov_ = strata_coverage(s_, PRIMARY_A)
    for meth in ["P", "B5_jackknife_plus", "B1_gauss_rmse"]:
        g_ = cov_[cov_.method == meth].set_index(["stratum", "group"]).coverage_sysavg
        mech.append({"scenario": tag, "method": meth, "share_outside_cal_range": float(reg.get("outside_cal_range", 0) / reg.sum()),
                     "cov_outside": float(g_.get(("region", "outside_cal_range"), np.nan)), "cov_inside": float(g_.get(("region", "inside_cal_range"), np.nan)),
                     "cov_late": float(g_.get(("rt_tertile", "late"), np.nan))})
if mech:
    mech = pd.DataFrame(mech); md(mech, os.path.join(TAB, "tableS16_robustness_mechanism")); metrics["robustness_mechanism"] = mech.to_dict("records")

# ---------------- generalisation ----------------
gen = []
for tag, name in [("gen_smrt", "SMRT-regression base model (RP)"), ("gen_hilic", "HILIC systems (order model)")]:
    r = load(tag, "rows")
    if r.empty:
        continue
    r = per_system(r); allps.append(r.assign(experiment=tag))
    for n in sorted(r.n.unique()):
        s = summary_table(r, n, PRIMARY_A, [m for m in ORDER if m in r.method.unique()], B=1000); s.insert(0, "n", n); s.insert(0, "setting", name); gen.append(s)
    pt = paired_tests(r, PRIMARY_N, PRIMARY_A, "P", [m for m in BASELINES if m in r.method.unique()], "winkler", B=1000)
    stats_out.setdefault("generalisation_tests", {})[tag] = pt.to_dict("records")
if gen:
    gen = pd.concat(gen, ignore_index=True); md(gen, os.path.join(TAB, "table6_generalisation")); metrics["generalisation"] = gen.to_dict("records")

# ---------------- sensitivity (difficulty-model hyper-parameters) ----------------
sens = []
for tag in sorted(glob.glob(os.path.join(a.raw, "sens_*_rows.parquet"))):
    t = os.path.basename(tag).replace("_rows.parquet", ""); r = per_system(pd.read_parquet(tag))
    s = summary_table(r, PRIMARY_N, PRIMARY_A, ["P", "B5_jackknife_plus"], B=1000); s.insert(0, "setting", t); sens.append(s)
if sens:
    sens = pd.concat(sens, ignore_index=True); md(sens, os.path.join(TAB, "tableS10_sensitivity")); metrics["sensitivity"] = sens.to_dict("records")

# ---------------- heterogeneity of the gain across source folds and its driver ----------------
from scipy import stats as _st
sub20 = ps[(ps.n == PRIMARY_N) & np.isclose(ps.alpha, PRIMARY_A)]
Pw = sub20[sub20.method == "P"].set_index("system")
het = []
for comp in ["B5_jackknife_plus", "B1_gauss_rmse", "B3_gp_PI"]:
    Cw = sub20[sub20.method == comp].set_index("system")
    for f, g in Pw.groupby("fold"):
        j = g.index.intersection(Cw.index); d = g.loc[j, "winkler"] - Cw.loc[j, "winkler"]; rel = g.loc[j, "winkler"] / Cw.loc[j, "winkler"] - 1
        het.append({"comparator": comp, "fold": int(f), "fold_name": FOLD_NAMES[int(f)], "systems": len(j), "median_diff_min": float(d.median()),
                    "median_rel_change": float(rel.median()), "share_P_better": float((d < 0).mean()), "rank_biserial": rank_biserial(d.values),
                    "p_wilcoxon": float(_st.wilcoxon(d).pvalue) if (d != 0).any() else 1.0})
    j = Pw.index.intersection(Cw.index); rel = Pw.loc[j, "winkler"] / Cw.loc[j, "winkler"] - 1
    het.append({"comparator": comp, "fold": -1, "fold_name": "all", "systems": len(j), "median_diff_min": float((Pw.loc[j, "winkler"] - Cw.loc[j, "winkler"]).median()),
                "median_rel_change": float(rel.median()), "rel_q1": float(rel.quantile(.25)), "rel_q3": float(rel.quantile(.75)),
                "share_P_better": float((rel < 0).mean())})
het = pd.DataFrame(het); md(het, os.path.join(TAB, "table11_fold_heterogeneity")); stats_out["fold_heterogeneity"] = het.to_dict("records")
Cj = sub20[sub20.method == "B5_jackknife_plus"].set_index("system"); j = Pw.index.intersection(Cj.index)
relj = Pw.loc[j, "winkler"] / Cj.loc[j, "winkler"] - 1
span = (meta.set_index("system").loc[j, "rt_max"] - meta.set_index("system").loc[j, "rt_min"]).astype(float)
rho = _st.spearmanr(relj.values, span.values)
fold_j = Pw.loc[j, "fold"].astype(int).values
within = {FOLD_NAMES[int(f)]: float(_st.spearmanr(relj.values[fold_j == f], span.values[fold_j == f]).statistic) for f in np.unique(fold_j)}
# fold-adjusted partial Spearman correlation: correlate within-fold-demeaned ranks
r1 = _st.rankdata(relj.values); r2 = _st.rankdata(span.values)
for f in np.unique(fold_j):
    r1[fold_j == f] -= r1[fold_j == f].mean(); r2[fold_j == f] -= r2[fold_j == f].mean()
pr = float(np.corrcoef(r1, r2)[0, 1]); dfree = len(r1) - 2 - (len(np.unique(fold_j)) - 1)
pp = float(2 * _st.t.sf(abs(pr) * np.sqrt(dfree / (1 - pr ** 2)), dfree))
stats_out["gain_vs_rt_span"] = {"spearman_rho": float(rho.statistic), "p": float(rho.pvalue), "n_systems": int(len(j)),
                                "within_fold_rho": within, "within_fold_rho_min": min(within.values()), "within_fold_rho_max": max(within.values()),
                                "partial_rho_fold_adjusted": pr, "partial_p": pp,
                                "median_rel_change_span_below_median": float(relj[span <= span.median()].median()),
                                "median_rel_change_span_above_median": float(relj[span > span.median()].median()), "median_span_min": float(span.median())}

# ---------------- decoy elimination at matched true-candidate retention (interpolated on the alpha grid) ----------------
if not dec.empty:
    mr = []
    for m_, g in t7.groupby("method"):
        g = g.sort_values("true_retained_mean")
        for target in [0.85, 0.90]:
            x, y_ = g.true_retained_mean.values, g.decoy_elim_mean.values
            val = float(np.interp(target, x, y_)) if x.min() <= target <= x.max() else np.nan
            mr.append({"method": m_, "label": L.get(m_, m_), "retention": target, "decoy_elim_interp": val})
    mr = pd.DataFrame(mr); md(mr, os.path.join(TAB, "tableS14_matched_retention")); metrics["matched_retention"] = mr.to_dict("records")

# ---------------- fixed-tolerance windows (Addendum B; same calibrant draws) ----------------
fw = load("fixedwin", "rows")
if not fw.empty:
    fwp = per_system(fw); allps.append(fwp.assign(experiment="fixedwin"))
    comb = pd.concat([ps[(ps.method.isin(["P", "B5_jackknife_plus", "B1_gauss_rmse"])) & np.isclose(ps.alpha, PRIMARY_A)], fwp], ignore_index=True)
    rows_fw = []
    for n in sorted(fwp.n.unique()):
        s_ = summary_table(comb, n, PRIMARY_A, ["P", "B5_jackknife_plus", "B1_gauss_rmse", "RMSE_x2", "RMSE_x3", "FIX_0.25min", "FIX_0.5min", "FIX_1min", "FIX_2min", "FIX_5pct", "FIX_10pct"], B=1000)
        s_.insert(0, "n", n); rows_fw.append(s_)
    tfw = pd.concat(rows_fw, ignore_index=True); md(tfw, os.path.join(TAB, "table9_fixed_windows")); metrics["fixed_windows"] = tfw.to_dict("records")
    pt = paired_tests(comb, PRIMARY_N, PRIMARY_A, "P", ["RMSE_x2", "RMSE_x3", "FIX_0.25min", "FIX_0.5min", "FIX_1min", "FIX_2min", "FIX_5pct", "FIX_10pct"], "winkler", B=1000)
    stats_out["fixed_window_tests"] = pt.to_dict("records")
    fwd = load("fixedwin", "decoys")
    if not fwd.empty:
        dd = fwd.groupby(["method"]).agg(decoy_elim_mean=("decoy_elim_rate", "mean"), true_retained_mean=("true_retained", "mean")).reset_index()
        md(dd, os.path.join(TAB, "tableS11_fixed_window_decoys")); metrics["fixed_window_decoys"] = dd.to_dict("records")

# ---------------- very small calibrant sets (exploratory) ----------------
sn = load("smalln", "rows")
if not sn.empty:
    snp = per_system(sn); allps.append(snp.assign(experiment="smalln"))
    rows_sn = []
    for n in sorted(snp.n.unique()):
        for al in sorted(snp.alpha.unique()):
            s_ = summary_table(snp, n, al, [m for m in ORDER if m in snp.method.unique()], B=500); s_.insert(0, "alpha", al); s_.insert(0, "n", n); rows_sn.append(s_)
    tsn = pd.concat(rows_sn, ignore_index=True); md(tsn, os.path.join(TAB, "tableS12_small_n")); metrics["small_n"] = tsn.to_dict("records")

# ---------------- NORMAN external validation (Addendum C) ----------------
nr = load("norman", "rows")
if not nr.empty:
    nr = nr.assign(system=nr.cs)
    out_nm = []
    for scen, g in nr.groupby("scenario"):
        gp = per_system(g)
        for n in sorted(gp.n.unique()) if scen != "fixed_mix" else [None]:
            gg = gp if n is None else gp[gp.n == n]
            if n is None:
                gg = gg.assign(n=0)
            for al in sorted(gg.alpha.unique()):
                s_ = summary_table(gg, 0 if n is None else n, al, [m for m in ORDER + ["P+SASC", "J+SASC"] if m in gg.method.unique()], B=1000)
                s_.insert(0, "alpha", al); s_.insert(0, "n", "lab mix" if n is None else n); s_.insert(0, "scenario", scen); out_nm.append(s_)
        tst = paired_tests(gp.assign(n=0) if scen == "fixed_mix" else gp[gp.n == 20], 0 if scen == "fixed_mix" else 20, PRIMARY_A, "P",
                           [m for m in BASELINES if m in gp.method.unique()], "winkler", B=1000)
        stats_out.setdefault("norman_tests", {})[scen] = tst.to_dict("records")
    tnm = pd.concat(out_nm, ignore_index=True); md(tnm, os.path.join(TAB, "table10_norman_external")); metrics["norman"] = tnm.to_dict("records")
    nst = load("norman", "strata")
    if not nst.empty:
        # strata are recorded at the primary calibrant number of each calibration scenario -> report per scenario
        nsc = pd.concat([strata_coverage(g_.assign(system=g_.cs), PRIMARY_A).assign(scenario=sc_) for sc_, g_ in nst.groupby("scenario")], ignore_index=True)
        nsc = nsc[["scenario"] + [c for c in nsc.columns if c != "scenario"]]
        md(nsc, os.path.join(TAB, "tableS13_norman_strata")); stats_out["norman_strata"] = nsc.to_dict("records")
    nmeta = os.path.join(a.out, "results", "experiment_logs", "norman_meta.json")
    if os.path.exists(nmeta):
        metrics["norman_meta"] = json.load(open(nmeta))

# ---------------- SASC (Addendum D) ----------------
sasc = os.path.join(a.raw, "sasc_rows.parquet")
if os.path.exists(sasc):
    sr = pd.read_parquet(sasc)
    sg = sr.groupby(["rule", "method", "system"], as_index=False).agg(coverage=("coverage", "mean"), width=("width_mean", "mean"), winkler=("winkler", "mean"))  # same width statistic as elsewhere
    rows_s = []
    for (rule, meth), g in sg.groupby(["rule", "method"]):
        lo, hi = boot_ci(g.coverage.values, B=2000)
        rows_s.append({"rule": rule, "method": meth, "systems": len(g), "coverage_mean": g.coverage.mean(), "coverage_ci_lo": lo, "coverage_ci_hi": hi,
                       "share_severe_under": float((g.coverage < 0.85).mean()), "width_median_min": g.width.median(), "winkler_mean": g.winkler.mean()})
    ts = pd.DataFrame(rows_s); md(ts, os.path.join(TAB, "table12_sasc")); metrics["sasc"] = ts.to_dict("records")
    cv = pd.read_csv(os.path.join(a.raw, "sasc_c_values.csv")); md(cv, os.path.join(TAB, "tableS15_sasc_c_values")); metrics["sasc_c"] = cv.to_dict("records")
    sp = []
    for rule in sorted(sg.rule.unique()):
        A_ = sg[(sg.rule == rule) & (sg.method == "P+SASC")].set_index("system"); B_ = sg[(sg.rule == rule) & (sg.method == "J+SASC")].set_index("system")
        j_ = A_.index.intersection(B_.index); d_ = A_.loc[j_, "winkler"] - B_.loc[j_, "winkler"]
        sp.append({"rule": rule, "comparison": "P+SASC vs J+SASC (interval score)", "systems": len(j_), "median_diff": float(d_.median()),
                   "rank_biserial": rank_biserial(d_.values), "p_wilcoxon": float(_st.wilcoxon(d_).pvalue) if (d_ != 0).any() else 1.0})
    stats_out["sasc_tests"] = sp
    ncf = os.path.join(a.raw, "sasc_norman_c.json")
    if os.path.exists(ncf):
        metrics["sasc_norman_c"] = json.load(open(ncf))

# ---------------- training-set sizes per fold (RP and HILIC) ----------------
ms_ = meta.copy(); ms_["fold"] = ms_.source.map(fold_of_source)
tsz = []
for mt_ in ["RP", "HILIC"]:
    mm = ms_[ms_.method_type == mt_]
    tgt = mm[(mm.n_compounds >= 60) & (~mm.source.str.contains("Domingo", na=False))]
    for f_ in sorted(tgt.fold.unique()):
        tr_ = mm[mm.fold != f_]
        tsz.append({"mode": mt_, "fold": int(f_), "target_systems": int((tgt.fold == f_).sum()), "training_systems_all": int(len(tr_)),
                    "training_systems_ge30": int((tr_.n_compounds >= 30).sum())})
tsz = pd.DataFrame(tsz); md(tsz, os.path.join(TAB, "tableS17_training_sizes")); metrics["training_sizes"] = tsz.to_dict("records")
# missing chromatographic metadata among RP systems (RepoRT encodes missing numeric metadata as 0; NORMAN methods enter as NaN,
# which LightGBM routes like 0 for splits learned without NaN, so both are treated identically by the difficulty model)
mrp_ = meta[meta.method_type == "RP"]
_miss = {c: (pd.to_numeric(mrp_[c], errors="coerce").fillna(0) == 0) for c in ["temp", "flow", "t0"]}
_all3 = _miss["temp"] & _miss["flow"] & _miss["t0"]
metrics["metadata_missing_rp"] = {"n_systems": int(len(mrp_)), **{c: int(v.sum()) for c, v in _miss.items()}, "all3": int(_all3.sum()),
                                  "all3_sources": sorted(set(mrp_.loc[_all3, "source"].astype(str).map(lambda s_: FOLD_NAMES[fold_of_source(s_)])))}
# per-fold size of the difficulty-model training set in the main run (parsed from the committed run log)
import re as _re
_lg = os.path.join(a.out, "results", "experiment_logs", "final_logs", "full.log.txt")
if os.path.exists(_lg):
    _rows = {int(f_): int(r_) for f_, r_ in _re.findall(r"\[main\] fold (\d) \(.*?sigma rows (\d+)", open(_lg).read())}
    if _rows:
        metrics["sigma_rows_main"] = {"per_fold": _rows, "min": min(_rows.values()), "max": max(_rows.values())}

# ---------------- review cycle 1 analyses (source-level inference, per-draw coverage, strata estimators, span model) ----------------
from conformrt.statistics import revision as RV
REV = {}
sub20_ = ps[(ps.n == PRIMARY_N) & np.isclose(ps.alpha, PRIMARY_A)]
def _lr(df_, ref, comp):
    A_ = df_[df_.method == ref].set_index("system"); B_ = df_[df_.method == comp].set_index("system"); j_ = A_.index.intersection(B_.index)
    return np.log(A_.loc[j_, "winkler"] / B_.loc[j_, "winkler"]), A_.loc[j_, "fold"]
sm_rows = []
for comp in [c for c in BASELINES if c in sub20_.method.unique()]:
    lr_, src_ = _lr(sub20_, "P", comp); r_ = RV.source_meta(lr_.values, src_.values); r_.update({"setting": "main", "ref": "P", "comparator": comp}); sm_rows.append(r_)
for abl in [m for m in ABLATIONS if m in sub20_.method.unique()]:
    lr_, src_ = _lr(sub20_, abl, "P"); r_ = RV.source_meta(lr_.values, src_.values); r_.update({"setting": "ablation", "ref": abl, "comparator": "P"}); sm_rows.append(r_)
for tag_, set_ in [("gen_smrt", "SMRT base model"), ("gen_hilic", "HILIC")]:
    pth_ = os.path.join(a.raw, f"{tag_}_rows.parquet")
    if os.path.exists(pth_):
        g_ = per_system(pd.read_parquet(pth_)); g_ = g_[(g_.n == PRIMARY_N) & np.isclose(g_.alpha, PRIMARY_A)]
        lr_, src_ = _lr(g_, "P", "B5_jackknife_plus"); r_ = RV.source_meta(lr_.values, src_.values); r_.update({"setting": set_, "ref": "P", "comparator": "B5_jackknife_plus"}); sm_rows.append(r_)
smt = pd.DataFrame([{k: v for k, v in r_.items() if k != "per_source_ratio"} | {"per_source_ratio": "; ".join(f"{FOLD_NAMES.get(int(float(k)), k)}: {v:.3f}" for k, v in r_["per_source_ratio"].items())} for r_ in sm_rows])
md(smt, os.path.join(TAB, "tableS18_source_meta")); stats_out["source_meta"] = sm_rows
# calibration-conditional (per-draw) coverage; exact split-conformal reference under exchangeability
dc = RV.draw_coverage(main_rows, ORDER); dc["split_reference_share_below_nom_minus_0.05"] = [
    RV.split_conformal_reference(int(n_), float(a_)).get("share_below_nom_minus_0.05", np.nan) if m_ == "B4_split_conformal" else np.nan
    for n_, a_, m_ in zip(dc.n, dc.alpha, dc.method)]
md(dc, os.path.join(TAB, "tableS19_draw_coverage")); REV["draw_coverage"] = dc.to_dict("records")
REV["split_reference"] = {f"n{n_}_a{a_:.2f}": RV.split_conformal_reference(n_, a_) for n_ in [10, 15, 20, 30, 40] for a_ in [0.05, 0.10, 0.20]}
# stratum coverage: system-average and pooled estimators with bootstrap CIs over systems
st_ = load("main", "strata")
if not st_.empty:
    se_ = RV.strata_estimators(st_[st_.stratum != "superclass"], PRIMARY_A, ["P", "B5_jackknife_plus", "B1_gauss_rmse", "B3_gp_PI"], B=min(a.B, 2000))
    md(se_, os.path.join(TAB, "tableS20_strata_estimators")); REV["strata_estimators"] = se_.to_dict("records")
# gradient length: continuous, source-adjusted model on RepoRT (exploratory) and its check on HILIC and NORMAN
lr_, src_ = _lr(sub20_, "P", "B5_jackknife_plus")
span_ = (meta.set_index("system").loc[lr_.index, "rt_max"] - meta.set_index("system").loc[lr_.index, "rt_min"]).astype(float)
spm = RV.span_model(lr_.values, np.log(span_.values), src_.values); REV["span_model"] = spm
ext = {}
pth_ = os.path.join(a.raw, "gen_hilic_rows.parquet")
if os.path.exists(pth_):
    g_ = per_system(pd.read_parquet(pth_)); g_ = g_[(g_.n == PRIMARY_N) & np.isclose(g_.alpha, PRIMARY_A)]
    lrh, _ = _lr(g_, "P", "B5_jackknife_plus"); sph = (meta.set_index("system").loc[lrh.index, "rt_max"] - meta.set_index("system").loc[lrh.index, "rt_min"]).astype(float)
    ext["HILIC"] = RV.external_span_check(lrh.values, np.log(sph.values), spm)
nr_ = load("norman", "rows")
if not nr_.empty:
    for scen_ in ["random_all", "fixed_mix", "pool_random"]:
        q_ = nr_[(nr_.scenario == scen_) & np.isclose(nr_.alpha, PRIMARY_A) & nr_.method.isin(["P", "B5_jackknife_plus"])]
        q_ = q_[q_.n.astype(str) == ("20" if scen_ != "fixed_mix" else q_.n.astype(str).iloc[0])] if len(q_) else q_
        if q_.empty:
            continue
        w_ = q_.groupby(["cs", "method"]).winkler.mean().unstack()
        spn = (q_.width_mean / q_.width_rel_mean).groupby(q_.cs).median()
        j_ = w_.dropna().index
        ext[f"NORMAN {scen_}"] = RV.external_span_check(np.log(w_.loc[j_, "P"] / w_.loc[j_, "B5_jackknife_plus"]).values, np.log(spn.loc[j_].values), spm)
REV["span_external"] = ext
md(pd.DataFrame([{"set": k_, **v_} for k_, v_ in ext.items()]), os.path.join(TAB, "tableS21_span_external"))
# SASC per held-out source group
sasc_ = os.path.join(a.raw, "sasc_rows.parquet")
if os.path.exists(sasc_):
    sr_ = pd.read_parquet(sasc_).groupby(["rule", "method", "fold", "system"], as_index=False).coverage.mean()
    sps = sr_.groupby(["rule", "method", "fold"]).agg(systems=("system", "nunique"), coverage=("coverage", "mean"),
                                                       share_below_085=("coverage", lambda x: float((x < 0.85).mean()))).reset_index()
    sps["fold_name"] = sps.fold.map(lambda f_: FOLD_NAMES.get(int(f_), f_))
    md(sps, os.path.join(TAB, "tableS22_sasc_per_source")); REV["sasc_per_source"] = sps.to_dict("records")
# HILIC per fold and relative widths
if os.path.exists(os.path.join(a.raw, "gen_hilic_rows.parquet")):
    g_ = per_system(pd.read_parquet(os.path.join(a.raw, "gen_hilic_rows.parquet"))); g_ = g_[(g_.n == PRIMARY_N) & np.isclose(g_.alpha, PRIMARY_A)]
    hpf = g_.groupby(["fold", "method"]).agg(systems=("system", "nunique"), coverage=("coverage", "mean"), width_median=("width", "median"),
                                             width_rel_median=("width_rel", "median")).reset_index()
    hpf["fold_name"] = hpf.fold.map(lambda f_: FOLD_NAMES.get(int(f_), f_))
    md(hpf, os.path.join(TAB, "tableS23_hilic_per_fold")); REV["hilic_per_fold"] = hpf.to_dict("records")
# windows with a distribution-free guarantee of at least 0.90: jackknife+-type at alpha = 0.05 versus split conformal at 0.10
gw = []
for n_ in [20, 30, 40]:
    for m_, a_ in [("P", 0.05), ("B5_jackknife_plus", 0.05), ("B4_split_conformal", 0.10), ("P", 0.10), ("B5_jackknife_plus", 0.10)]:
        g_ = ps[(ps.n == n_) & np.isclose(ps.alpha, a_) & (ps.method == m_)]
        if len(g_):
            gw.append({"n": n_, "method": m_, "alpha": a_, "guaranteed_level": (1 - 2 * a_) if m_ != "B4_split_conformal" else (1 - a_),
                       "coverage": float(g_.coverage.mean()), "width_median": float(g_.width.median())})
gw = pd.DataFrame(gw); md(gw, os.path.join(TAB, "tableS24_guaranteed_windows")); REV["guaranteed_windows"] = gw.to_dict("records")
stats_out["revision1"] = REV

# ---------------- machine-readable outputs ----------------
pd.concat(allps, ignore_index=True).to_csv(os.path.join(RES, "raw_results.csv"), index=False)
def _clean(o):
    if isinstance(o, dict): return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, list): return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)): return None if (np.isnan(o) or np.isinf(o)) else float(o)
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.bool_,)): return bool(o)
    return o
json.dump(_clean(metrics), open(os.path.join(RES, "metrics.json"), "w"), indent=1)
json.dump(_clean(stats_out), open(os.path.join(RES, "statistics.json"), "w"), indent=1)
F.fig_workflow(FIG)
print("analysis complete:", len(os.listdir(TAB)), "table files;", len(os.listdir(FIG)), "figure files")

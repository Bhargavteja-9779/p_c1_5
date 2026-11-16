"""Generate all tables, figures and statistics from raw experiment outputs (no hand-typed numbers).
Usage: python experiments/analyze.py --raw ../results/raw --out ../ [--B 10000]"""
import argparse, os, sys, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from conformrt.statistics.analysis import (per_system, summary_table, paired_tests, strata_coverage, boot_ci, rank_biserial,
                                           METHOD_LABELS, BASELINES, ABLATIONS)
from conformrt.visualization import figures as F

ap = argparse.ArgumentParser(); ap.add_argument("--raw", default="../results/raw"); ap.add_argument("--out", default="..")
ap.add_argument("--processed", default="data/processed"); ap.add_argument("--B", type=int, default=10000)
a = ap.parse_args()
RES, TAB, FIG = (os.path.join(a.out, d) for d in ("results", "tables", "figures"))
for d in (RES, TAB, FIG):
    os.makedirs(d, exist_ok=True)
PRIMARY_N, PRIMARY_A = 20, 0.1
ORDER = ["P", "B6_jk+_ensemble", "B5_jackknife_plus", "B4_split_conformal", "B3_gp_PI", "B2_pspline_PI", "B1_gauss_rmse"]
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
if main_rows.empty:
    sys.exit("no main results yet")
ps = per_system(main_rows)
allps = [ps.assign(experiment="main", base="order", strategy="random", perturb="none")]

# ---------------- Table 1: dataset ----------------
meta = pd.read_csv(os.path.join(a.processed, "systems.csv"), dtype={"system": str})
from conformrt.evaluation.protocol import fold_of_source, FOLD_NAMES
tsys = ps.system.unique()
mt = meta[meta.system.isin(tsys)].copy(); mt["fold"] = mt.source.map(fold_of_source); mt["fold_name"] = mt.fold.map(FOLD_NAMES)
t1 = mt.groupby("fold_name").agg(systems=("system", "size"), compounds_median=("n_compounds", "median"),
                                  compounds_total=("n_compounds", "sum"), rt_span_median_min=("rt_max", "median"),
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
F.fig_alpha_calibration(full, FIG, [m for m in ["P", "B5_jackknife_plus", "B4_split_conformal", "B3_gp_PI", "B2_pspline_PI", "B1_gauss_rmse"] if m in ps.method.unique()], L, n=PRIMARY_N)

# coverage validity (system-bootstrap and source-cluster bootstrap) for every n/alpha
val = []
for (n, al, m), s in ps.groupby(["n", "alpha", "method"]):
    lo, hi = boot_ci(s.coverage.values, B=min(a.B, 4000))
    clo, chi = boot_ci(s.coverage.values, B=min(a.B, 4000), clusters=s.fold.values)
    val.append({"n": n, "alpha": al, "method": m, "mean_cov": s.coverage.mean(), "ci_lo": lo, "ci_hi": hi, "cluster_ci_lo": clo,
                "cluster_ci_hi": chi, "nominal": 1 - al, "mean_minus_nominal": s.coverage.mean() - (1 - al),
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
    F.fig_strata(sc, FIG, [m for m in ["P", "B5_jackknife_plus", "B1_gauss_rmse"] if m in sc.method.unique()], L, PRIMARY_A, order)
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
                 [m for m in ["P", "B5_jackknife_plus", "B3_gp_PI", "B2_pspline_PI", "B1_gauss_rmse"] if m in t7.method.unique()], L)
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
    F.fig_robustness(rob, FIG, [m for m in ["P", "B5_jackknife_plus", "B3_gp_PI", "B1_gauss_rmse"] if m in rob.method.unique()], L)

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
                s_ = summary_table(gg, 0 if n is None else n, al, [m for m in ORDER if m in gg.method.unique()], B=1000)
                s_.insert(0, "alpha", al); s_.insert(0, "n", "lab mix" if n is None else n); s_.insert(0, "scenario", scen); out_nm.append(s_)
        tst = paired_tests(gp.assign(n=0) if scen == "fixed_mix" else gp[gp.n == 20], 0 if scen == "fixed_mix" else 20, PRIMARY_A, "P",
                           [m for m in BASELINES if m in gp.method.unique()], "winkler", B=1000)
        stats_out.setdefault("norman_tests", {})[scen] = tst.to_dict("records")
    tnm = pd.concat(out_nm, ignore_index=True); md(tnm, os.path.join(TAB, "table10_norman_external")); metrics["norman"] = tnm.to_dict("records")
    nst = load("norman", "strata")
    if not nst.empty:
        nsc = strata_coverage(nst.assign(system=nst.cs), PRIMARY_A); md(nsc, os.path.join(TAB, "tableS13_norman_strata")); stats_out["norman_strata"] = nsc.to_dict("records")
    nmeta = os.path.join(a.out, "results", "experiment_logs", "norman_meta.json")
    if os.path.exists(nmeta):
        metrics["norman_meta"] = json.load(open(nmeta))

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

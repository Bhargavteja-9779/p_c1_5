"""Build the Supplementary Material (Sections S1-S6, Tables S1-S26, Figs. S1-S2) from machine-generated result files.
Every number comes from tables/*.csv or results/*.json; text tokens {{key|fmt}} are resolved as in build_manuscript.py
and unresolved tokens abort the build. Supplementary tables are numbered in the order in which they appear in this
document; supplementary/table_index.csv maps each number to the analysis output it is built from.
Outputs: supplementary/Supplementary_Material.{md,docx,pdf}, supplementary/tables/Table_S<k>_<name>.csv, supplementary/table_index.csv
Usage: python scripts/build_supplementary.py [--no-pdf]"""
import argparse, json, os, re, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "src"))
from docbuild import resolve, pandoc_docx, postprocess_docx, docx_to_pdf, mdtab  # noqa: E402
from manuscript_facts import build_facts  # noqa: E402
from conformrt.statistics.analysis import METHOD_LABELS  # noqa: E402
from conformrt.evaluation.protocol import FOLD_NAMES  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--no-pdf", action="store_true"); a = ap.parse_args()
RES, TAB, FIG, SUP = (os.path.join(ROOT, d) for d in ("results", "tables", "figures", "supplementary"))
OUT_TAB = os.path.join(SUP, "tables"); os.makedirs(OUT_TAB, exist_ok=True)

LAB = dict(METHOD_LABELS)
LAB.update({"P+SASC": "Normalised jackknife+ + SASC", "J+SASC": "Jackknife+ + SASC", "J+": "Jackknife+",
            "RMSE_x2": "±2·RMSE", "RMSE_x3": "±3·RMSE", "FIX_0.25min": "±0.25 min", "FIX_0.5min": "±0.5 min", "FIX_1min": "±1 min",
            "FIX_2min": "±2 min", "FIX_5pct": "±5 % of predicted RT", "FIX_10pct": "±10 % of predicted RT"})
MAIN7 = ["P", "B5_jackknife_plus", "B6_jk+_ensemble", "B4_split_conformal", "B1_gauss_rmse", "B2_pspline_PI", "B3_gp_PI", "B7_gp_meta"]
SCEN = {"fixed_mix": "Laboratory calibrant mix", "pool_random": "Random calibration chemicals", "random_all": "Random compounds"}
STRAT = {("seen", "seen"): "Seen by base model", ("seen", "unseen"): "Not seen by base model",
         ("similarity", "sim>=0.7"): "Tanimoto ≥ 0.7", ("similarity", "sim0.4-0.7"): "Tanimoto 0.4–0.7", ("similarity", "sim<0.4"): "Tanimoto < 0.4",
         ("ionisable", "neutral"): "Neutral", ("ionisable", "ionisable"): "Ionisable",
         ("region", "inside_cal_range"): "Inside calibrant score range", ("region", "outside_cal_range"): "Outside calibrant score range",
         ("rt_tertile", "early"): "First third of RT range", ("rt_tertile", "middle"): "Middle third of RT range", ("rt_tertile", "late"): "Last third of RT range",
         ("sample_type", "cal"): "Calibration chemicals", ("sample_type", "sus"): "Suspects"}
SETTING = {"sens_l15_e400": "15 leaves, 400 trees", "sens_l63_e400": "63 leaves, 400 trees", "sens_l31_e200": "31 leaves, 200 trees",
           "sens_l31_e800": "31 leaves, 800 trees"}
MECH = {"main": "Random calibrants (main)", "robust_span_none": "RT-spanning calibrants", "robust_common_none": "Commonly measured calibrants",
        "robust_random_noise": "Calibrant RT noise", "robust_random_outlier": "One outlying calibrant"}


# ------------------------------------------------------------------ formatting helpers
def num(x, d=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    x = float(x)
    return "∞" if np.isinf(x) else f"{x:.{d}f}"


def pct(x, d=0):
    return "–" if x is None or pd.isna(x) else f"{100 * float(x):.{d}f}"


def pval(x):
    return "–" if pd.isna(x) else ("< 0.001" if x < 0.001 else f"{x:.3f}")


def ci(m, lo, hi, d=3):
    return f"{num(m, d)} [{num(lo, d)}, {num(hi, d)}]"


def lab(m):
    return LAB.get(m, m)


SUMMARY_HEAD = ["Method", "Systems", "Mean coverage [95 % CI]", "10th percentile of coverage", "Systems > 0.05 below nominal (%)",
                "Median width (min)", "Mean interval score (min)"]


def summary_cells(r):
    return [lab(r.method), int(r.n_systems), ci(r.coverage_mean, r.coverage_ci_lo, r.coverage_ci_hi), num(r.coverage_p10),
            pct(r.share_severe_under), num(r.width_median_min, 2), num(r.winkler_mean, 2)]


def rd(name):
    return pd.read_csv(os.path.join(TAB, name + ".csv"))


# ------------------------------------------------------------------ tables (in document order)
TABLES = []  # (number, source, slug, caption, markdown, dataframe for csv export)


def add(source, slug, caption, header, rows, df):
    k = len(TABLES) + 1
    TABLES.append((k, source, slug, caption, mdtab(header, rows), df))
    return k


def t_training_sizes():
    d = rd("tableS17_training_sizes")
    rows = [[r["mode"], FOLD_NAMES.get(int(r.fold), r.fold), int(r.target_systems), int(r.training_systems_all), int(r.training_systems_ge30)] for _, r in d.iterrows()]
    add("tableS17_training_sizes", "training_sizes",
        "Target systems and training systems per leave-source-out fold. Training systems are all systems of the same separation mode in the "
        "other four source groups; those with at least 30 compounds (last column) were used to simulate calibrant protocols for training the "
        "difficulty model. METLIN SMRT belongs to the 'other sources' group; a random subsample of 15,000 of its retention times was added to "
        "the base-model training data in every fold.",
        ["Mode", "Held-out source group", "Target systems", "Training systems", "Training systems with ≥ 30 compounds"], rows, d)


def t_all_n_alpha():
    d = rd("tableS1_all_n_alpha"); d = d[d.method.isin(MAIN7)]
    rows = [[int(r.n), num(r.alpha, 2)] + summary_cells(r) for _, r in d.iterrows()]
    add("tableS1_all_n_alpha", "all_n_alpha",
        "Coverage, width and interval score of all window constructions for every calibrant number n and level α (RP target systems; "
        "per-system averages over 20 calibrant draws, then averaged over systems; 95 % bootstrap CI over systems). ∞: unbounded window "
        "(too few calibrants for the requested level; coverage is then trivially 1). The ablated variants are given in Table S17.",
        ["n", "α"] + SUMMARY_HEAD, rows, d)


def t_validity():
    d = rd("tableS2_coverage_validity"); d = d[d.method.isin(MAIN7)]
    rows = [[int(r.n), num(r.alpha, 2), lab(r.method), num(r.mean_cov), f"[{num(r.ci_lo)}, {num(r.ci_hi)}]", f"[{num(r.cluster_ci_lo)}, {num(r.cluster_ci_hi)}]",
             f"{float(r.mean_minus_nominal):+.3f}", pct(r["share_systems_below_nominal_minus_0.05"])] for _, r in d.iterrows()]
    add("tableS2_coverage_validity", "coverage_validity",
        "Mean coverage with 95 % bootstrap confidence intervals from resampling systems and from a cluster bootstrap over the five source "
        "folds, deviation from the nominal level 1 − α and share of systems whose coverage is more than 0.05 below nominal. A coverage of 1.000 "
        "with interval [1.000, 1.000] indicates unbounded windows (Table S2).",
        ["n", "α", "Method", "Mean coverage", "95 % CI (systems)", "95 % CI (cluster bootstrap)", "Mean − nominal", "Systems > 0.05 below nominal (%)"], rows, d)


def t_paired_all():
    d = rd("tableS3_paired_winkler_all_n_alpha")
    rows = [[int(r.n), num(r.alpha, 2), lab(r.comparator), int(r.n_systems), num(r.median_ref, 2), num(r.median_comp, 2),
             f"{float(r.median_diff):+.2f} [{float(r.diff_ci_lo):+.2f}, {float(r.diff_ci_hi):+.2f}]", pct(r.share_ref_better),
             f"{float(r.rank_biserial):+.2f}", pval(r.p_holm)] for _, r in d.iterrows()]
    add("tableS3_paired_winkler_all_n_alpha", "paired_interval_score_all_n_alpha",
        "Paired comparison of per-system interval scores of the normalised jackknife+ windows with each comparator for every n and α (two-sided "
        "Wilcoxon signed-rank tests with Holm correction; negative differences favour the normalised jackknife+ windows). Comparators whose windows "
        "are unbounded at a given (n, α) are omitted.",
        ["n", "α", "Comparator", "Systems", "Median score, normalised jackknife+ (min)", "Median score, comparator (min)", "Median paired difference [95 % CI] (min)",
         "Systems where normalised jackknife+ is better (%)", "Rank-biserial r", "Holm-adjusted p"], rows, d)


def t_small_n():
    d = rd("tableS12_small_n"); d = d[d.method.isin(MAIN7)]
    rows = [[int(r.n), num(r.alpha, 2)] + summary_cells(r) for _, r in d.iterrows()]
    add("tableS12_small_n", "small_n",
        "Exploratory results with very small calibrant sets (n = 5 and 8; post hoc). ∞: unbounded window (n < 1/α − 1).",
        ["n", "α"] + SUMMARY_HEAD, rows, d)


def t_fixed():
    d = rd("table9_fixed_windows")
    rows = [[int(r.n)] + summary_cells(r) for _, r in d.iterrows()]
    add("table9_fixed_windows", "fixed_tolerance_windows",
        "Fixed-tolerance windows (±0.25 to ±2 min, ±5 % and ±10 % of the predicted RT) and residual-based windows (±2·RMSE, ±3·RMSE of the "
        "calibrant residuals) evaluated on the same calibrant draws as the main study (post hoc), with the normalised jackknife+ windows, jackknife+ and the "
        "Gaussian window (α = 0.10) for reference. The column 'Systems > 0.05 below nominal' refers to the nominal level 0.90.",
        ["n"] + SUMMARY_HEAD, rows, d)


def t_per_fold():
    d = rd("tableS4_per_fold_n20_a010"); d = d[d.method.isin(MAIN7)]
    rows = [[r.fold_name, lab(r.method), int(r.systems), num(r.coverage), num(r.winkler, 2), num(r.width_median, 2)] for _, r in d.iterrows()]
    add("tableS4_per_fold_n20_a010", "per_fold",
        "Mean coverage, mean interval score and median width per held-out source group (n = 20, α = 0.10).",
        ["Source group", "Method", "Systems", "Mean coverage", "Mean interval score (min)", "Median width (min)"], rows, d)


def t_fold_het():
    d = rd("table11_fold_heterogeneity")
    rows = [[lab(r.comparator), "All folds" if r.fold_name == "all" else r.fold_name, int(r.systems), f"{float(r.median_diff_min):+.2f}",
             f"{100 * float(r.median_rel_change):+.0f}", pct(r.share_P_better), f"{float(r.rank_biserial):+.2f}", pval(r.p_wilcoxon)] for _, r in d.iterrows()]
    add("table11_fold_heterogeneity", "per_fold_paired_tests",
        "Paired comparison of per-system interval scores of the normalised jackknife+ windows with three comparators within each held-out source group "
        "(n = 20, α = 0.10; negative differences favour the normalised jackknife+ windows; unadjusted two-sided Wilcoxon signed-rank tests).",
        ["Comparator", "Source group", "Systems", "Median paired difference (min)", "Median relative change (%)", "Systems where normalised jackknife+ is better (%)",
         "Rank-biserial r", "p"], rows, d)


def strata_pivot(d, methods, extra_index=()):
    d = d[d.method.isin(methods) & (d.stratum != "superclass")].copy()
    d["Stratum"] = [STRAT.get((s, g), f"{s}: {g}") for s, g in zip(d.stratum, d.group)]
    idx = list(extra_index) + ["stratum", "group", "Stratum"]
    pv = d.pivot_table(index=idx, columns="method", values="coverage_sysavg", sort=False)
    ev = d.groupby(idx, sort=False).total.max()  # evaluations per method are identical within a stratum
    return pv, ev


def t_strata():
    d = rd("tableS5_strata_n20_a010"); pv, ev = strata_pivot(d, MAIN7)
    meths = [m for m in MAIN7 if m in pv.columns]
    rows = [[k[-1], f"{int(ev.loc[k]):,}"] + [num(pv.loc[k, m]) for m in meths] for k in pv.index]
    add("tableS5_strata_n20_a010", "coverage_by_stratum",
        "Mean coverage (per-system coverage within the stratum, averaged over systems) by stratum of test compounds (n = 20, α = 0.10, 20 draws). "
        "Seen/not seen: test compound present in/absent from the base-model training data; Tanimoto: maximum similarity to the training "
        "compounds; calibrant score range: range of the base-model scores of the calibrants; RT range thirds: thirds of each system's retention-time range "
        "(not thirds of the compounds).",
        ["Stratum", "Compound–draw evaluations"] + [lab(m) for m in meths], rows, d)


def t_superclass():
    d = rd("tableS6_superclass_coverage"); meths = [m for m in MAIN7 if m in d.columns]
    rows = [[r.group] + [num(r[m]) for m in meths] for _, r in d.iterrows()]
    add("tableS6_superclass_coverage", "coverage_by_superclass",
        "Mean coverage by ClassyFire superclass (n = 20, α = 0.10).", ["Superclass"] + [lab(m) for m in meths], rows, d)


def t_worst():
    d = rd("tableS7_worst_stratum"); d = d[d.method.isin(MAIN7)]
    rows = [[lab(r.method), num(r.worst_stratum_coverage)] for _, r in d.iterrows()]
    add("tableS7_worst_stratum", "worst_stratum",
        "Lowest stratum-level coverage of each method over the strata of Table S9 (n = 20, α = 0.10; nominal 0.90).",
        ["Method", "Lowest stratum coverage"], rows, d)


def t_errors():
    d = rd("tableS9_error_classes")
    rows = [[r.cls, int(r.misses), pct(r.miss_rate, 1)] for _, r in d.iterrows()]
    add("tableS9_error_classes", "error_classes",
        "ClassyFire classes with the largest number of misses of the normalised jackknife+ 90 % windows (n = 20, first calibrant draw of each system); "
        "miss rate = misses divided by the number of test compounds of the class.", ["ClassyFire class", "Misses", "Miss rate (%)"], rows, d)


def t_decoys():
    d = rd("table7_candidate_filtering")
    rows = [[lab(r.method), num(r.alpha, 2), int(r.systems), pct(r.true_retained_mean, 1), pct(r.decoy_elim_mean, 1), pct(r.decoy_elim_median, 1),
             f"{r.parents_per_system:.0f}", f"{r.decoys_per_system:.0f}"] for _, r in d.iterrows()]
    add("table7_candidate_filtering", "candidate_filtering",
        "Candidate filtering with isomeric decoys (same molecular formula, different InChIKey connectivity block; n = 20, calibrant draws 0–4): "
        "mean retention of the correct structure, mean and median share of decoys eliminated, and mean numbers of test compounds with at least "
        "one decoy and of decoys per system and draw.",
        ["Method", "α", "Systems", "Correct structure retained (%)", "Decoys eliminated, mean (%)", "Decoys eliminated, median (%)",
         "Test compounds with decoys", "Decoys"], rows, d)


def t_decoy_tests():
    d = rd("tableS8_decoy_tests")
    rows = [[num(r.alpha, 2), lab(r.comparator), int(r.n_systems), f"{100 * float(r.median_diff_elim):+.1f}", f"{float(r.rank_biserial):+.2f}", pval(r.p_wilcoxon)]
            for _, r in d.iterrows()]
    add("tableS8_decoy_tests", "decoy_paired_tests",
        "Paired comparison of per-system decoy elimination of the normalised jackknife+ windows with comparators at the same α (positive differences: the "
        "normalised jackknife+ windows eliminate more decoys; unadjusted two-sided Wilcoxon signed-rank tests). Differences at the same nominal level partly "
        "reflect different retention of the correct structure (Table S13); see Table S15 for a comparison at matched retention.",
        ["α", "Comparator", "Systems", "Median paired difference (percentage points)", "Rank-biserial r", "p"], rows, d)


def t_matched():
    d = rd("tableS14_matched_retention")
    rows = [[r.label, pct(r.retention), pct(r.decoy_elim_interp, 1)] for _, r in d.iterrows()]
    add("tableS14_matched_retention", "matched_retention",
        "Share of isomeric decoys eliminated at matched retention of the correct structure (85 % and 90 %), obtained by linear interpolation of "
        "each method's retention–elimination curve over α = 0.05, 0.10 and 0.20; '–': retention level not reached at any tested α.",
        ["Method", "Retention of correct structure (%)", "Decoys eliminated (%)"], rows, d)


def t_fixed_decoys():
    d = rd("tableS11_fixed_window_decoys")
    rows = [[lab(r.method), pct(r.true_retained_mean, 1), pct(r.decoy_elim_mean, 1)] for _, r in d.iterrows()]
    add("tableS11_fixed_window_decoys", "fixed_window_decoys",
        "Retention of the correct structure and elimination of isomeric decoys by fixed-tolerance and residual-based windows (n = 20, "
        "calibrant draws 0–4; post hoc).", ["Window", "Correct structure retained (%)", "Decoys eliminated (%)"], rows, d)


def t_ablation():
    d = rd("table4_ablation_n20_a010")
    rows = [[r.label, f"{float(r.median_diff):+.2f} [{float(r.diff_ci_lo):+.2f}, {float(r.diff_ci_hi):+.2f}]", f"{float(r.rank_biserial):+.2f}",
             pval(r.p_holm), num(r.coverage_mean)] for _, r in d.iterrows()]
    add("table4_ablation_n20_a010", "ablation",
        "Ablation of the difficulty model and of the construction (n = 20, α = 0.10): median paired change in the per-system interval score "
        "relative to the full normalised jackknife+ model (positive = worse) with 95 % bootstrap CI, rank-biserial correlation, Holm-adjusted p of the "
        "Wilcoxon signed-rank test, and mean coverage. Feature-group ablations retrain the difficulty model without that group.",
        ["Variant", "Median paired change in interval score [95 % CI] (min)", "Rank-biserial r", "Holm-adjusted p", "Mean coverage"], rows, d)


def t_sensitivity():
    d = rd("tableS10_sensitivity")
    rows = [[SETTING.get(r.setting, r.setting)] + summary_cells(r) for _, r in d.iterrows()]
    add("tableS10_sensitivity", "hyperparameter_sensitivity",
        "Sensitivity to the hyper-parameters of the gradient-boosted difficulty model (default: 31 leaves, 400 trees; n = 20, α = 0.10, "
        "10 calibrant draws per system).", ["Difficulty model"] + SUMMARY_HEAD, rows, d)


def t_robustness():
    d = rd("table5_robustness")
    rows = [[r.scenario[0].upper() + r.scenario[1:]] + summary_cells(r) for _, r in d.iterrows()]
    add("table5_robustness", "robustness",
        "Robustness to calibrant selection and calibrant errors (n = 20, α = 0.10; 10 calibrant draws per system in the perturbed scenarios). "
        "Noise: Gaussian noise with a standard deviation of 2 % of the RT span added to all calibrant RTs; outlier: one calibrant shifted by "
        "10 % of the span.", ["Scenario"] + SUMMARY_HEAD, rows, d)


def t_mechanism():
    d = rd("tableS16_robustness_mechanism")
    rows = [[MECH.get(r.scenario, r.scenario), lab(r.method), pct(r.share_outside_cal_range, 1), num(r.cov_inside), num(r.cov_outside), num(r.cov_late)]
            for _, r in d.iterrows()]
    add("tableS16_robustness_mechanism", "robustness_mechanism",
        "Where coverage is lost under non-random calibrant selection and calibrant errors (n = 20, α = 0.10): share of test compounds whose "
        "base score lies outside the calibrant score range, and coverage inside and outside that range and for compounds eluting in the last third of the "
        "system's RT range.",
        ["Scenario", "Method", "Test compounds outside calibrant range (%)", "Coverage inside range", "Coverage outside range", "Coverage, last third of RT range"], rows, d)


def t_generalisation():
    d = rd("table6_generalisation"); d = d[d.method.isin(MAIN7)]
    names = {"SMRT-regression base model (RP)": "SMRT-trained regressor (RP)", "HILIC systems (order model)": "HILIC systems"}
    rows = [[names.get(r.setting, r.setting), int(r.n)] + summary_cells(r) for _, r in d.iterrows()]
    add("table6_generalisation", "generalisation",
        "Generalisation to a second base model (gradient-boosted regressor trained on METLIN SMRT only, applied to the RP target systems) and to "
        "HILIC systems (retention-order and difficulty models trained on HILIC systems of the other source groups); α = 0.10.",
        ["Setting", "n"] + SUMMARY_HEAD, rows, d)


def t_norman():
    d = rd("table10_norman_external")
    rows = [[SCEN.get(r.scenario, r.scenario), "lab mix" if str(r.n) == "lab mix" else int(float(r.n)), num(r.alpha, 2)] + summary_cells(r) for _, r in d.iterrows()]
    add("table10_norman_external", "norman_external",
        "External validation on the NORMAN interlaboratory NTS data ({{norman.n_systems|int}} chromatographic methods). Laboratory calibrant mix: "
        "all detected calibration chemicals of each method as calibrants, suspects as test compounds; random calibration chemicals: random subsets "
        "of the calibration chemicals as calibrants, suspects as test compounds; random compounds: calibrants drawn at random from all compounds. "
        "SASC rows apply the selection-aware factor estimated by replaying the calibrant-mix rule on RepoRT systems (post hoc).",
        ["Calibration", "n", "α"] + SUMMARY_HEAD, rows, d)


def t_norman_strata():
    d = rd("tableS13_norman_strata")
    if "scenario" not in d.columns:
        raise SystemExit("tableS13_norman_strata.csv lacks the scenario column: rerun experiments/analyze.py")
    meths = [m for m in MAIN7 if m in set(d.method)]
    pv, ev = strata_pivot(d, meths, extra_index=("scenario",))
    rows = [[SCEN.get(k[0], k[0]), k[-1], f"{int(ev.loc[k]):,}"] + [num(pv.loc[k, m]) for m in meths if m in pv.columns] for k in pv.index]
    add("tableS13_norman_strata", "norman_coverage_by_stratum",
        "Mean coverage by stratum of test compounds on the NORMAN data (α = 0.10) for each calibration scenario (laboratory calibrant mix: all "
        "detected calibration chemicals; random scenarios: n = 20).",
        ["Calibration", "Stratum", "Compound–draw evaluations"] + [lab(m) for m in meths if m in pv.columns], rows, d)


def t_sasc():
    d = rd("table12_sasc")
    rule = {"common": "Commonly measured compounds", "random": "Random", "span": "RT-spanning"}
    rows = [[rule.get(r.rule, r.rule), lab(r.method), int(r.systems), ci(r.coverage_mean, r.coverage_ci_lo, r.coverage_ci_hi), pct(r.share_severe_under),
             num(r.width_median_min, 2), num(r.winkler_mean, 2)] for _, r in d.iterrows()]
    add("table12_sasc", "sasc",
        "Selection-aware system-level calibration (SASC) on the RepoRT target systems for three calibrant-selection rules (n = 20, α = 0.10; post hoc). "
        "The factor c is estimated per fold by replaying the rule on the training systems (Table S25).",
        ["Calibrant selection rule", "Method", "Systems", "Mean coverage [95 % CI]", "Systems > 0.05 below nominal (%)", "Median width (min)",
         "Mean interval score (min)"], rows, d)


def t_sasc_c():
    d = rd("tableS15_sasc_c_values")
    rows = [[FOLD_NAMES.get(int(r.fold), r.fold), r.rule, int(r.n_train_systems), num(r.c_P, 2), num(r.c_J, 2)] for _, r in d.iterrows()]
    nc = json.load(open(os.path.join(RES, "raw", "sasc_norman_c.json"))) if os.path.exists(os.path.join(RES, "raw", "sasc_norman_c.json")) else None
    if nc:
        rows.append(["NORMAN laboratory mixes", "NORMAN calibration chemicals", int(nc["n_systems"]), num(nc["c_P"], 2), num(nc["c_J"], 2)])
    add("tableS15_sasc_c_values", "sasc_factors",
        "SASC factors c by which the scaled scores are multiplied, estimated on the training systems of each fold for the normalised jackknife+ windows (c_P) "
        "and jackknife+ (c_J); c < 1 narrows and c > 1 widens the windows. Last row: factor applied to the NORMAN laboratory mixes, estimated on "
        "the RepoRT RP systems that contain at least 10 of the NORMAN calibration chemicals and at least 10 other compounds, with these "
        "chemicals as calibrants (each system scored with the models of its own leave-source-out fold).",
        ["Held-out source group", "Rule", "Training systems", "c_P", "c_J"], rows, d)


def t_compute():
    cb = json.load(open(os.path.join(RES, "experiment_logs", "compute_benchmark.json")))
    rows = [[k, f"{v['calibration_ms']:.0f}", f"{v['windows_us_per_candidate']:.0f}"] for k, v in cb["per_system"].items()]
    d = pd.DataFrame([{"n": int(k), **v} for k, v in cb["per_system"].items()])
    add("compute_benchmark.json", "compute",
        f"Computation time on one CPU thread (benchmark system {cb['benchmark_system']}): calibration of a new system (leave-one-out monotone maps "
        f"and difficulty scales of the calibrants) and window computation per candidate structure. One-off training on public data: base model "
        f"{cb['base_train_s_1thread'] / 60:.1f} min ({cb['base_train_rows']:,} retention times), difficulty model {cb['sigma_train_s_1thread']:.1f} s; "
        f"peak memory {cb['peak_rss_MB']:.0f} MB.", ["Calibrants n", "Calibration (ms)", "Window per candidate (µs)"], rows, d)


# ------------------------------------------------------------------ revision: cross-system conformal calibration (XSC)
XL = {"B1_gauss_rmse": "Gaussian ±z·RMSE", "B5_jackknife_plus": "Jackknife+", "P": "Jackknife+, meta-learned σ", "JP_lab": "Jackknife+ × XSC multiplier (B8)",
      "XSC_flat_lab": "XSC, σ ≡ 1 (B9)", "XSC_sigs_avg": "XSC, average target", "XSC_sigs_lab": "XSC, laboratory-level (proposed)",
      "XSC_flat_avg": "XSC, σ ≡ 1, average target", "XSC_sig_lab": "XSC, σ without shrinkage", "XSC_sig_avg": "XSC, σ without shrinkage, average"}
RULE = {"random": "Random", "span": "RT-spanning", "common": "Most common", "polar": "Most polar", "norman": "NORMAN mix"}


def t_xsc_rules():
    d = rd("xsc_table_rules_n")
    rows = [[RULE.get(r.rule, r.rule), int(r.n), XL.get(r.method, r.method), num(r.coverage_mean), pct(r.fail, 1), pct(r.fail_085, 1),
             pct(r.share_infinite), num(r.width_median, 2), num(r.is_median, 2)] for _, r in d.iterrows()]
    add("xsc_table_rules_n", "xsc_rules_n",
        "Laboratory-level retention for every calibrant-selection rule and calibrant number (target retention 0.90, δ = 0.10; 228 held-out "
        "reversed-phase systems, leave-source-out; 20 calibrant draws per system for random calibrants, 10 otherwise). Failing: share of "
        "laboratories (system × draw) with retention below 0.90 (or 0.85). Jackknife+-based windows are not defined for n < 9.",
        ["Rule", "n", "Method", "Mean retention", "Failing < 0.90 (%)", "Failing < 0.85 (%)", "Infinite windows (%)", "Median width (min)",
         "Median interval score (min)"], rows, d)


def t_xsc_mult():
    d = rd("xsc_table_multipliers")
    fam = {"flat": "σ ≡ 1", "sigs": "shrunk σ (proposed)", "jp": "jackknife+ scores"}
    rows = [[FOLD_NAMES.get(int(r.fold), r.fold), RULE.get(r.rule, r.rule), int(r.n), fam.get(r.family, r.family), int(r.S), num(r["c_avg_0.1"], 2),
             num(r["c_lab_0.1"], 2), num(r["c_labx_0.1"], 2)] for _, r in d.iterrows()]
    add("xsc_table_multipliers", "xsc_multipliers",
        "XSC multipliers per held-out source group, selection rule, calibrant number and window family, calibrated on S training systems: "
        "average target (α = 0.10), laboratory-level target (ε = δ = 0.10; draw-averaged distribution function, used in the main text) and the "
        "exact one-replay-per-system order statistic.",
        ["Held-out source", "Rule", "n", "Window family", "S", "c (average)", "c (laboratory)", "c (laboratory, exact)"], rows, d)


def t_conventional():
    d = rd("table2_main_n20_a010")
    rows = [summary_cells(r) for _, r in d.iterrows()]
    add("table2_main_n20_a010", "conventional_windows",
        "Average coverage of conventional and calibrant-only conformal windows (n = 20, α = 0.10; 228 systems, leave-source-out; per-system "
        "means over 20 draws).", SUMMARY_HEAD, rows, d)


def t_xsc_source():
    d = rd("xsc_table_per_source")
    rows = [[r["Held-out source"], int(r["Systems"]), r["Method"], num(r["Mean retention"]), num(r["Labs < 0.90 (%)"], 1), num(r["Median width (min)"], 2)] for _, r in d.iterrows()]
    add("xsc_table_per_source", "xsc_per_source",
        "Laboratory-level failure rate (retention < 0.90) and median width per held-out source group (n = 20 random calibrants).",
        ["Held-out source", "Systems", "Method", "Mean retention", "Failing laboratories (%)", "Median width (min)"], rows, d)


def t_xsc_norman():
    if not os.path.exists(os.path.join(TAB, "xsc_table_norman.csv")):
        return
    d = rd("xsc_table_norman")
    rows = [[r["Scenario"], r["n"], r["Method"], int(r["Methods (labs)"]), num(r["Mean retention"]), num(r["Labs < 0.90 (%)"], 1),
             num(r["Labs < 0.85 (%)"], 1), num(r["Median width (min)"], 2)] for _, r in d.iterrows()]
    add("xsc_table_norman", "xsc_norman",
        "External validation of XSC on the 37 NORMAN interlaboratory methods (multipliers calibrated on RepoRT only). Random scenarios: 20 "
        "calibrant draws per method; laboratory mix: each laboratory's calibration chemicals as calibrants and its suspects as analytes.",
        ["Scenario", "n", "Method", "Laboratories", "Mean retention", "Failing < 0.90 (%)", "Failing < 0.85 (%)", "Median width (min)"], rows, d)


# ------------------------------------------------------------------ document
facts, _ = build_facts(RES, TAB)
for f in [t_training_sizes, t_all_n_alpha, t_validity, t_paired_all, t_small_n, t_fixed, t_per_fold, t_fold_het, t_strata, t_superclass,
          t_worst, t_errors, t_decoys, t_decoy_tests, t_matched, t_fixed_decoys, t_ablation, t_sensitivity, t_robustness, t_mechanism,
          t_generalisation, t_norman, t_norman_strata, t_sasc, t_sasc_c, t_compute, t_xsc_rules, t_xsc_mult, t_conventional, t_xsc_source, t_xsc_norman]:
    f()
TB = {t[1]: t for t in TABLES}


def table_block(source):
    k, _, slug, cap, md, _ = TB[source]
    return f"\n{md}\n\nTable: Table S{k}. {cap}\n"


def fig_block(name, cap):
    return f"\n![{cap}]({os.path.join(FIG, name + '.png')}){{width=16cm}}\n"


proof = open(os.path.join(SUP, "S2_coverage_proof.md")).read()
xsc_proof = open(os.path.join(SUP, "S7_laboratory_level_guarantee.md")).read()
doc = f"""**Supplementary Material**

**Laboratory-level guarantees for retention-time filtering in non-target screening by cross-system conformal calibration**

Contents: S1 Dataset details · S2 Coverage proof for the normalised jackknife+ window · S3 Results for all calibrant numbers and α levels ·
S4 Per-fold results · S5 Conditional coverage and error analysis · S6 Candidate filtering, ablation, sensitivity, robustness, generalisation,
external validation, selection-aware calibration and computational cost · S7 Laboratory-level guarantee of XSC (proof) · S8 XSC results for all rules, calibrant numbers, multipliers and external laboratories. All tables are generated by `code/scripts/build_supplementary.py`
from the analysis outputs in `tables/` and `results/`; `supplementary/table_index.csv` lists the source file of every table. Sections S2–S6 analyse the average coverage of calibrant-only windows, including the normalised jackknife+ window with the meta-learned difficulty model (a comparator in the main article) and its selection-aware calibration (SASC; the special case of cross-system calibration with a jackknife+-shaped window and an average target). The proposed method, cross-system conformal calibration, is covered in Sections S7–S8.

## S1. Dataset details

Retention times were taken from RepoRT v1.0 (Zenodo 10.5281/zenodo.16267768). For each chromatographic system, entries without a valid
retention time, standardised SMILES or standard InChIKey were removed. Compounds were identified by the first 14 characters (connectivity block)
of the standard InChIKey; entries of the same connectivity block (stereoisomers and duplicate entries) were averaged if their retention times
differed by at most 0.2 min and removed otherwise. The mobile-phase pH class was derived from the reported eluent additives (acidic if any acid
additive was present in eluent A or B, basic if any basic additive and no acid, buffered if only ammonium salts, none otherwise) and the organic
modifier from the reported solvent fractions (acetonitrile, methanol, mixed or other). Column length, particle size, temperature, flow rate and
dead time were taken from the metadata and the gradient duration from the gradient table. Missing numeric metadata were encoded as 0 (temperature, flow rate and dead time were missing for {{{{mdmiss.temp|int}}}}, {{{{mdmiss.flow|int}}}} and {{{{mdmiss.t0|int}}}} of the {{{{mdmiss.n_systems|int}}}} RP systems). For the NORMAN methods, which do not report these three descriptors, they were passed as missing values (NaN); because the gradient-boosted difficulty model never saw NaN for these features during training, LightGBM routes NaN like 0 (verified with lightgbm 4.7.0: identical predictions for NaN and 0), so that both data sets are encoded identically. Reversed-phase systems
with at least 60 compounds after cleaning were target systems, except METLIN SMRT, which served only as training data; HILIC systems with at
least 60 compounds formed the out-of-mechanism test set. Fig. S1 summarises the target systems and Table S1 the numbers of training systems
available in each leave-source-out fold.
{table_block("tableS17_training_sizes")}
{fig_block("figS1_dataset_overview", "Fig. S1. Reversed-phase target systems: (a) number of systems per source group, (b) number of compounds per system and (c) retention time of the last-eluting compound.")}

{proof}

## S3. Results for all calibrant numbers and α levels
{table_block("tableS1_all_n_alpha")}
{table_block("tableS2_coverage_validity")}
{table_block("tableS3_paired_winkler_all_n_alpha")}
{table_block("tableS12_small_n")}
{table_block("table9_fixed_windows")}
{fig_block("figS2_alpha_calibration", "Fig. S2. Mean empirical coverage against the nominal coverage 1 − α (n = 20, α = 0.05, 0.10 and 0.20); dashed line: identity.")}

## S4. Per-fold results
{table_block("tableS4_per_fold_n20_a010")}
{table_block("table11_fold_heterogeneity")}

## S5. Conditional coverage and error analysis
{table_block("tableS5_strata_n20_a010")}
{table_block("tableS6_superclass_coverage")}
{table_block("tableS7_worst_stratum")}
{table_block("tableS9_error_classes")}

## S6. Candidate filtering, ablation, sensitivity, robustness, generalisation, external validation, selection-aware calibration and computational cost
{table_block("table7_candidate_filtering")}
{table_block("tableS8_decoy_tests")}
{table_block("tableS14_matched_retention")}
{table_block("tableS11_fixed_window_decoys")}
{table_block("table4_ablation_n20_a010")}
{table_block("tableS10_sensitivity")}
{table_block("table5_robustness")}
{table_block("tableS16_robustness_mechanism")}
{table_block("table6_generalisation")}
{table_block("table10_norman_external")}
{table_block("tableS13_norman_strata")}
{table_block("table12_sasc")}
{table_block("tableS15_sasc_c_values")}
{table_block("compute_benchmark.json")}

{xsc_proof}

## S8. Cross-system conformal calibration: all rules, calibrant numbers, multipliers and external laboratories
{table_block("xsc_table_rules_n")}
{table_block("xsc_table_multipliers")}
{table_block("table2_main_n20_a010")}
{table_block("xsc_table_per_source")}
{table_block("xsc_table_norman") if "xsc_table_norman" in TB else ""}
{fig_block("figX3_lab_calibration", "Fig. S3. Laboratory-level calibration of XSC: observed share of laboratories whose retention falls below 1 − ε against the nominal δ (random calibrants, all calibrant numbers and all five held-out source groups pooled); dashed line: identity.")}
"""
order_in_doc = [int(m) for m in re.findall(r"Table: Table S(\d+)\.", doc)]
assert order_in_doc == sorted(order_in_doc) == list(range(1, len(TABLES) + 1)), order_in_doc
doc, missing = resolve(doc, facts)
if re.search(r"(?<!\{)\{[a-zA-Z_.]+\|[a-z0-9.]+\}(?!\})", doc):
    sys.exit("single-brace token left in the document (f-string escaping error)")
if missing:
    sys.exit("UNRESOLVED TOKENS: " + ", ".join(sorted(set(missing))))
md_path = os.path.join(SUP, "Supplementary_Material.md"); open(md_path, "w").write(doc)
idx = []
for k, src, slug, cap, _, df in TABLES:
    fn = f"Table_S{k}_{slug}.csv"; df.to_csv(os.path.join(OUT_TAB, fn), index=False)
    idx.append({"table": f"S{k}", "file": f"supplementary/tables/{fn}", "source": src, "caption": resolve(cap, facts)[0]})
pd.DataFrame(idx).to_csv(os.path.join(SUP, "table_index.csv"), index=False)
docx = os.path.join(SUP, "Supplementary_Material.docx")
pandoc_docx(md_path, docx, os.path.join(ROOT, "manuscript", "build", "reference_aca.docx"), ROOT)
postprocess_docx(docx, line_numbers=False, table_pt=8)
print(f"supplementary: {len(TABLES)} tables; docx written: {docx}")
if not a.no_pdf:
    ok, err = docx_to_pdf(docx, SUP)
    print("pdf:", "ok" if ok else err)

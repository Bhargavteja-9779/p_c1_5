"""Derive every number quoted in the manuscript from machine-generated result files.
build_facts(results_dir, tables_dir) -> (facts: dict[str, value], tables: dict[str, markdown_table])"""
import os, json, re
import numpy as np, pandas as pd

LAB = {"B1_gauss_rmse": "Gaussian ±z·RMSE", "B2_pspline_PI": "P-spline PI", "B3_gp_PI": "GP projection PI", "B7_gp_meta": "GP projection PI, meta-learned prior",
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
    _g = lambda v: f"{v:.1f}".rstrip("0").rstrip(".")  # medians of an even count may end in .5: show it, do not round
    F["median_compounds"] = _g(d['median_compounds_per_system'])
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
    # supplementary tables flattened ---------------------------------------------------
    s1 = pd.read_csv(os.path.join(TAB, "tableS1_all_n_alpha.csv"))
    for _, r in s1.iterrows():
        for c in ["coverage_mean", "width_median_min", "winkler_mean", "share_infinite", "share_severe_under", "share_below_nominal"]:
            F[f"s1.{r.method}.n{int(r.n)}.a{r.alpha:.2f}.{c}"] = r[c]
    st_ = pd.read_csv(os.path.join(TAB, "tableS5_strata_n20_a010.csv"))
    tot = st_[st_.stratum == "seen"].groupby("group").total.max()
    F["share_seen_test"] = f"{100 * tot.get('seen', 0) / tot.sum():.0f} %"
    het = pd.DataFrame(st.get("fold_heterogeneity", []))
    for _, r in het.iterrows():
        for c in ["median_rel_change", "share_P_better", "p_wilcoxon", "median_diff_min", "rel_q1", "rel_q3"]:
            if c in r and pd.notna(r[c]):
                F[f"het.{r.comparator}.{r.fold_name}.{c}"] = r[c]
    for k, v in st.get("gain_vs_rt_span", {}).items():
        F[f"span.{k}"] = v
    for r in m.get("matched_retention", []):
        F[f"mr.{r['method']}.{r['retention']:.2f}"] = r["decoy_elim_interp"]
    for r in m.get("compute_per_draw", []):
        F[f"cpu.n{int(r['n'])}.map_ms"] = r["t_map_loo_ms"]; F[f"cpu.n{int(r['n'])}.sigma_ms"] = r["t_sigma_ms"]
    for r in st.get("decoy_tests", []):
        F[f"dectest.{r['comparator']}.{r['alpha']:.2f}.p"] = r["p_wilcoxon"]; F[f"dectest.{r['comparator']}.{r['alpha']:.2f}.diff"] = r["median_diff_elim"]
    ea = m.get("error_analysis", {})
    if ea:
        F["err.miss_rate"] = ea["P_misses_draw0"] / ea["P_evaluated_draw0"]; F["err.late_share"] = ea["side"].get("eluted_later_than_window", 0) / ea["P_misses_draw0"]
    s9 = pd.read_csv(os.path.join(TAB, "tableS9_error_classes.csv"))
    for _, r in s9.iterrows():
        F[f"errcls.{r.cls}"] = r.miss_rate
    t7 = pd.read_csv(os.path.join(TAB, "table7_candidate_filtering.csv"))
    F["dec.parents_per_system"] = t7.parents_per_system.iloc[0]; F["dec.decoys_per_system"] = t7.decoys_per_system.iloc[0]
    F["gp_miss_factor"] = (1 - F["main.B3_gp_PI.coverage_mean"]) / 0.10
    _hil = [r for r in m.get("generalisation", []) if r["setting"].startswith("HILIC") and int(r["n"]) == 20 and r["method"] == "P"]
    if _hil:
        F["n_hilic_systems"] = int(_hil[0]["n_systems"])
    SC = {"random calibrants (main)": "main", "RT-spanning calibrants": "span", "commonly measured calibrants": "common",
          "calibrant RT noise (2% of span)": "noise", "one outlying calibrant (+10% span)": "outlier"}
    for r in m.get("robustness", []):
        k = SC.get(r["scenario"], r["scenario"])
        for c in ["coverage_mean", "width_median_min", "width_mean_min", "winkler_mean", "share_severe_under"]:
            F[f"rob.{k}.{r['method']}.{c}"] = r[c]
    for r in m.get("robustness_mechanism", []):
        k = {"main": "main", "robust_span_none": "span", "robust_common_none": "common", "robust_random_noise": "noise", "robust_random_outlier": "outlier"}.get(r["scenario"], r["scenario"])
        for c in ["share_outside_cal_range", "cov_outside", "cov_inside", "cov_late"]:
            F[f"mech.{k}.{r['method']}.{c}"] = r[c]
    for k, v in m.get("metadata_missing_rp", {}).items():
        F[f"mdmiss.{k}"] = v if not isinstance(v, list) else ", ".join(v)
    if m.get("sigma_rows_main"):
        F["cb.sigma_rows_fold_min"] = m["sigma_rows_main"]["min"]; F["cb.sigma_rows_fold_max"] = m["sigma_rows_main"]["max"]
    for r in st.get("coverage_validity", []):  # pre-specified one-sided test of under-coverage
        F[f"valp.{r['method']}.n{int(r['n'])}.a{r['alpha']:.2f}"] = r.get("p_under_one_sided")
        F[f"valpc.{r['method']}.n{int(r['n'])}.a{r['alpha']:.2f}"] = r.get("cluster_p_under_one_sided")
    for scen_, rows_ in st.get("norman_tests", {}).items():
        for r in rows_:
            for k in ["median_diff", "p_wilcoxon", "p_holm", "share_ref_better"]:
                F[f"nortest.{scen_}.{r['metric']}.{r['comparator']}.{k}"] = r[k]
    for set_, rows_ in st.get("generalisation_tests", {}).items():
        for r in rows_:
            for k in ["median_diff", "p_wilcoxon", "p_holm", "share_ref_better"]:
                F[f"gentest.{set_}.{r['metric']}.{r['comparator']}.{k}"] = r[k]
    _pf = pd.DataFrame(st.get("per_fold", []))
    if len(_pf):
        for m_, g_ in _pf.groupby("method"):
            F[f"foldcov.{m_}.min"] = g_.coverage.min(); F[f"foldcov.{m_}.max"] = g_.coverage.max()
    tsz = pd.DataFrame(m.get("training_sizes", []))
    for mt_ in ["RP", "HILIC"]:
        t_ = tsz[tsz["mode"] == mt_] if len(tsz) else tsz
        if len(t_):
            F[f"tsz.{mt_}.ge30.min"] = int(t_.training_systems_ge30.min()); F[f"tsz.{mt_}.ge30.max"] = int(t_.training_systems_ge30.max())
    for r in m.get("fixed_windows", []):
        for c in ["coverage_mean", "width_median_min", "winkler_mean", "share_severe_under", "coverage_p10"]:
            F[f"fw.{r['method']}.n{int(r['n'])}.{c}"] = r[c]
    for r in m.get("fixed_window_decoys", []):
        F[f"fwdec.{r['method']}.elim"] = r["decoy_elim_mean"]; F[f"fwdec.{r['method']}.ret"] = r["true_retained_mean"]
    for r in m.get("small_n", []):
        for c in ["coverage_mean", "width_median_min", "share_infinite"]:
            F[f"sn.{r['method']}.n{int(r['n'])}.a{r['alpha']:.2f}.{c}"] = r[c]
    for r in m.get("sensitivity", []):
        F[f"sens.{r['setting']}.{r['method']}.winkler"] = r["winkler_mean"]; F[f"sens.{r['setting']}.{r['method']}.coverage"] = r["coverage_mean"]
    for r in m.get("sasc", []):
        for c in ["coverage_mean", "width_median_min", "winkler_mean", "share_severe_under", "coverage_ci_lo", "coverage_ci_hi"]:
            F[f"sasc.{r['rule']}.{r['method']}.{c}"] = r[c]
    cv_ = pd.DataFrame(m.get("sasc_c", []))
    for rule_, g_ in (cv_.groupby("rule") if len(cv_) else []):
        F[f"sascc.{rule_}.P.min"] = g_.c_P.min(); F[f"sascc.{rule_}.P.max"] = g_.c_P.max()
        F[f"sascc.{rule_}.J.min"] = g_.c_J.min(); F[f"sascc.{rule_}.J.max"] = g_.c_J.max()
        F[f"sascc.{rule_}.ntrain.min"] = int(g_.n_train_systems.min()); F[f"sascc.{rule_}.ntrain.max"] = int(g_.n_train_systems.max())
    for k_, v_ in (m.get("sasc_norman_c") or {}).items():
        F[f"sascN.{k_}"] = v_
    for r in m.get("norman", []):
        nk = "mix" if r["n"] == "lab mix" else str(r["n"])
        for c in ["coverage_mean", "coverage_ci_lo", "coverage_ci_hi", "width_median_min", "winkler_mean", "share_severe_under", "n_systems"]:
            F[f"nor.{r['scenario']}.{nk}.a{float(r['alpha']):.2f}.{r['method']}.{c}"] = r[c]
    if m.get("norman_meta"):
        F["norman.n_systems"] = m["norman_meta"]["n_systems"]; F["norman.n_compounds"] = m["norman_meta"]["n_compounds"]; F["norman.seen_share"] = m["norman_meta"]["seen_share"]
    GS = {"SMRT-regression base model (RP)": "smrt", "HILIC systems (order model)": "hilic"}
    for r in m.get("generalisation", []):
        k = GS.get(r["setting"], r["setting"])
        for c in ["coverage_mean", "width_median_min", "winkler_mean", "share_severe_under", "n_systems"]:
            F[f"gen.{k}.n{int(r['n'])}.{r['method']}.{c}"] = r[c]
    cb = m.get("compute_benchmark")
    if cb:
        F["cb.base_train_min"] = cb["base_train_s_1thread"] / 60; F["cb.sigma_train_s"] = cb["sigma_train_s_1thread"]; F["cb.sigma_rows"] = cb["sigma_train_rows"]
        F["cb.peak_rss_MB"] = cb["peak_rss_MB"]; F["cb.base_rows"] = cb["base_train_rows"]
        for n, v in cb["per_system"].items():
            F[f"cb.n{n}.cal_ms"] = v["calibration_ms"]; F[f"cb.n{n}.us_per_cand"] = v["windows_us_per_candidate"]
    nr_ = pd.read_parquet(os.path.join(os.path.dirname(RES), "results", "raw", "norman_rows.parquet"), columns=["scenario", "n", "cs"])
    nn_ = nr_[nr_.scenario == "fixed_mix"].groupby("cs").n.first()
    F["nor.mix_n.min"] = int(nn_.min()); F["nor.mix_n.max"] = int(nn_.max())
    # share of test-compound evaluations in each group of a stratum (e.g. thirds of the RT range)
    for st_ in {k.split(".")[1] for k in F if k.startswith("strata_n.")}:
        grp = {k.split(".", 2)[2]: v for k, v in F.items() if k.startswith(f"strata_n.{st_}.")}
        tot_ = sum(grp.values())
        for g_, v_ in grp.items():
            F[f"strata_share.{st_}.{g_}"] = v_ / tot_ if tot_ else float("nan")
    if "main.B7_gp_meta.width_median_min" in F:
        F["ratio.b7_p_width"] = F["main.B7_gp_meta.width_median_min"] / F["main.P.width_median_min"]
    # coverage range over the conformal constructions (split conformal, jackknife+, ensemble-normalised jackknife+, proposed)
    CONF = ["B4_split_conformal", "B5_jackknife_plus", "B6_jk+_ensemble", "P"]
    for tag, pat in [("main", "main.{m}.coverage_mean"), ("hilic", "gen.hilic.n20.{m}.coverage_mean"), ("norrand", "nor.random_all.20.a0.10.{m}.coverage_mean")]:
        vals = [F[pat.format(m=m_)] for m_ in CONF if pat.format(m=m_) in F]
        if vals:
            F[f"confrange.{tag}.min"] = min(vals); F[f"confrange.{tag}.max"] = max(vals)
    # tables (markdown) ------------------------------------------------------------
    from docbuild import mdtab  # typographic minus, content-proportional column widths
    # Table 1 — data
    t1 = pd.read_csv(os.path.join(TAB, "table1_dataset.csv"))
    T["table1"] = mdtab(["Source group (fold)", "Target systems", "Median compounds per system", "Retention times", "Median RT span (min)", "Acidic / non-acidic mobile phase", "ACN / MeOH"],
                        [[r.fold_name, int(r.systems), f"{r.compounds_median:.1f}".rstrip("0").rstrip("."), f"{int(r.compounds_total):,}", f"{r.rt_span_median_min:.1f}", f"{int(r.acidic)} / {int(r.non_acidic)}", f"{int(r.acn)} / {int(r.meoh)}"] for _, r in t1.iterrows()])
    # Table 2 — main comparison incl. fixed tolerances (n = 20, alpha = 0.10)
    t2 = pd.read_csv(os.path.join(TAB, "table2_main_n20_a010.csv"))
    fw = pd.read_csv(os.path.join(TAB, "table9_fixed_windows.csv")); fw = fw[(fw.n == 20) & fw.method.isin(["RMSE_x2", "FIX_0.5min", "FIX_1min", "FIX_2min", "FIX_10pct"])]
    rows = []
    for _, r in pd.concat([t2, fw], ignore_index=True).iterrows():
        rows.append([LAB.get(r.method, r.method), f"{_f(r.coverage_mean)} [{_f(r.coverage_ci_lo)}, {_f(r.coverage_ci_hi)}]", f"{100*r.share_severe_under:.0f}",
                     _f(r.width_median_min, 2), _f(r.winkler_mean, 2)])
    T["table2"] = mdtab(["Window construction", "Mean coverage [95 % CI]", "Systems with coverage < 0.85 (%)", "Median width (min)", "Mean interval score (min)"], rows)
    # Table 3 — paired comparisons of the proposed windows (interval score), overall and per source fold
    t3 = pd.read_csv(os.path.join(TAB, "table3_paired_tests_n20_a010.csv")); t3 = t3[t3.metric == "winkler"]
    het = pd.read_csv(os.path.join(TAB, "table11_fold_heterogeneity.csv"))
    rows = []
    for _, r in t3.iterrows():
        hh = het[(het.comparator == r.comparator) & (het.fold_name != "all")]
        order = ["BGC", "Folberth", "PredRet", "Souihi+Stoffel", "Other sources"]
        hh = hh.set_index("fold_name").reindex([o for o in order if o in set(hh.fold_name)]) if len(hh) else hh
        per = " / ".join(f"{100*h.median_rel_change:+.0f}" for _, h in hh.iterrows()) if len(hh) else "–"
        rows.append([LAB.get(r.comparator, r.comparator), f"{r.median_diff:+.2f} [{r.diff_ci_lo:+.2f}, {r.diff_ci_hi:+.2f}]", f"{100*r.share_ref_better:.0f}",
                     f"{r.rank_biserial:+.2f}", _p(r.p_holm), per])
    T["table3"] = mdtab(["Comparator", "Median paired difference in interval score (min) [95 % CI]", "Systems where proposed is better (%)", "Rank-biserial r", "Holm-adjusted p", "Median relative change per fold (%)"], rows)
    # Table 4 — robustness, generalisation, external validation and SASC (n = 20, alpha = 0.10)
    def cov(key):
        v = F.get(key); return "–" if v is None else _f(v)
    scen = [("Random calibrants (main)", "main.{m}.coverage_mean"), ("RT-spanning calibrants", "rob.span.{m}.coverage_mean"),
            ("Commonly measured calibrants", "rob.common.{m}.coverage_mean"), ("Calibrant RT noise (SD 2 % of span)", "rob.noise.{m}.coverage_mean"),
            ("One outlying calibrant (+10 % of span)", "rob.outlier.{m}.coverage_mean"), ("SMRT-regression base model", "gen.smrt.n20.{m}.coverage_mean"),
            ("HILIC systems", "gen.hilic.n20.{m}.coverage_mean"), ("NORMAN, random calibrants", "nor.random_all.20.a0.10.{m}.coverage_mean"),
            ("NORMAN, laboratory calibrant mix", "nor.fixed_mix.mix.a0.10.{m}.coverage_mean")]
    rows = [[name] + [cov(k.format(m=m)) for m in ["P", "B5_jackknife_plus", "B1_gauss_rmse", "B3_gp_PI"]] for name, k in scen]
    rows.append(["Commonly measured calibrants + SASC", cov("sasc.common.P+SASC.coverage_mean"), cov("sasc.common.J+SASC.coverage_mean"), "–", "–"])
    rows.append(["NORMAN laboratory mix + SASC", cov("nor.fixed_mix.mix.a0.10.P+SASC.coverage_mean"), cov("nor.fixed_mix.mix.a0.10.J+SASC.coverage_mean"), "–", "–"])
    T["table4"] = mdtab(["Scenario (n = 20, nominal 0.90)", "Proposed", "Jackknife+", "Gaussian ±z·RMSE", "GP projection PI"], rows)
    xsc_facts(F, T, RES, TAB)
    return F, T


XLAB = {"B1_gauss_rmse": "Gaussian ±z·RMSE", "B5_jackknife_plus": "Jackknife+", "P": "Jackknife+, meta-learned σ",
        "JP_lab": "Jackknife+ × XSC multiplier (B8)", "XSC_flat_avg": "XSC, σ ≡ 1, average target", "XSC_flat_lab": "XSC, σ ≡ 1 (B9)",
        "XSC_sig_avg": "XSC, σ, no shrinkage, average", "XSC_sig_lab": "XSC, σ, no shrinkage", "XSC_sigs_avg": "XSC, average target",
        "XSC_sigs_lab": "**XSC, laboratory-level (proposed)**"}


def xsc_facts(F, T, RES, TAB):
    """Facts and tables of Addendum G (cross-system conformal calibration)."""
    path = os.path.join(RES, "xsc_metrics.json")
    if not os.path.exists(path):
        return
    M = json.load(open(path))
    S = list(M["S_training_systems"].values()); F["x.S_min"], F["x.S_max"] = min(S), max(S)
    F["x.n_systems"] = M["n_target_systems"]
    for k, v in M["summary"].items():
        rule, n, meth, lev = k.split("|")
        pre = "x" if abs(float(lev) - 0.1) < 1e-9 else f"xa{float(lev):.2f}"
        for s, val in v.items():
            F[f"{pre}.{rule}.{n}.{meth}.{s}"] = val
    for meth, (lo, hi) in M.get("main_ci", {}).items():
        F[f"x.ci.{meth}.lo"], F[f"x.ci.{meth}.hi"] = lo, hi
    for k, v in M.items():
        if k.startswith("ratio_"):
            for s, val in v.items():
                F[f"x.{k}.{s}"] = val
    for k, v in M.get("per_source", {}).items():
        for s, val in v.items():
            F[f"x.src.{k}.{s}"] = val
    srcfail = [v["fail"] for k, v in M.get("per_source", {}).items() if k.endswith(".XSC_sigs_lab")]
    if srcfail:
        F["x.src.max_fail_proposed"], F["x.src.min_fail_proposed"] = max(srcfail), min(srcfail)
    jf = [v["fail"] for k, v in M.get("per_source", {}).items() if k.endswith(".B5_jackknife_plus")]
    if jf:
        F["x.src.max_fail_jp"], F["x.src.min_fail_jp"] = max(jf), min(jf)
    off = os.path.join(RES, "xsc_offline_metrics.json")
    if os.path.exists(off):
        O = json.load(open(off))
        for k, (fa, fe) in O.get("lab_calibration", {}).items():
            rule, eps, dl = k.split("|"); F[f"xo.cal.{rule}.{float(eps):.2f}.{float(dl):.2f}"] = fa; F[f"xo.calx.{rule}.{float(eps):.2f}.{float(dl):.2f}"] = fe
        for k, v in O.get("local", {}).items():
            rule, n, meth = k.split("|")
            for s, val in v.items():
                F[f"xo.loc.{rule}.{n}.{meth}.{s}"] = val
    F["x.mult_diff_q90"] = M.get("multiplier_exact_vs_avg_maxreldiff", float("nan"))
    for meth, r in M.get("candidates", {}).items():
        for s, val in r.items():
            F[f"xc.{meth}.{s}"] = val
    for meth, v in M.get("candidates_lab_fail", {}).items():
        F[f"xc.labfail.{meth}"] = v
    for k in ("n_pcl", "median_list"):
        if k in M:
            F[f"xc.{k}"] = M[k]
    for k, r in M.get("norman", {}).items():
        sc, n, meth = k.split("|")
        for s, val in r.items():
            F[f"xn.{sc}.{n.replace(' ', '_')}.{meth}.{s}"] = val
    for meth, r in M.get("case", {}).items():
        for s, val in r.items():
            F[f"xk.{meth}.{s}"] = val
    if "case_median_list" in M:
        F["xk.median_list"] = M["case_median_list"]
    if "norman_S_mix" in M:
        F["xn.S_mix"] = M["norman_S_mix"]
    from docbuild import mdtab

    def tab(name, header=None):
        f = os.path.join(TAB, name + ".csv")
        if not os.path.exists(f):
            return
        d = pd.read_csv(f)
        cols = []
        for c in d.columns:
            v = d[c]
            if v.dtype.kind == "f":
                mx = v.replace([np.inf, -np.inf], np.nan).abs().max()
                dd = 3 if mx < 1.5 else (2 if mx < 20 else 1)
                cols.append(v.map(lambda x, dd=dd: _f(x, dd)))
            else:
                cols.append(v.astype(str))
        rows = [list(r) for r in zip(*cols)]
        T[name] = mdtab(header or list(d.columns), rows)
    tab("xsc_table_main_n20"); tab("xsc_table_per_source"); tab("xsc_table_candidates"); tab("xsc_table_norman"); tab("xsc_table_case")

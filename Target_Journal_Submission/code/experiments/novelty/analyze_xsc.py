"""Analysis of Addendum G (XSC): tables, statistics (system-cluster bootstrap) and figures.
Usage: python experiments/novelty/analyze_xsc.py   (from code/; reads ../results/raw/xsc, writes ../results/xsc_metrics.json,
../tables/xsc_*.{csv,md}, ../figures/figX*_*.{png,pdf})"""
import os, sys, json, glob
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
os.chdir(os.path.join(HERE, "..", ".."))
from conformrt.visualization.figures import plt, INK, INK2, MUTED, GRID, W1, W15, W2, _save as _fsave  # noqa
from conformrt.evaluation.protocol import FOLD_NAMES

def _save(fig, path):
    _fsave(fig, os.path.dirname(path), os.path.basename(path)); plt.close(fig)


RAW = "../results/raw/xsc"; TAB = "../tables"; FIG = "../figures"
B = int(os.environ.get("XSC_BOOT", 2000)); RNG = np.random.default_rng(0)
LABEL = {"B1_gauss_rmse": "±z·RMSE (Gaussian)", "B5_jackknife_plus": "Jackknife+", "P": "Jackknife+, meta-learned σ",
         "JP_lab": "Jackknife+, lab-calibrated", "XSC_flat_avg": "XSC flat (average)", "XSC_flat_lab": "XSC flat (lab-level)",
         "XSC_sig_avg": "XSC σ (average)", "XSC_sig_lab": "XSC σ (lab-level)", "XSC_sigs_avg": "XSC (average)",
         "XSC_sigs_lab": "XSC (lab-level, proposed)"}
COL = {"XSC_sigs_lab": "#2a78d6", "B1_gauss_rmse": "#eb6834", "B5_jackknife_plus": "#1baf7a", "JP_lab": "#eda100", "P": "#e87ba4",
       "XSC_sigs_avg": "#2a78d6", "XSC_flat_lab": "#52514e"}
MRK = {"XSC_sigs_lab": "o", "B1_gauss_rmse": "s", "B5_jackknife_plus": "D", "JP_lab": "^", "P": "v", "XSC_sigs_avg": "o", "XSC_flat_lab": "P"}
MAIN = ["B1_gauss_rmse", "B5_jackknife_plus", "P", "JP_lab", "XSC_flat_avg", "XSC_flat_lab", "XSC_sigs_avg", "XSC_sigs_lab"]


def boot_ci(df, stat, by="system"):
    """System-cluster bootstrap CI of stat(df) (stat: function of a frame)."""
    groups = df.groupby(by).indices; keys = list(groups)
    vals = []
    for _ in range(B):
        pick = RNG.choice(len(keys), len(keys), replace=True)
        idx = np.concatenate([groups[keys[k]] for k in pick])
        vals.append(stat(df.iloc[idx]))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def summarise(d, eps):
    fin = np.isfinite(d.width_median)
    return {"n_draws": int(len(d)), "n_systems": int(d.system.nunique()), "coverage_mean": float(d.coverage.mean()),
            "fail": float((d.coverage < 1 - eps).mean()), "fail_085": float((d.coverage < 0.85).mean()),
            "q05": float(d.coverage.quantile(0.05)), "share_infinite": float((~fin).mean()),
            "width_median": float(d.width_median.median()), "width_rel_median": float((d.width_median / d.span).median()),
            "is_median": float(d.winkler.median())}


def md_table(df, path, floatfmt=3):
    df.to_csv(path + ".csv", index=False)
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        lines.append("| " + " | ".join((f"{v:.{floatfmt}f}" if isinstance(v, float) and np.isfinite(v) else ("∞" if isinstance(v, float) else str(v))) for v in r) + " |")
    open(path + ".md", "w").write("\n".join(lines) + "\n")


def main():
    M = {}
    ev = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f"{RAW}/fold*_eval.parquet"))], ignore_index=True)
    M["folds_done"] = sorted(int(x) for x in ev.fold.unique())
    M["n_target_systems"] = int(ev.system.nunique())
    mult = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{RAW}/fold*_multipliers.csv"))], ignore_index=True)
    M["S_training_systems"] = {str(k): int(v) for k, v in mult[(mult.rule == "random") & (mult.n == 20)].groupby("fold").S.first().items()}
    # ---------------- per (rule, n, method, level) summaries ----------------
    rows = []
    for (rule, n, meth, lev), d in ev.groupby(["rule", "n", "method", "level"]):
        rows.append(dict(rule=rule, n=int(n), method=meth, level=float(lev), **summarise(d, lev)))
    S = pd.DataFrame(rows); S.to_csv(f"{TAB}/xsc_all_summaries.csv", index=False)
    M["summary"] = {f"{r.rule}|{r.n}|{r.method}|{r.level}": {k: getattr(r, k) for k in ["coverage_mean", "fail", "fail_085", "q05", "share_infinite", "width_median", "width_rel_median", "is_median", "n_systems", "n_draws"]} for r in S.itertuples()}
    # ---------------- main table: random, n = 20, α = ε = 0.10, with bootstrap CIs ----------------
    main = ev[(ev.rule == "random") & (ev.n == 20) & (ev.level == 0.10)]
    T = []
    for meth in MAIN:
        d = main[main.method == meth]
        if not len(d):
            continue
        s = summarise(d, 0.10)
        lo, hi = boot_ci(d, lambda x: (x.coverage < 0.90).mean())
        T.append({"Method": LABEL[meth], "Mean retention": s["coverage_mean"], "Labs < 0.90 (%)": 100 * s["fail"],
                  "95% CI (%)": f"{100*lo:.1f}–{100*hi:.1f}", "Labs < 0.85 (%)": 100 * s["fail_085"], "5th pct retention": s["q05"],
                  "Median width (min)": s["width_median"], "Width / RT span": s["width_rel_median"], "Median IS": s["is_median"]})
        M.setdefault("main_ci", {})[meth] = [lo, hi]
    md_table(pd.DataFrame(T), f"{TAB}/xsc_table_main_n20")
    # ---------------- per source group (held-out fold) failure rates ----------------
    P = []
    for f, d in main.groupby("fold"):
        for meth in ["B1_gauss_rmse", "B5_jackknife_plus", "P", "JP_lab", "XSC_sigs_avg", "XSC_sigs_lab"]:
            x = d[d.method == meth]
            if len(x):
                M.setdefault("per_source", {})[f"f{int(f)}.{meth}"] = {"fail": float((x.coverage < 0.9).mean()), "width": float(x.width_median.median()), "systems": int(x.system.nunique())}
                P.append({"Held-out source": FOLD_NAMES[int(f)], "Systems": int(x.system.nunique()), "Method": LABEL[meth],
                          "Mean retention": x.coverage.mean(), "Labs < 0.90 (%)": 100 * (x.coverage < 0.9).mean(),
                          "Median width (min)": x.width_median.median()})
    Pf = pd.DataFrame(P); md_table(Pf, f"{TAB}/xsc_table_per_source")

    # ---------------- width of proposed vs lab-calibrated jackknife+ (paired by system) ----------------
    pv = main[main.method.isin(["XSC_sigs_lab", "JP_lab", "XSC_flat_lab"])].groupby(["system", "method"]).width_median.median().unstack()
    for a, b_ in [("XSC_sigs_lab", "JP_lab"), ("XSC_sigs_lab", "XSC_flat_lab")]:
        r = (pv[a] / pv[b_]).replace([np.inf, -np.inf], np.nan).dropna()
        M[f"ratio_{a}_vs_{b_}"] = {"median": float(r.median()), "share_narrower": float((r < 1).mean()), "n": int(len(r))}
    # ---------------- n dependence and rules ----------------
    G = S[(S.level == 0.10) & S.method.isin(["B1_gauss_rmse", "B5_jackknife_plus", "P", "JP_lab", "XSC_flat_lab", "XSC_sigs_avg", "XSC_sigs_lab"])]
    md_table(G[["rule", "n", "method", "coverage_mean", "fail", "fail_085", "share_infinite", "width_median", "is_median"]]
             .sort_values(["rule", "n", "method"]), f"{TAB}/xsc_table_rules_n")
    # ---------------- multipliers ----------------
    md_table(mult[mult.family.isin(["flat", "sigs", "jp"])][["fold", "rule", "n", "family", "S", "c_avg_0.1", "c_lab_0.1", "c_labx_0.1"]],
             f"{TAB}/xsc_table_multipliers")
    M["multiplier_exact_vs_avg_maxreldiff"] = float(((mult["c_labx_0.1"] - mult["c_lab_0.1"]).abs() / mult["c_lab_0.1"]).replace([np.inf], np.nan).dropna().quantile(0.9))
    # ---------------- candidate lists ----------------
    cf = sorted(glob.glob(f"{RAW}/fold*_candidates.parquet"))
    if cf:
        c = pd.concat([pd.read_parquet(f) for f in cf], ignore_index=True)
        C = []
        for meth, d in c.groupby("method"):
            C.append({"Method": LABEL.get(meth, meth), "key": meth, "Features": int(len(d)), "Correct retained": d.true_kept.mean(),
                      "Candidates eliminated": d.elim.mean(), "Median candidates kept": float(d.n_kept.median()),
                      "Median list size": float(d.n_cand.median()),
                      "Top-1 (retained)": float(((d["rank"] == 1) & d.true_kept).mean()),
                      "Top-3 (retained)": float(((d["rank"] <= 3) & d.true_kept).mean())})
        M["median_list"] = float(c.n_cand.median())
        z = np.load(os.path.join("cache", "pubchemlite", "cand_features.npz"), allow_pickle=True) if os.path.exists(os.path.join("cache", "pubchemlite", "cand_features.npz")) else None
        M["n_features_evaluated"] = int(c[["system", "par"]].drop_duplicates().shape[0])
        pcl = "data/raw/pubchemlite/PubChemLite_exposomics_20260925.csv"
        if os.path.exists(pcl):
            M["n_pcl"] = f"{int(pd.read_csv(pcl, usecols=['FirstBlock']).FirstBlock.nunique()):,}"
        if z is not None:
            M["n_pcl_added"] = int(len(z["ik14"]))
        Cd = pd.DataFrame(C)
        M["candidates"] = {r["key"]: {"kept": r["Correct retained"], "elim": r["Candidates eliminated"], "nkept": r["Median candidates kept"],
                                      "top1": r["Top-1 (retained)"], "top3": r["Top-3 (retained)"], "features": r["Features"]} for r in C}; md_table(Cd.drop(columns="key"), f"{TAB}/xsc_table_candidates")
        # per-laboratory: share of laboratories (system x draw) retaining < 90 % of correct structures
        lab = c.groupby(["method", "system", "draw"]).true_kept.mean().reset_index()
        M["candidates_lab_fail"] = {m: float((g.true_kept < 0.9).mean()) for m, g in lab.groupby("method")}
    # ---------------- NORMAN ----------------
    if os.path.exists(f"{RAW}/norman_eval.parquet"):
        nm = pd.read_parquet(f"{RAW}/norman_eval.parquet"); nm["system"] = nm.cs
        Nn = []
        for (sc, n, meth), d in nm.assign(n=np.where(nm.scenario == "fixed_mix", -1, nm.n)).groupby(["scenario", "n", "method"]):
            s = summarise(d, 0.10)
            Nn.append({"Scenario": sc, "n": "lab mix" if n == -1 else int(n), "Method": LABEL.get(meth, meth), "key": meth,
                       "Methods (labs)": s["n_systems"], "Mean retention": s["coverage_mean"], "Labs < 0.90 (%)": 100 * s["fail"],
                       "Labs < 0.85 (%)": 100 * s["fail_085"], "Median width (min)": s["width_median"], "Median IS": s["is_median"]})
        Nd = pd.DataFrame(Nn)
        M["norman"] = {f"{r['Scenario']}|{r['n']}|{r['key']}": {"cov": r["Mean retention"], "fail": r["Labs < 0.90 (%)"] / 100,
                       "fail085": r["Labs < 0.85 (%)"] / 100, "width": r["Median width (min)"], "labs": r["Methods (labs)"]} for r in Nn}
        md_table(Nd.drop(columns="key"), f"{TAB}/xsc_table_norman")
        nmult = pd.read_csv(f"{RAW}/norman_multipliers.csv"); M["norman_S_mix"] = int(nmult[nmult.rule == "norman"].S.max())
        md_table(nmult, f"{TAB}/xsc_table_norman_multipliers")
    if os.path.exists(f"{RAW}/norman_case.parquet"):
        cs = pd.read_parquet(f"{RAW}/norman_case.parquet")
        K = []
        for meth, d in cs.groupby("method"):
            lab = d.groupby("cs").true_kept.mean()
            K.append({"Method": LABEL.get(meth, {"RMSE_x2": "±2·RMSE"}.get(meth, meth)), "key": meth, "Laboratories": int(d.cs.nunique()),
                      "Suspects × labs": int(len(d)), "Correct retained": d.true_kept.mean(), "Labs < 0.90 (%)": 100 * float((lab < 0.9).mean()),
                      "Candidates eliminated": d.elim.mean(), "Median candidates kept": float(d.n_kept.median()),
                      "Top-1 (retained)": float(((d["rank"] == 1) & d.true_kept).mean()), "Top-3 (retained)": float(((d["rank"] <= 3) & d.true_kept).mean())})
        M["case"] = {r["key"]: {"kept": r["Correct retained"], "labfail": r["Labs < 0.90 (%)"] / 100, "elim": r["Candidates eliminated"],
                                "nkept": r["Median candidates kept"], "top1": r["Top-1 (retained)"], "top3": r["Top-3 (retained)"],
                                "labs": r["Laboratories"], "n": r["Suspects × labs"]} for r in K}
        M["case_median_list"] = float(cs[cs.method == cs.method.iloc[0]].n_cand.median())
        md_table(pd.DataFrame(K).drop(columns="key"), f"{TAB}/xsc_table_case")
    json.dump(M, open("../results/xsc_metrics.json", "w"), indent=1, default=float)
    figures(ev, S)
    print("analysis done; folds", M["folds_done"])


def figures(ev, S):
    main = ev[(ev.rule == "random") & (ev.n == 20) & (ev.level == 0.10)]
    meths = ["B1_gauss_rmse", "B5_jackknife_plus", "P", "JP_lab", "XSC_sigs_lab"]
    fig, ax = plt.subplots(1, 3, figsize=(W2, 62 / 25.4))
    # (a) ECDF of laboratory-level retention
    for m in meths:
        x = np.sort(main[main.method == m].coverage.values); yv = np.arange(1, len(x) + 1) / len(x)
        ax[0].step(x, yv, where="post", color=COL[m], lw=1.6 if m == "XSC_sigs_lab" else 1.1, label=LABEL[m])
    ax[0].axvline(0.90, color=INK2, lw=0.8, ls=(0, (3, 2))); ax[0].axhline(0.10, color=MUTED, lw=0.6, ls=(0, (1, 2)))
    ax[0].set_xlim(0.6, 1.0); ax[0].set_ylim(0, 1); ax[0].set_xlabel("Retention of correct structures in one laboratory")
    ax[0].set_ylabel("Share of laboratories (cumulative)"); ax[0].set_title("a  n = 20 random calibrants", loc="left")
    ax[0].text(0.902, 0.95, "target 0.90", fontsize=6.5, color=INK2, va="top")
    ax[0].legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(0.0, 0.88))
    # (b) failure rate vs n ; (c) width vs n
    for m in ["B1_gauss_rmse", "B5_jackknife_plus", "P", "JP_lab", "XSC_sigs_lab"]:
        d = S[(S.rule == "random") & (S.level == 0.10) & (S.method == m)].sort_values("n")
        if not len(d):
            continue
        ax[1].plot(d.n, 100 * d.fail, marker=MRK[m], ms=4, color=COL[m], lw=1.2)
        dd = d[d.share_infinite < 0.5]
        ax[2].plot(dd.n, dd.width_median, marker=MRK[m], ms=4, color=COL[m], lw=1.2)
    ax[1].axhline(10, color=INK2, lw=0.8, ls=(0, (3, 2))); ax[1].text(41, 11, "δ = 10 %", fontsize=6.5, color=INK2, ha="right")
    ax[1].set_xscale("log"); ax[2].set_xscale("log")
    for a in ax[1:]:
        a.set_xticks([3, 5, 10, 20, 40]); a.set_xticklabels(["3", "5", "10", "20", "40"]); a.set_xlabel("Number of calibrants n")
    ax[1].set_ylabel("Laboratories with retention < 0.90 (%)"); ax[1].set_title("b  laboratory-level failures", loc="left")
    ax[2].set_ylabel("Median window width (min)"); ax[2].set_title("c  price: window width", loc="left"); ax[2].set_yscale("log")
    for a in ax:
        a.grid(True, color=GRID, lw=0.5); a.set_axisbelow(True)
    fig.tight_layout(); _save(fig, os.path.join(FIG, "figX1_lab_level"))
    # rules figure
    rules = ["random", "span", "common", "polar"]; meths = ["B1_gauss_rmse", "B5_jackknife_plus", "P", "XSC_sigs_avg", "XSC_sigs_lab"]
    fig, ax = plt.subplots(1, 2, figsize=(W2, 58 / 25.4), sharey=True)
    for k, n in enumerate([5, 20]):
        w = 0.16
        for i, m in enumerate(meths):
            v = [S[(S.rule == r) & (S.n == n) & (S.level == 0.10) & (S.method == m)].fail.values for r in rules]
            v = [100 * x[0] if len(x) else np.nan for x in v]
            ax[k].bar(np.arange(4) + (i - 2) * w, v, width=w * 0.88, color=COL.get(m, MUTED), alpha=0.55 if m == "XSC_sigs_avg" else 1.0,
                      label=LABEL[m], hatch="///" if m == "XSC_sigs_avg" else None, edgecolor="white", lw=0.5)
        ax[k].axhline(10, color=INK2, lw=0.8, ls=(0, (3, 2)))
        ax[k].set_xticks(range(4)); ax[k].set_xticklabels(["random", "RT-spanning", "most common", "most polar"])
        ax[k].set_title(f"{'ab'[k]}  n = {n}" + (" (jackknife+ not finite)" if n < 9 else ""), loc="left")
        ax[k].grid(True, axis="y", color=GRID, lw=0.5); ax[k].set_axisbelow(True)
    ax[0].set_ylabel("Laboratories with retention < 0.90 (%)"); ax[1].legend(frameon=False, fontsize=6, loc="upper left")
    fig.tight_layout(); _save(fig, os.path.join(FIG, "figX2_selection_rules"))


if __name__ == "__main__":
    main()

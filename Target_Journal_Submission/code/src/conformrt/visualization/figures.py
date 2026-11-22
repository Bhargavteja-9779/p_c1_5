"""Publication figures (matplotlib; PNG 300 dpi + vector PDF). All figures are drawn from machine-generated result files.
Colour rules (validated with the dataviz palette validator, light mode): categorical slots in fixed order
#2a78d6, #eb6834, #1baf7a, #eda100, #e87ba4, #008300, #6250d6; proposed method always slot 1 (blue);
baselines in neutral ink when >3 methods share a dot/strip plot; identity never by colour alone (direct labels + markers)."""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8f8d87", "#e4e2dc"
MARKERS = ["o", "s", "D", "^", "v", "P", "X"]
# Colour and marker follow the METHOD (entity), identical in every figure (validated: adjacent pairs for line charts;
# all-pairs only for the first-three trio P / Gaussian / jackknife+ used in dot plots).
METHOD_COLORS = {"P": "#2a78d6", "B1_gauss_rmse": "#eb6834", "B5_jackknife_plus": "#1baf7a", "B4_split_conformal": "#eda100",
                 "B3_gp_PI": "#e87ba4", "B2_pspline_PI": "#008300", "B6_jk+_ensemble": "#6250d6",
                 "B7_gp_meta": "#e87ba4"}  # GP family: same hue as B3, distinguished by a dashed line and its own marker
METHOD_MARKERS = {"P": "o", "B1_gauss_rmse": "s", "B5_jackknife_plus": "D", "B4_split_conformal": "^", "B3_gp_PI": "v",
                  "B2_pspline_PI": "P", "B6_jk+_ensemble": "X", "B7_gp_meta": "*"}
METHOD_LINESTYLES = {"B7_gp_meta": (0, (3, 1.5))}


def mls(m):
    return METHOD_LINESTYLES.get(m, "-")


def mcol(m, i=0):
    return METHOD_COLORS.get(m, SLOTS[i % len(SLOTS)])


def mmark(m, i=0):
    return METHOD_MARKERS.get(m, MARKERS[i % len(MARKERS)])
MM = 1 / 25.4
W1, W15, W2 = 90 * MM, 140 * MM, 190 * MM

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.8, "axes.edgecolor": INK2,
                     "axes.linewidth": 0.6, "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK,
                     "text.color": INK, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300,
                     "figure.dpi": 150, "lines.linewidth": 1.6, "grid.color": GRID, "grid.linewidth": 0.5,
                     "pdf.fonttype": 42})


def _save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, name + ".png"), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(out, name + ".pdf"), bbox_inches="tight", metadata={"CreationDate": None})  # deterministic PDF
    plt.close(fig)


def _panel(ax, letter):
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=9, ha="right", va="bottom")


# ------------------------------------------------------------------ Fig 1: workflow schematic (no data)
def fig_workflow(out):
    fig, ax = plt.subplots(figsize=(W2, 70 * MM)); ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis("off")
    fs_ = 6.2
    boxes = [
        (0.5, 19, 22, 13, "Public RT library\n(RepoRT; other laboratories'\nsystems, leave-source-out)", "#eef4fc"),
        (25.5, 19, 22, 13, "Base model:\nretention-order score s(x)\n(pairwise ranking)", "#eef4fc"),
        (50.5, 19, 23, 13, "Target system:\nn calibrants (sᵢ, tᵢ)\nmonotone map g,\nleave-one-out refits g₋ᵢ", "#fdf0ea"),
        (76.5, 19, 23, 13, "Meta-learned difficulty\nσ(x | system, calibrants)\nfitted on simulated\ntransfers", "#eaf7f1"),
        (12.5, 1.5, 43, 12.5, "Normalised jackknife+ window C(x); lower / upper end:\nkₗ-th / kᵤ-th smallest of g₋ᵢ(s(x)) ∓ σ₋ᵢ(x)Rᵢ\nRᵢ = |tᵢ − g₋ᵢ(sᵢ)| / σ₋ᵢ(xᵢ);  P(t ∈ C) ≥ 1 − 2α", "#f3f2ef"),
        (61, 1.5, 38.5, 12.5, "Candidate filtering in NTS:\nkeep a structure if the observed RT ∈ C(x)\n→ correct structure kept with prob. ≈ 1 − α", "#f3f2ef"),
    ]
    for x, y, w, h, txt, fc in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.0", fc=fc, ec=INK2, lw=0.7))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs_, color=INK, linespacing=1.3)
    arr = dict(arrowstyle="-|>", mutation_scale=8, lw=0.9, color=INK2)
    for (x0, y0, x1, y1) in [(22.9, 25.5, 25.1, 25.5), (47.9, 25.5, 50.1, 25.5), (62, 18.7, 40, 14.3), (88, 18.7, 46, 14.3),
                             (56.0, 7.75, 60.7, 7.75)]:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), **arr))
    # inner cross-fitting path: library -> difficulty model (drawn above the boxes)
    ax.plot([11.5, 11.5, 88, 88], [32.4, 36.2, 36.2, 33.2], color=INK2, lw=0.9)
    ax.add_patch(FancyArrowPatch((88, 34.0), (88, 32.4), **arr))
    ax.text(50, 37.0, "inner cross-fitting: the calibrant protocol is simulated on training systems scored by models that never saw them",
            ha="center", va="bottom", fontsize=5.9, color=INK2)
    _save(fig, out, "fig1_workflow")


# ------------------------------------------------------------------ Fig 2: main comparison (strip + box)
def fig_main(ps, out, n=20, alpha=0.1, order=None, labels=None):
    sub = ps[(ps.n == n) & np.isclose(ps.alpha, alpha)]
    order = [m for m in order if m in sub.method.unique()]
    fig, axes = plt.subplots(1, 3, figsize=(W2, 72 * MM), sharey=True, gridspec_kw={"wspace": 0.12})
    rng = np.random.default_rng(0)
    for k, (col, xl) in enumerate([("coverage", "Empirical coverage per system"), ("width", "Mean window width per system (min)"),
                                   ("winkler", "Mean interval (Winkler) score (min)")]):
        ax = axes[k]
        for i, m in enumerate(order):
            v = sub[sub.method == m][col].values; v = v[np.isfinite(v)]
            c = SLOTS[0] if m == "P" else INK2
            y = len(order) - 1 - i
            ax.scatter(v, y + rng.uniform(-0.18, 0.18, len(v)), s=4, color=c, alpha=0.35, lw=0)
            q1, q2, q3 = np.percentile(v, [25, 50, 75])
            ax.plot([q1, q3], [y, y], color=c, lw=2.6, solid_capstyle="round")
            ax.plot([q2], [y], marker="|", ms=9, mew=1.6, color="white" if m == "P" else INK)
            ax.plot([np.mean(v)], [y], marker="D", ms=3.2, color=c, mec="white", mew=0.5)
        if col == "coverage":
            ax.axvline(1 - alpha, color=INK, lw=0.8, ls=(0, (3, 2)))
            ax.text(1 - alpha, len(order) - 0.35, f"nominal {1-alpha:.2f}", ha="center", fontsize=6.3, color=INK)
            ax.set_xlim(max(0.4, sub.coverage.min() - 0.02), 1.0)
        else:
            ax.set_xscale("log")
        ax.set_xlabel(xl); ax.grid(axis="x"); _panel(ax, "abc"[k])
    axes[0].set_yticks(range(len(order))); axes[0].set_yticklabels([labels.get(m, m) for m in order][::-1])
    for t in axes[0].get_yticklabels():
        if t.get_text() == labels.get("P"):
            t.set_color(SLOTS[0]); t.set_fontweight("bold")
    _save(fig, out, "fig2_main_comparison")


# ------------------------------------------------------------------ Fig 3: calibrant number
def fig_calibrants(summ_n, out, methods, labels, alpha=0.1):
    """summ_n: DataFrame with columns method, n, alpha, coverage_mean, winkler_mean, share_infinite.
    Cells whose windows are infinite (conformal methods with too few calibrants) are omitted, not drawn as coverage 1."""
    s = summ_n[np.isclose(summ_n.alpha, alpha)]
    fig, axes = plt.subplots(1, 2, figsize=(W2, 72 * MM), gridspec_kw={"wspace": 0.28})
    for k, (col, yl) in enumerate([("coverage_mean", "Mean coverage over systems"), ("winkler_mean", "Mean interval score (min, log scale)")]):
        ax = axes[k]
        for i, m in enumerate(methods):
            d = s[s.method == m].sort_values("n")
            d = d[np.isfinite(d["winkler_mean"]) & (d.get("share_infinite", 0) == 0)]
            if d.empty:
                continue
            ax.plot(d.n, d[col], color=mcol(m, i), marker=mmark(m, i), ls=mls(m), ms=3.8, mec="white", mew=0.5, label=labels.get(m, m))
        if col == "coverage_mean":
            ax.axhline(1 - alpha, color=INK, lw=0.8, ls=(0, (3, 2)))
            ax.text(s.n.min(), 1 - alpha - 0.004, f"nominal {1 - alpha:.2f}", ha="left", va="top", fontsize=6, color=INK)
        else:
            ax.set_yscale("log")
        ax.set_xlabel("Number of calibrants n"); ax.set_ylabel(yl); ax.grid(axis="y"); _panel(ax, "ab"[k])
        ax.set_xticks(sorted(s.n.unique()))
    axes[0].legend(frameon=False, loc="lower right", fontsize=5.8)
    inf_m = s[(s.get("share_infinite", 0) > 0)][["method", "n"]]
    if len(inf_m):
        note = "; ".join(f"{labels.get(m, m)}: n={', '.join(str(int(x)) for x in g.n)}" for m, g in inf_m.groupby("method"))
        fig.text(0.5, -0.09, "Infinite (non-informative) windows omitted — " + note, ha="center", fontsize=5.6, color=INK2)
    _save(fig, out, "fig3_calibrant_number")


# ------------------------------------------------------------------ Fig 4: conditional coverage (dot plot, <=3 methods)
STRATUM_LABELS = {"seen": "seen by base model", "unseen": "not seen by base model", "sim>=0.7": "Tanimoto ≥ 0.7",
                  "sim0.4-0.7": "Tanimoto 0.4–0.7", "sim<0.4": "Tanimoto < 0.4", "neutral": "neutral", "ionisable": "ionisable",
                  "inside_cal_range": "inside calibrant range", "outside_cal_range": "outside calibrant range",
                  "early": "first third of RT range", "middle": "middle third of RT range", "late": "last third of RT range"}


def fig_strata(sc, out, methods, labels, alpha=0.1, strata_order=None):
    s = sc.copy()
    rows = [(st, g) for st, gs in strata_order for g in gs if ((s.stratum == st) & (s.group == g)).any()]
    fig, ax = plt.subplots(figsize=(W15, (18 + 5.2 * len(rows)) * MM))
    for i, m in enumerate(methods):
        for j, (st, g) in enumerate(rows):
            d = s[(s.stratum == st) & (s.group == g) & (s.method == m)]
            if d.empty:
                continue
            y = len(rows) - 1 - j + (i - (len(methods) - 1) / 2) * 0.22
            ax.plot(d.coverage_sysavg.values[0], y, marker=mmark(m, i), ms=4.2, color=mcol(m, i), mec="white", mew=0.5,
                    ls="none", label=labels.get(m, m) if j == 0 else None)
    ax.axvline(1 - alpha, color=INK, lw=0.8, ls=(0, (3, 2)))
    prev = None
    for j, (st, g) in enumerate(rows):
        if prev is not None and st != prev:
            ax.axhline(len(rows) - 1 - j + 0.5, color=GRID, lw=0.8)
        prev = st
    yt = []
    for (st, g) in rows:
        tot = s[(s.stratum == st) & (s.group == g)].total.max()
        yt.append(f"{STRATUM_LABELS.get(g, g)}  (n={int(tot):,})")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(yt[::-1]); ax.set_xlabel("Coverage within stratum (system average)")
    ax.grid(axis="x"); ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3)
    _save(fig, out, "fig4_conditional_coverage")


# ------------------------------------------------------------------ Fig 5: candidate filtering operating curves
def fig_decoys(dec, out, methods, labels):
    fig, ax = plt.subplots(figsize=(W1 * 1.25, 72 * MM))
    for i, m in enumerate(methods):
        d = dec[dec.method == m].sort_values("alpha")
        if d.empty:
            continue
        ax.plot(d.true_retained, d.decoy_elim_rate, color=mcol(m, i), marker=mmark(m, i), ms=3.8, mec="white", mew=0.5, label=labels.get(m, m))
        for _, r in d.iterrows():
            if m == methods[0]:
                ax.annotate(f"α={r.alpha:g}", (r.true_retained, r.decoy_elim_rate), xytext=(3, -8), textcoords="offset points", fontsize=5.6, color=INK2)
    ax.set_xlabel("True-candidate retention (coverage)"); ax.set_ylabel("Isomeric decoys eliminated (fraction)")
    ax.grid(); ax.legend(frameon=False, fontsize=5.8, loc="upper right")
    _save(fig, out, "fig5_candidate_filtering")


# ------------------------------------------------------------------ Fig 6: ablation
def fig_ablation(abl, out, labels):
    """abl: rows method, d_winkler (median paired diff vs P), lo, hi, coverage_mean."""
    abl = abl.reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(W2, (20 + 6 * len(abl)) * MM), sharey=True, gridspec_kw={"wspace": 0.1})
    y = np.arange(len(abl))[::-1]
    axes[0].errorbar(abl.median_diff, y, xerr=[abl.median_diff - abl.diff_ci_lo, abl.diff_ci_hi - abl.median_diff], fmt="o",
                     color=SLOTS[0], ms=4, mec="white", mew=0.5, elinewidth=1.2, capsize=0)
    axes[0].axvline(0, color=INK, lw=0.8); axes[0].set_xlabel("Δ interval score vs full model (min)\n(median over systems, 95% CI; >0 = worse)")
    axes[1].plot(abl.coverage_mean, y, "o", color=SLOTS[0], ms=4, mec="white", mew=0.5)
    axes[1].set_xlabel("Mean coverage over systems")
    for ax, l in zip(axes, "ab"):
        ax.grid(axis="x"); _panel(ax, l)
    axes[0].set_yticks(y); axes[0].set_yticklabels([labels.get(m, m) for m in abl.method])
    _save(fig, out, "fig6_ablation")


# ------------------------------------------------------------------ Fig 7: robustness small multiples
def fig_robustness(rob, out, methods, labels, alpha=0.1):
    scen = list(dict.fromkeys(rob.scenario))
    fig, axes = plt.subplots(1, 2, figsize=(W2, (22 + 7 * len(scen)) * MM), sharey=True, gridspec_kw={"wspace": 0.1})
    for k, (col, xl) in enumerate([("coverage_mean", "Mean coverage over systems"), ("winkler_mean", "Mean interval score (min)")]):
        ax = axes[k]
        for i, m in enumerate(methods):
            for j, sc in enumerate(scen):
                d = rob[(rob.scenario == sc) & (rob.method == m)]
                if d.empty:
                    continue
                y = len(scen) - 1 - j + (i - (len(methods) - 1) / 2) * 0.13
                ax.plot(d[col].values[0], y, marker=mmark(m, i), ms=4, color=mcol(m, i), mec="white", mew=0.5, ls="none",
                        label=labels.get(m, m) if (j == 0 and k == 0) else None)
        if col == "coverage_mean":
            ax.axvline(1 - alpha, color=INK, lw=0.8, ls=(0, (3, 2)))
        ax.set_xlabel(xl); ax.grid(axis="x"); _panel(ax, "ab"[k])
    axes[0].set_yticks(range(len(scen))); axes[0].set_yticklabels(scen[::-1])
    axes[0].legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, fontsize=5.8)
    _save(fig, out, "fig7_robustness")


# ------------------------------------------------------------------ Fig: dataset overview
def fig_dataset(meta_targets, out):
    fig, axes = plt.subplots(1, 3, figsize=(W2, 58 * MM), gridspec_kw={"wspace": 0.45})
    src = meta_targets.fold_name.value_counts()
    axes[0].barh(src.index[::-1], src.values[::-1], color=SLOTS[0], height=0.6)
    axes[0].set_xlabel("Target systems (count)")
    axes[1].hist(np.log10(meta_targets.n_compounds), bins=25, color=SLOTS[0], rwidth=0.9)
    axes[1].set_xlabel("Compounds per system (log10)"); axes[1].set_ylabel("Systems")
    axes[2].hist(meta_targets.rt_max, bins=25, color=SLOTS[0], rwidth=0.9)
    axes[2].set_xlabel("Last eluting compound (min)"); axes[2].set_ylabel("Systems")
    for ax, l in zip(axes, "abc"):
        _panel(ax, l); ax.grid(axis="x" if ax is axes[0] else "y")
    _save(fig, out, "figS1_dataset_overview")


# ------------------------------------------------------------------ calibration curve across alpha
def fig_alpha_calibration(summ_alpha, out, methods, labels, n=20):
    s = summ_alpha[summ_alpha.n == n]
    fig, ax = plt.subplots(figsize=(W1, 70 * MM))
    for i, m in enumerate(methods):
        d = s[s.method == m].sort_values("alpha")
        ax.plot(1 - d.alpha, d.coverage_mean, color=mcol(m, i), marker=mmark(m, i), ls=mls(m), ms=3.8, mec="white", mew=0.5, label=labels.get(m, m))
    lim = [0.75, 1.0]
    ax.plot(lim, lim, color=INK, lw=0.8, ls=(0, (3, 2))); ax.set_xlim(lim); ax.set_ylim(0.6, 1.0)
    ax.set_xlabel("Nominal coverage 1 − α"); ax.set_ylabel("Mean empirical coverage"); ax.grid(); ax.legend(frameon=False, fontsize=5.6)
    _save(fig, out, "figS2_alpha_calibration")


# ------------------------------------------------------------------ Fig 1 (revision): XSC workflow schematic (no data)
def fig_workflow_xsc(out):
    fig, ax = plt.subplots(figsize=(W2, 78 * MM)); ax.set_xlim(0, 100); ax.set_ylim(0, 46); ax.axis("off")
    fs_ = 6.1
    boxes = [
        (0.5, 25, 23, 14, "Public RT data:\nhundreds of chromatographic\nsystems (RepoRT)", "#eef4fc"),
        (27, 25, 27, 14, "Replay the laboratory's protocol\n(n calibrants, selection rule ρ)\non every public system s;\nmodels never saw s", "#eef4fc"),
        (57.5, 25, 42, 14, "Per system: smallest multiplier Vₛ giving\nretention ≥ 1 − ε on its other compounds\nc = ⌈(1−δ)(S+1)⌉-th smallest Vₛ\n⇒ P{laboratory retention ≥ 1 − ε} ≥ 1 − δ", "#eaf7f1"),
        (0.5, 3, 31, 15, "Laboratory: n calibrants (sᵢ, tᵢ)\nmonotone map g; LOO residuals rᵢ\nscale τ̂ shrunk towards the\nmeta-learned difficulty σ(x)", "#fdf0ea"),
        (35, 3, 30, 15, "Window\nC(x) = g(s(x)) ± c · τ̂ · σ(x)\nfinite for any n ≥ 3", "#f3f2ef"),
        (68.5, 3, 31, 15, "Candidate filtering in NTS:\nkeep a structure if observed RT ∈ C(x)\n(e.g. PubChemLite isomers)", "#f3f2ef"),
    ]
    for x, y, w, h, txt, fc in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.0", fc=fc, ec=INK2, lw=0.7))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs_, color=INK, linespacing=1.3)
    arr = dict(arrowstyle="-|>", mutation_scale=8, lw=0.9, color=INK2)
    for (x0, y0, x1, y1) in [(24.3, 32, 26.7, 32), (54.3, 32, 57.2, 32), (32.0, 10.5, 34.7, 10.5), (65.8, 10.5, 68.2, 10.5),
                             (78, 24.7, 55, 18.3)]:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), **arr))
    ax.add_patch(FancyArrowPatch((16, 18.3), (35, 24.7), **dict(arr, linestyle=(0, (3, 2)))))
    ax.text(14, 21.0, "protocol (n, ρ)", fontsize=5.8, color=INK2, rotation=17)
    ax.text(70, 21.2, "multiplier c", fontsize=5.8, color=INK2)
    ax.text(50, 41.5, "systems, not compounds, are the exchangeable calibration units", ha="center", fontsize=6.4, color=INK2, style="italic")
    _save(fig, out, "fig1_workflow")

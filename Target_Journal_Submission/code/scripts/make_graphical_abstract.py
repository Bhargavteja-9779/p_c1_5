"""Graphical abstract (Elsevier: 13 x 5 cm at 300 dpi = 1535 x 591 px, above the 1328 x 531 px minimum).
Numbers are read from tables/table2_main_n20_a010.csv. Output: figures/graphical_abstract.{png,pdf,tiff}
Usage: python scripts/make_graphical_abstract.py"""
import os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from conformrt.visualization.figures import METHOD_COLORS, INK, INK2, MUTED, GRID  # noqa: E402  (also sets rcParams)

t2 = pd.read_csv(os.path.join(ROOT, "tables", "table2_main_n20_a010.csv")).set_index("method")
ORDER = [("B3_gp_PI", "GP projection"), ("B2_pspline_PI", "P-spline"), ("B1_gauss_rmse", "Gaussian ±z·RMSE"),
         ("B5_jackknife_plus", "Jackknife+"), ("P", "Proposed")]
CM = 1 / 2.54
fig = plt.figure(figsize=(13 * CM, 5 * CM))
ax = fig.add_axes([0.0, 0.0, 0.57, 1.0]); ax.set_xlim(0, 57); ax.set_ylim(0, 25); ax.axis("off")
fs = 5.2
boxes = [(0.8, 14.4, 15.8, 8.4, "Public RT data\n(other laboratories'\nLC systems)", "#eef4fc"),
         (19.2, 14.4, 18.6, 8.4, "Retention-order score\n+ difficulty model\nlearned across systems", "#eef4fc"),
         (40.2, 14.4, 16.0, 8.4, "New LC system:\nn ≥ 1/α − 1\ncalibrants", "#fdf0ea"),
         (3.0, 1.8, 31.0, 8.6, "Monotone calibrant map +\nnormalised jackknife+  →  RT window\nfor every candidate structure", "#eaf7f1"),
         (37.5, 1.8, 18.7, 8.6, "Keep a candidate\nif its observed RT\nlies in the window", "#f3f2ef")]
for x, y, w, h, txt, fc in boxes:
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=0.8", fc=fc, ec=INK2, lw=0.6))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, color=INK, linespacing=1.25)
arr = dict(arrowstyle="-|>", mutation_scale=7, lw=0.8, color=INK2)
for (x0, y0, x1, y1) in [(17.0, 18.6, 18.8, 18.6), (28.5, 14.0, 22.0, 10.8), (48.1, 14.0, 30.0, 10.8), (34.4, 6.1, 37.1, 6.1)]:
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), **arr))

bx = fig.add_axes([0.765, 0.25, 0.215, 0.62])
for i, (m, lab) in enumerate(ORDER):
    r = t2.loc[m]; c = METHOD_COLORS[m] if m in ("P", "B5_jackknife_plus") else MUTED
    bx.plot([r.coverage_ci_lo, r.coverage_ci_hi], [i, i], color=c, lw=1.6)
    bx.plot([r.coverage_mean], [i], marker="o", ms=4.0, color=c, mec="white", mew=0.5)
bx.axvline(0.90, color=INK, lw=0.7, ls=(0, (3, 2)))
bx.set_yticks(range(len(ORDER)))
bx.set_yticklabels([f"{l}  {t2.loc[m].coverage_mean:.3f}" for m, l in ORDER], fontsize=5.2)
bx.set_xlim(0.80, 0.93); bx.set_ylim(-0.6, len(ORDER) - 0.4)
bx.set_xticks([0.80, 0.85, 0.90]); bx.tick_params(axis="x", labelsize=5.2); bx.grid(axis="x", color=GRID, lw=0.5)
bx.set_xlabel("Mean coverage (95 % CI)\n228 new RP systems, n = 20", fontsize=5.2)
bx.text(0.90, len(ORDER) - 0.45, "nominal", fontsize=4.8, color=INK, ha="center", va="bottom")
out = os.path.join(ROOT, "figures")
for ext in ("png", "tiff"):
    fig.savefig(os.path.join(out, f"graphical_abstract.{ext}"), dpi=300)
fig.savefig(os.path.join(out, "graphical_abstract.pdf"), metadata={"CreationDate": None})
print("graphical abstract written")

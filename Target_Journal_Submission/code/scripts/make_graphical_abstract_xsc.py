"""Graphical abstract (13 x 5 cm, 300 dpi) for the XSC article, drawn from results/xsc_metrics.json."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "src"))
from conformrt.visualization.figures import plt, INK, INK2, GRID
from matplotlib.patches import FancyBboxPatch
M = json.load(open(os.path.join(HERE, "..", "..", "results", "xsc_metrics.json")))
S = M["summary"]; C = M["case"]
rows = [("±z·RMSE", S["random|20|B1_gauss_rmse|0.1"]["fail"], C["B1_gauss_rmse"]["labfail"], "#eb6834"),
        ("±2·RMSE", M["extra"]["RMSE_x2"]["fail"], C["RMSE_x2"]["labfail"], "#e87ba4"),
        ("Jackknife+", S["random|20|B5_jackknife_plus|0.1"]["fail"], C["B5_jackknife_plus"]["labfail"], "#1baf7a"),
        ("XSC (this work)", S["random|20|XSC_sigs_lab|0.1"]["fail"], C["XSC_flat_lab"]["labfail"], "#2a78d6")]
fig = plt.figure(figsize=(13 / 2.54, 5 / 2.54))
ax0 = fig.add_axes([0.0, 0.0, 0.42, 1.0]); ax0.axis("off"); ax0.set_xlim(0, 10); ax0.set_ylim(0, 10)
for (y, txt, fc) in [(6.6, "Lab's protocol: n calibrants,\nselection rule", "#fdf0ea"), (3.6, "Replayed on hundreds of\npublic LC systems", "#eef4fc"),
                     (0.6, "RT window: ≥90 % of correct\nstructures kept in ≥90 % of labs", "#eaf7f1")]:
    ax0.add_patch(FancyBboxPatch((0.4, y), 9.0, 2.4, boxstyle="round,pad=0.15,rounding_size=0.4", fc=fc, ec=INK2, lw=0.6))
    ax0.text(4.9, y + 1.2, txt, ha="center", va="center", fontsize=6.2, color=INK)
for y in (6.5, 3.5):
    ax0.annotate("", xy=(4.9, y - 0.4), xytext=(4.9, y), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.8))
ax = fig.add_axes([0.53, 0.2, 0.44, 0.66])
import numpy as np
y = np.arange(len(rows)); h = 0.36
for i, (name, a, b, col) in enumerate(rows):
    ax.barh(i + h / 2, 100 * a, height=h * 0.9, color=col); ax.barh(i - h / 2, 100 * b, height=h * 0.9, color=col, alpha=0.5, hatch="///", edgecolor="white", lw=0)
    ax.text(100 * a + 1.5, i + h / 2, f"{100*a:.0f} %", va="center", fontsize=5.5, color=INK)
    ax.text(100 * b + 1.5, i - h / 2, f"{100*b:.0f} %", va="center", fontsize=5.5, color=INK)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=6); ax.invert_yaxis()
ax.axvline(10, color=INK2, lw=0.7, ls=(0, (3, 2))); ax.set_xlim(0, 75)
ax.set_xlabel("Laboratories keeping < 90 % of correct structures (%)", fontsize=5.8)
ax.tick_params(axis="x", labelsize=5.5); ax.grid(True, axis="x", color=GRID, lw=0.4); ax.set_axisbelow(True)
ax.set_title("solid: 228 held-out systems · hatched: 37 NORMAN labs,\nown mixes (XSC with σ ≡ 1)", fontsize=5.4, color=INK2, loc="left")
out = os.path.join(HERE, "..", "..", "figures")
for ext, kw in (("png", {"dpi": 300}), ("tiff", {"dpi": 300}), ("pdf", {"metadata": {"CreationDate": None}})):
    fig.savefig(os.path.join(out, f"graphical_abstract.{ext}"), **kw)
print("graphical abstract written")

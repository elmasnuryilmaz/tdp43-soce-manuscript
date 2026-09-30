#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Journal graphical abstract (Elsevier format: 5:2, at least 1328 x 531 px at 300 dpi, large type).

The detailed four-panel summary of the repository is figures/graphical_abstract.* (code/fig_graphical_abstract.py).
This is the simplified version for the journal: one message per column, left to right, few words, type sized so
that it stays legible when ScienceDirect scales the image to a 500 x 200 pixel window. Every value is read from
files of this repository:

  A  RT-qPCR fold changes and the three Fura-2 wells   supplementary/S1_laboratory_source_data.xlsx
  B  CBARP exon-4 junction usage per library            source_data/CBARP_locus/CBARP_junction_b_share_per_library.tsv
     number of datasets with CBARP events               supplementary/S3_rMATS_significant_events.csv.gz
  C  ALS tissue, one dot per brain or cord region       tables/Table5_cross_disease_comparison.csv

No regional effect is averaged and no diseases are ranked. Colour: blue = control, orange = knockdown.
Outputs: figures/graphical_abstract_journal.{png,tiff,pdf,svg}

Run from the repository root:  /usr/bin/python3 code/fig_graphical_abstract_journal.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrow, FancyBboxPatch
from PIL import Image

P = Path(__file__).resolve().parents[1]
OUT = P / "figures" / "graphical_abstract_journal"

W, H = 13.0, 5.2                       # inches; 13 x 5.2 = 5:2, 3900 x 1560 px at 300 dpi
BLUE, ORANGE, INK, GREY, PALE = "#0072B2", "#D55E00", "#1f1f1f", "#5c5c5c", "#f0f0f0"

plt.rcParams.update({
    "font.family": "Arial", "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.facecolor": "white",
    "figure.facecolor": "white", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

# ------------------------------------------------------------------ source data
wb = P / "supplementary" / "S1_laboratory_source_data.xlsx"
rel = pd.read_excel(wb, sheet_name="Target_qPCR_rel")
fura = pd.read_excel(wb, sheet_name="Fura2_per_culture")
NT = "Non-targeting shRNA control"
GENES = ["TRPC1", "STIM1", "ORAI1", "ATP2A3"]
share = pd.read_csv(P / "source_data" / "CBARP_locus" / "CBARP_junction_b_share_per_library.tsv", sep="\t")
t5 = pd.read_csv(P / "tables" / "Table5_cross_disease_comparison.csv")
t5 = t5[(t5.group == "ALS") & t5.gene.isin(["TRPC1", "SARAF", "CBARP"])]

s3 = pd.read_csv(P / "supplementary" / "S3_rMATS_significant_events.csv.gz")
cb = s3[(s3.gene.astype(str).str.upper() == "CBARP") & (s3.model == "JC") & (s3.passes_coverage_filter == True)]  # noqa: E712
N_DATASETS_WITH_CBARP, N_DATASETS = cb.dataset.nunique(), s3.dataset.nunique()
assert (N_DATASETS_WITH_CBARP, N_DATASETS) == (5, 6), (N_DATASETS_WITH_CBARP, N_DATASETS)

fig = plt.figure(figsize=(W, H))
bg = fig.add_axes([0, 0, 1, 1])
bg.set_xlim(0, W)
bg.set_ylim(0, H)
bg.axis("off")
renderer = fig.canvas.get_renderer()


def axes_in(x, y, w, h, **kw):
    """Axes placed in inch coordinates of the canvas."""
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def text(x, y, s, size, weight="normal", color=INK, ha="center", va="center", style="normal", max_w=None, floor=12.5, **kw):
    """Text on the canvas; the size is reduced (never below ``floor``) until it fits in ``max_w`` inches."""
    t = bg.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va, fontstyle=style, **kw)
    if max_w:
        while size > floor:
            bb = t.get_window_extent(renderer)
            if bb.width / fig.dpi <= max_w:
                break
            size -= 0.25
            t.set_fontsize(size)
    return t


COLS = {"A": (0.30, 3.80), "B": (4.60, 3.80), "C": (8.90, 3.80)}    # x, width (inches)
GX = {k: v[0] + v[1] / 2 for k, v in COLS.items()}

# ------------------------------------------------------------------ headers (one line each)
HEAD_Y, TITLE_Y, CAP_Y = 4.95, 4.52, (1.36, 1.10)
for k, h1 in (("A", "TDP-43 knockdown cells"), ("B", "Public RNA-seq"), ("C", "ALS tissue")):
    text(GX[k], HEAD_Y, h1, 22, "bold", max_w=3.7)

# flow arrows between the columns
for x in (4.20, 8.50):
    bg.add_patch(FancyArrow(x, 3.05, 0.30, 0, width=0.20, head_width=0.46, head_length=0.20,
                            length_includes_head=True, color="#c8c8c8", lw=0))


def caption(col, line1, line2):
    text(GX[col], CAP_Y[0], line1, 14.5, color=GREY, max_w=3.75)
    text(GX[col], CAP_Y[1], line2, 14.5, color=GREY, max_w=3.75)


# ------------------------------------------------------------------ A1: RT-qPCR
text(GX["A"], TITLE_Y, "mRNA, RT-qPCR, n = 4", 16, "bold", max_w=3.75)
ax = axes_in(0.45, 3.22, 3.50, 0.84)
folds = []
for i, g in enumerate(GENES):
    c = rel[(rel.gene == g) & (rel.group == NT)].rel_expression.to_numpy(float)
    k = rel[(rel.gene == g) & (rel.group == "shTDP-43")].rel_expression.to_numpy(float)
    assert len(c) == len(k) == 4
    fold = k.mean() / c.mean()
    folds.append(fold)
    ax.bar(i - 0.17, 1.0, 0.30, color=BLUE)
    ax.bar(i + 0.17, fold, 0.30, color=ORANGE)
    ax.text(i + 0.17, fold + 0.12, f"{fold:.1f}×", ha="center", va="bottom", fontsize=17, fontweight="bold", color=INK)
    ax.text(i, -0.30, g, ha="center", va="top", fontsize=14.5, fontstyle="italic", color=INK)
assert [round(f, 1) for f in folds] == [1.8, 1.9, 1.7, 3.2], folds
ax.set_xlim(-0.55, 3.55)
ax.set_ylim(0, 3.45)
ax.axis("off")

# ------------------------------------------------------------------ A2: Fura-2, one plate
text(GX["A"], 2.62, r"$\mathbf{Ca^{2+}}$ readdition", 16, "bold", max_w=3.75)
rows = [(NT, "control", BLUE), ("shTDP-43", "knockdown", ORANGE)]
ax = axes_in(1.52, 1.52, 2.40, 0.86)
for j, (grp, lab, colr) in enumerate(rows):
    v = fura.loc[fura.group == grp, "readdition"].to_numpy(float)
    assert len(v) == 3
    y = 1 - j
    ax.barh(y, v.mean(), 0.62, color=colr)
    ax.scatter(v, y + np.array([-0.16, 0.0, 0.16]), s=44, facecolor="white", edgecolor=INK, lw=1.1, zorder=3)
    bg.text(1.42, 1.52 + 0.86 * (y + 0.5) / 2.0, lab, ha="right", va="center", fontsize=15, fontweight="bold", color=colr)
ax.set_xlim(0, 2.3)
ax.set_ylim(-0.5, 1.5)
ax.axis("off")

# ------------------------------------------------------------------ B: CBARP exon-4 junction usage
text(GX["B"], TITLE_Y, "CBARP splicing changed", 17, "bold", max_w=3.75)
text(GX["B"], TITLE_Y - 0.27, f"in {N_DATASETS_WITH_CBARP} of {N_DATASETS} datasets", 17, "bold", max_w=3.75)
ax = axes_in(4.85, 2.22, 3.30, 1.62)
xs = {("SH-SY5Y", "Control"): 0.0, ("SH-SY5Y", "TDP-43 KD"): 1.0, ("iPSC colonies", "Control"): 2.5, ("iPSC colonies", "TDP-43 KD"): 3.5}
means = {}
for (ds, grp), x0 in xs.items():
    v = share[(share.dataset == ds) & (share.group == grp)].share_b.to_numpy(float)
    colr = BLUE if grp == "Control" else ORANGE
    m = float(v.mean())
    means[(ds, grp)] = m
    ax.bar(x0, m, 0.78, color=colr)
    ax.scatter(x0 + np.linspace(-0.17, 0.17, len(v)), v, s=40, facecolor="white", edgecolor=INK, lw=1.0, zorder=3)
    ax.text(x0, min(m, 1.0) + 0.11, f"{100 * m:.0f}%", ha="center", va="bottom", fontsize=16, fontweight="bold", color=INK)
    ax.text(x0, -0.05, "control" if grp == "Control" else "KD", ha="center", va="top", fontsize=13.5, color=colr, fontweight="bold")
for xc, lab in ((0.5, "SH-SY5Y"), (3.0, "iPSC colonies")):
    ax.text(xc, -0.30, lab, ha="center", va="top", fontsize=15, color=INK)
ax.set_xlim(-0.55, 4.05)
ax.set_ylim(0, 1.25)
ax.axis("off")
assert round(means[("SH-SY5Y", "Control")], 2) == 0.17 and round(means[("SH-SY5Y", "TDP-43 KD")], 2) == 0.86
assert round(means[("iPSC colonies", "Control")], 2) == 0.01 and round(means[("iPSC colonies", "TDP-43 KD")], 2) == 0.78
caption("B", "share of exon-4 reads using the", "alternative 3′ splice site")

# ------------------------------------------------------------------ C: ALS tissue, one dot per region
text(GX["C"], TITLE_Y, "TRPC1, SARAF up; CBARP down", 16, "bold", max_w=3.75)
ax = axes_in(9.00, 2.12, 3.05, 1.90)
genes_c = ["TRPC1", "SARAF", "CBARP"]
for xi, g in enumerate(genes_c):
    rows_g = t5[t5.gene == g].sort_values("cliffs_delta")
    assert len(rows_g) == 10, (g, len(rows_g))
    for j, (_, r) in enumerate(rows_g.iterrows()):
        sig = r.q_value < 0.05
        ax.scatter(xi + (j % 3 - 1) * 0.20, r.cliffs_delta, s=95, facecolor=INK if sig else "white",
                   edgecolor=INK, lw=1.3, zorder=3)
    ax.text(xi, -1.02, g, ha="center", va="top", fontsize=16, fontstyle="italic", fontweight="bold", color=INK)
ax.axhline(0, color=GREY, lw=1.3, ls=(0, (4, 3)))
ax.set_xlim(-0.55, 2.55)
ax.set_ylim(-0.95, 1.10)
ax.axis("off")
bg.text(12.62, 3.98, "↑ higher\nin ALS", ha="right", va="top", fontsize=14, color=GREY, linespacing=1.05)
bg.text(12.62, 2.38, "↓ lower", ha="right", va="bottom", fontsize=14, color=GREY)
caption("C", "one dot per region, Cliff's δ;", "filled dots: q < 0.05")

# ------------------------------------------------------------------ bottom line
bg.add_patch(FancyBboxPatch((0.30, 0.14), W - 0.60, 0.66, boxstyle="round,pad=0,rounding_size=0.12",
                            facecolor=PALE, edgecolor="none"))
text(W / 2, 0.47, "TDP-43 loss: calcium-regulatory RNA candidates", 22, "bold", max_w=12.0)

fig.savefig(str(OUT) + ".png", dpi=300)
fig.savefig(str(OUT) + ".pdf")
fig.savefig(str(OUT) + ".svg")
plt.close(fig)

im = Image.open(str(OUT) + ".png").convert("RGB")      # no alpha channel in the files sent to the journal
im.save(str(OUT) + ".png", dpi=(300, 300))
im.save(str(OUT) + ".tiff", dpi=(300, 300), compression="tiff_lzw")
assert im.size == (3900, 1560), im.size
assert im.size[0] >= 1328 and im.size[1] >= 531 and abs(im.size[0] / im.size[1] - 2.5) < 1e-9
print("written:", str(OUT) + ".{png,tiff,pdf,svg}", im.size)

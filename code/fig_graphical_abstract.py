#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Graphical abstract (repository summary and optional journal abstract image), release v1.0.8.

Every value is read from files in this repository, not typed in:
  A  RT-qPCR fold changes and the Fura-2 wells   supplementary/S1_laboratory_source_data.xlsx
  B  change in transcript-family composition     tables/Table1_transcript_family_abundance.csv
  C  CBARP exon-4 junction usage, per library    source_data/CBARP_locus/CBARP_junction_b_share_per_library.tsv
  D  ALS tissue, one dot per brain or cord region tables/Table5_cross_disease_comparison.csv

Panel D shows each region's own Cliff's delta and its own Benjamini-Hochberg decision from Table 5.
No averaged effect and no disease ranking is drawn: averaging regional deltas would not be a
pooled effect or a tested quantity.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[1]
OUT = P / "figures" / "graphical_abstract"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 9.5, "axes.titleweight": "bold", "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight", "savefig.facecolor": "white",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})
BLUE, ORANGE, TEAL, RED, GREY, INK = "#0072B2", "#D55E00", "#0d6259", "#c0392b", "#7f8c8d", "#222222"

# ------------------------------------------------------------------ source data
wb = P / "supplementary" / "S1_laboratory_source_data.xlsx"
rel = pd.read_excel(wb, sheet_name="Target_qPCR_rel")
fura = pd.read_excel(wb, sheet_name="Fura2_per_culture")
NT = "Non-targeting shRNA control"
genes = ["TRPC1", "STIM1", "ORAI1", "ATP2A3"]
t1 = pd.read_csv(P / "tables" / "Table1_transcript_family_abundance.csv")
t1 = t1[t1.Gene != "FAMILY TOTAL"].set_index("Gene")
share = pd.read_csv(P / "source_data" / "CBARP_locus" / "CBARP_junction_b_share_per_library.tsv", sep="\t")
t5 = pd.read_csv(P / "tables" / "Table5_cross_disease_comparison.csv")
t5 = t5[(t5.group == "ALS") & t5.gene.isin(["TRPC1", "SARAF", "CBARP"])]

fig = plt.figure(figsize=(15.5, 5.6))
gs = fig.add_gridspec(1, 4, width_ratios=[1.05, 1.0, 0.85, 1.2], wspace=0.55)

# ---------------------------------------------------- A: laboratory measurements
sub = gs[0, 0].subgridspec(2, 1, height_ratios=[1.15, 1.0], hspace=0.62)
ax = fig.add_subplot(sub[0])
for i, g in enumerate(genes):
    c = rel[(rel.gene == g) & (rel.group == NT)].rel_expression.to_numpy()
    k = rel[(rel.gene == g) & (rel.group == "shTDP-43")].rel_expression.to_numpy()
    fold = k.mean() / c.mean()
    ax.bar(i, fold, 0.6, color=RED, alpha=0.85)
    ax.scatter(i + np.linspace(-0.15, 0.15, len(k)), k / c.mean(), s=14, facecolor="white", edgecolor=INK, lw=0.7, zorder=3)
    ax.text(i, 4.05, f"{fold:.1f}×", ha="center", fontsize=8.5, weight="bold", color="#333")
ax.axhline(1, color="#666", lw=1, ls="--")
ax.set_xticks(range(4), genes, fontsize=8.5)
ax.set_ylim(0, 4.6)
ax.set_ylabel("mRNA, knockdown / control", fontsize=8.5)
ax.set_title("A · RT-qPCR, 4 biological replicates/group", loc="left", fontsize=9)
bx = fig.add_subplot(sub[1])
groups = [NT, "shTDP-43"]
for i, (g, colr) in enumerate(zip(groups, (BLUE, RED))):
    v = fura.loc[fura.group == g, "readdition"].to_numpy(float)
    bx.bar(i, v.mean(), 0.55, color=colr, alpha=0.85)
    bx.scatter(i + np.linspace(-0.12, 0.12, 3), v, s=16, facecolor="white", edgecolor=INK, lw=0.7, zorder=3)
bx.set_xticks([0, 1], ["non-targeting", "shTDP-43"], fontsize=8.5)
bx.set_ylabel("Ca$^{2+}$ readdition, Δ(F340/F380)", fontsize=8)
bx.set_ylim(0, 2.3)
bx.set_title("Fura-2, one plate: 3 wells/group\n(descriptive; not paired with RT-qPCR)", loc="left", fontsize=8.5)

# ---------------------------------------------------- B: composition within families
ax = fig.add_subplot(gs[0, 1])
sel = ["ATP2A2", "MCU", "MCUB", "ORAI2", "ORAI1", "ATP2A3", "ORAI3", "STIM1", "SARAF"]
d = [t1.loc[s, "delta_TPM_adjusted"] for s in sel]
order = np.argsort(d)
sel = [sel[i] for i in order]
d = [d[i] for i in order]
ax.barh(np.arange(len(sel)), d, color=[RED if v < 0 else TEAL for v in d], alpha=0.9)
ax.axvline(0, color="#333", lw=1)
ax.set_yticks(np.arange(len(sel)), [f"{s}  ({t1.loc[s, 'share_of_family_control_pct']:.0f}%)" for s in sel], fontsize=8.2)
ax.set_xlabel("change in transcript abundance (ΔTPM, adjusted\nfor library composition); % = share of its family\nin control cells", fontsize=8)
ax.set_xlim(-30, 88)
ax.set_title("B · Public SH-SY5Y RNA-seq: family shifts", loc="left", fontsize=9)

# ---------------------------------------------------- C: CBARP junction usage
ax = fig.add_subplot(gs[0, 2])
xs = {("SH-SY5Y", "Control"): 0, ("SH-SY5Y", "TDP-43 KD"): 1, ("iPSC colonies", "Control"): 2.6, ("iPSC colonies", "TDP-43 KD"): 3.6}
for (ds, grp), x0 in xs.items():
    v = share[(share.dataset == ds) & (share.group == grp)].share_b.to_numpy(float)
    colr = BLUE if grp == "Control" else ORANGE
    ax.scatter(x0 + np.linspace(-0.14, 0.14, len(v)), v, s=22, facecolor="white" if grp == "Control" else colr, edgecolor=colr if grp == "Control" else INK, lw=0.8, zorder=3)
    ax.plot([x0 - 0.28, x0 + 0.28], [v.mean()] * 2, color=colr, lw=1.8)
ax.set_xticks(list(xs.values()), ["Control", "KD", "Control", "KD"], fontsize=8.5)
for xc, lab in ((0.5, "SH-SY5Y"), (3.1, "iPSC colonies")):
    ax.text(xc, -0.085, lab, ha="center", va="top", fontsize=8.3, transform=ax.get_xaxis_transform())
ax.set_ylim(-0.03, 1.08)
ax.set_ylabel("share of CBARP exon-4 reads using\nthe alternative 3′ splice site", fontsize=8.3)
ax.set_title("C · CBARP splicing, two human models", loc="left", fontsize=9)

# ---------------------------------------------------- D: ALS tissue, one dot per region
ax = fig.add_subplot(gs[0, 3])
gene_rows = ["TRPC1", "SARAF", "CBARP"]
for yy, g in zip(range(len(gene_rows))[::-1], gene_rows):
    rows = t5[t5.gene == g].sort_values("cliffs_delta")
    for j, (_, r) in enumerate(rows.iterrows()):
        sig = r.q_value < 0.05
        colr = ORANGE if r.cliffs_delta > 0 else BLUE
        ax.scatter(r.cliffs_delta, yy + (j % 3 - 1) * 0.11, s=34, facecolor=colr if sig else "white", edgecolor=colr, lw=1.1, zorder=3)
ax.axvline(0, color="#333", lw=1)
ax.set_yticks(range(len(gene_rows))[::-1], gene_rows, fontsize=9, weight="bold")
ax.set_xlim(-0.85, 1.05)
ax.set_ylim(-0.6, 3.3)
ax.set_xlabel("Cliff's δ, ALS vs controls (one dot per region)")
ax.set_title("D · ALS tissue: regional differences", loc="left", fontsize=9)
ax.scatter([], [], s=34, facecolor=ORANGE, edgecolor=ORANGE, label="q < 0.05, higher in ALS")
ax.scatter([], [], s=34, facecolor=BLUE, edgecolor=BLUE, label="q < 0.05, lower in ALS")
ax.scatter([], [], s=34, facecolor="white", edgecolor=GREY, label="not significant")
ax.legend(frameon=False, fontsize=7.6, loc="upper left")

fig.suptitle("Calcium-regulatory RNA candidates after TDP-43 knockdown, and a one-plate Ca$^{2+}$-readdition observation",
             y=1.03, fontsize=12.5, weight="bold")
fig.text(0.5, -0.06,
         "The RNA analyses identify candidates for independently replicated experiments; they do not establish a mechanism for the Fura-2 observation.",
         ha="center", fontsize=9.3, color="#444")
fig.savefig(str(OUT) + ".png")
fig.savefig(str(OUT) + ".pdf")
fig.savefig(str(OUT) + ".svg")
plt.close(fig)
print("written:", str(OUT) + ".png")

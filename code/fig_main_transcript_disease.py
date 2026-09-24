#!/usr/bin/env python3
"""Main Figures 2 and 3, regenerated from Tables 4 and 5."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

PKG = Path("/Users/elmas/Desktop/MAKALE/11_NEUROCHEMISTRY_INTERNATIONAL_FIGURE_REVISION")
OUT = PKG / "figures" / "main"
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE = "#0072B2", "#D55E00"
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.labelsize": 8,
    "axes.titlesize": 9, "axes.titleweight": "bold",
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .7, "xtick.major.width": .7, "ytick.major.width": .7,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "savefig.dpi": 300, "figure.facecolor": "white", "savefig.facecolor": "white",
})

def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png")
    fig.savefig(OUT / f"{stem}.svg")
    plt.close(fig)

# Figure 2: adjusted TPM and DESeq2 gene-level change from six public libraries.
t = pd.read_csv(PKG / "tables" / "Table1_transcript_family_abundance.csv").set_index("Gene")
groups = [
    ["ORAI2", "ORAI1", "ORAI3"],
    ["STIM1", "STIM2"],
    ["TRPC1"],
    ["SARAF", "CBARP", "STIMATE"],
    ["ATP2A2", "ATP2A3"],
    ["MCU", "MICU2"],
]
genes, yy = [], []
y = 0
for group in groups:
    for g in group:
        genes.append(g)
        yy.append(y)
        y += 1
    y += .45
fig, (ax_a, ax_b) = plt.subplots(
    1, 2, figsize=(7.1, 5.8), sharey=True, layout="constrained",
    gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": .07})
for g, pos in zip(genes, yy):
    r = t.loc[g]
    c, k = float(r.TPM_control_adjusted), float(r.TPM_KD_adjusted)
    ax_a.plot([c, k], [pos, pos], color="#AAB2B7", lw=1.1, zorder=1)
    ax_a.scatter(c, pos, s=30, color=BLUE, zorder=3)
    ax_a.scatter(k, pos, s=30, color=ORANGE, zorder=3)
    sig = pd.notna(r.padj) and r.padj < .05
    ax_b.scatter(float(r.log2FC), pos, s=33, facecolor="#333333" if sig else "white",
                 edgecolor="#333333", linewidth=.9, zorder=3)
ax_a.set_yticks(yy, genes)
ax_a.invert_yaxis()
for tick in ax_a.get_yticklabels():
    if tick.get_text() in {"TRPC1", "STIM1", "ORAI1", "ATP2A3"}:
        tick.set_fontweight("bold")
ax_a.set_xscale("log")
ax_a.set_xlim(.8, 350)
ax_a.set_xlabel("Adjusted TPM (log scale)")
ax_a.set_title("A  Transcript abundance", loc="left")
ax_a.scatter([], [], color=BLUE, s=30, label="Control")
ax_a.scatter([], [], color=ORANGE, s=30, label="TDP-43 depletion")
ax_a.legend(loc="lower right", frameon=False, fontsize=7.5)
ax_b.axvline(0, color="#444444", lw=.8)
ax_b.set_xlim(-1.9, 2.3)
ax_b.set_xlabel("DESeq2 log$_2$ fold change")
ax_b.set_title("B  Differential expression", loc="left")
ax_b.tick_params(axis="y", left=False, labelleft=False)
ax_b.scatter([], [], s=30, facecolor="#333333", edgecolor="#333333", label="q < 0.05")
ax_b.scatter([], [], s=30, facecolor="white", edgecolor="#333333", label="q ≥ 0.05")
ax_b.legend(loc="lower right", frameon=False, fontsize=7.5)
save(fig, "Figure2_transcript_profile")

# Figure 3: ALS-only effect-size matrix. Other diseases are in Supplementary S7.
d = pd.read_csv(PKG / "tables" / "Table5_cross_disease_comparison.csv")
d = d[(d.group == "ALS") & d.gene.isin(["TRPC1", "SARAF", "CBARP"])].copy()
regions = [
    "Hippocampus", "Cerebellum", "Motor cortex (lateral)", "Frontal cortex",
    "Motor cortex (medial)", "Temporal cortex", "Occipital cortex",
    "Spinal cord (lumbar)", "Spinal cord (cervical)", "Spinal cord (thoracic)",
]
targets = ["TRPC1", "SARAF", "CBARP"]
effect = d.pivot(index="region", columns="gene", values="cliffs_delta").loc[regions, targets]
qval = d.pivot(index="region", columns="gene", values="q_value").loc[regions, targets]
assert effect.notna().all().all() and qval.notna().all().all()
fig, ax = plt.subplots(figsize=(7.1, 5.5), layout="constrained")
cmap = LinearSegmentedColormap.from_list("accessible_diverging",
                                         ["#95C4DF", "#FFFFFF", "#F1B481"])
im = ax.imshow(effect.to_numpy(), vmin=-1, vmax=1, cmap=cmap, aspect="auto")
for i, region in enumerate(regions):
    for j, gene in enumerate(targets):
        v = effect.loc[region, gene]
        mark = "*" if qval.loc[region, gene] < .05 else ""
        ax.text(j, i, f"{v:+.2f}{mark}", ha="center", va="center",
                fontsize=7.8, color="black", weight="bold" if mark else "normal")
ax.set_xticks(range(3), targets, fontsize=8, weight="bold")
ax.set_yticks(range(len(regions)), regions)
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, length=0)
ax.axhline(6.5, color="#555555", lw=.65)
ax.set_title("ALS tissue versus non-neurological controls", loc="left", pad=10)
cb = fig.colorbar(im, ax=ax, fraction=.042, pad=.035, ticks=[-1, -.5, 0, .5, 1])
cb.set_label("Cliff's δ (case − control)")
save(fig, "Figure3_ALS_expression")
print("Main Figures 2 and 3 written")

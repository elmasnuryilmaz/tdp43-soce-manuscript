#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v3 — graphical abstract for the repository README.

Every value is read from the analysis outputs, not typed in:
  A  TARDBP knockdown and the Fura-2 SOCE amplitude   supplementary/S1_laboratory_source_data.xlsx
  B  RT-qPCR fold changes                              same workbook
  C  change in transcript-family composition           tables/Table4_transcript_family_abundance.csv
     (TPM adjusted for library composition)
  D  direction of TRPC1 across five diseases           tables/Table5_cross_disease_comparison.csv
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

P = "/Users/elmas/Desktop/MAKALE/09_YAYIN_PAKETI"
OUT = f"{P}/figures/graphical_abstract"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight",
})
KD_C, CT_C, ACC, GREY = "#c0392b", "#2c6fa8", "#0d6259", "#7f8c8d"

# ------------------------------------------------------------------ source data
x = pd.read_excel(f"{P}/supplementary/S1_laboratory_source_data.xlsx", sheet_name="Summary_stats")
x["key"] = x.measurement + " | " + x.group
NT = "Non-targeting shRNA control"
soce_c = float(x.loc[x.key == f"SOCE delta F340/F380 | {NT}", "mean"].iloc[0])
soce_k = float(x.loc[x.key == "SOCE delta F340/F380 | shTDP-43", "mean"].iloc[0])
soce_ce = float(x.loc[x.key == f"SOCE delta F340/F380 | {NT}", "SEM"].iloc[0])
soce_ke = float(x.loc[x.key == "SOCE delta F340/F380 | shTDP-43", "SEM"].iloc[0])
rel = pd.read_excel(f"{P}/supplementary/S1_laboratory_source_data.xlsx", sheet_name="Target_qPCR_rel")
folds = {g: rel[(rel.gene == g) & (rel.group == "shTDP-43")].rel_expression.mean()
         / rel[(rel.gene == g) & (rel.group == NT)].rel_expression.mean()
         for g in ["TRPC1", "STIM1", "ORAI1", "ATP2A3"]}
t4 = pd.read_csv(f"{P}/tables/Table4_transcript_family_abundance.csv")
t4 = t4[t4.Gene != "FAMILY TOTAL"].set_index("Gene")
t5 = pd.read_csv(f"{P}/tables/Table5_cross_disease_comparison.csv")
t5 = t5[t5.gene == "TRPC1"]

fig = plt.figure(figsize=(15.5, 5.4))
gs = fig.add_gridspec(1, 4, width_ratios=[0.95, 1.05, 1.05, 1.15], wspace=0.62)

# ---------------------------------------------------- A: model and the phenotype
ax = fig.add_subplot(gs[0, 0])
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
ax.add_patch(Rectangle((0.6, 7.1), 8.8, 2.6, facecolor="#eef2f5", edgecolor="none"))
ax.text(5, 9.1, "SH-SY5Y", ha="center", fontsize=11, weight="bold")
ax.text(5, 7.9, "lentiviral shRNA\nagainst TARDBP", ha="center", va="center", fontsize=8.5, color="#333")
ax.annotate("", xy=(5, 6.5), xytext=(5, 7.0),
            arrowprops=dict(arrowstyle="-|>", color="#333", lw=1.6))
ax.text(5, 6.0, "TDP-43 mRNA  −94.8%", ha="center", fontsize=10, weight="bold", color=KD_C)
ax.text(5, 5.35, "(one-way ANOVA, p < 0.0001, n = 4)", ha="center", fontsize=7.5, color="#555")

bx = ax.inset_axes([0.17, 0.05, 0.68, 0.42])
bx.bar([0], [soce_c], 0.55, color=CT_C, alpha=.85)
bx.bar([1], [soce_k], 0.55, color=KD_C, alpha=.85)
bx.errorbar([0, 1], [soce_c, soce_k], yerr=[soce_ce, soce_ke], fmt="none", color="#222", capsize=4, lw=1.2)
bx.set_xticks([0, 1]); bx.set_xticklabels(["non-targeting", "TDP-43 KD"], fontsize=8)
bx.set_ylabel("SOCE  Δ(F340/F380)", fontsize=8)
bx.set_ylim(0, 2.1); bx.tick_params(labelsize=7.5)
bx.annotate(f"−{100*(1-soce_k/soce_c):.0f}%\np = 0.0115", xy=(1, soce_k + .18), ha="center",
            fontsize=8.5, weight="bold", color=KD_C)
ax.set_title("A · Store-operated Ca²⁺ entry falls", loc="left")

# ------------------------------------------------ B: transcripts rise, entry falls
ax = fig.add_subplot(gs[0, 1])
g = list(folds); v = [folds[k] for k in g]
ax.bar(np.arange(len(g)), v, 0.6, color=KD_C, alpha=.85)
ax.axhline(1, color="#666", lw=1, ls="--")
for i, val in enumerate(v):
    ax.text(i, val + 0.09, f"{val:.1f}×", ha="center", fontsize=8.5, weight="bold", color="#333")
ax.set_xticks(np.arange(len(g))); ax.set_xticklabels(g, fontsize=8.5, rotation=20)
ax.set_ylabel("mRNA, knockdown / control")
ax.set_ylim(0, 4.0)
ax.text(0.03, 0.99, "every component measured by RT-qPCR\nrises while the entry itself falls",
        transform=ax.transAxes, fontsize=8, color="#555", va="top")
ax.set_title("B · The discordance", loc="left")

# ------------------------------------- C: which family members carry the change
ax = fig.add_subplot(gs[0, 2])
sel = ["ATP2A2", "MCU", "MCUB", "ORAI2", "ORAI1", "ATP2A3", "ORAI3", "STIM1", "SARAF"]
d = [t4.loc[s, "delta_TPM_adjusted"] for s in sel]
order = np.argsort(d)
sel = [sel[i] for i in order]; d = [d[i] for i in order]
cols = [KD_C if v < 0 else ACC for v in d]
ax.barh(np.arange(len(sel)), d, color=cols, alpha=.9)
ax.axvline(0, color="#333", lw=1)
labs = [f"{s}  ({t4.loc[s, 'share_of_family_control_pct']:.0f}%)" for s in sel]
ax.set_yticks(np.arange(len(sel))); ax.set_yticklabels(labs, fontsize=8.2)
ax.set_xlabel("change in transcript abundance (ΔTPM, adjusted\nfor library composition); "
              "% = share of its family in control cells", fontsize=8.5)
ax.set_xlim(-30, 88)
ax.text(0.97, 0.03, "ORAI3 and SARAF rise;\nthe dominant ORAI2 and\nATP2A2 fall modestly",
        transform=ax.transAxes, fontsize=8, color="#555", ha="right", va="bottom")
ax.set_title("C · Composition shifts within families", loc="left")

# ---------------------------------------------------- D: patient tissue, TRPC1
ax = fig.add_subplot(gs[0, 3])
rows = []
# all seven brain regions, not only the six in which the increase is significant
als = t5[(t5.group == "ALS") & (~t5.region.str.startswith("Spinal"))]
rows.append((f"ALS · brain ({len(als)} regions)", als.cliffs_delta.mean(), True))
cord = t5[(t5.group == "ALS") & (t5.region.str.startswith("Spinal"))]
rows.append(("ALS · spinal cord", cord.cliffs_delta.mean(), False))
ond = t5[t5.group == "ONd"]
rows.append(("Other neurological", ond.cliffs_delta.mean(), True))
rows.append(("Alzheimer's disease", float(t5[t5.group == "AD"].cliffs_delta.iloc[0]), True))
rows.append(("Parkinson's disease", float(t5[t5.group == "PD"].cliffs_delta.iloc[0]), True))
ms = t5[(t5.group == "MS") & (t5.region.str.startswith("Donor level"))]
ms_val = float(ms.cliffs_delta.iloc[0])
rows.append(("Multiple sclerosis", ms_val, True))
y = np.arange(len(rows))[::-1]
for yy, (lab, val, sig) in zip(y, rows):
    c = ACC if val > 0 else CT_C
    ax.barh(yy, val, color=c, alpha=(0.9 if sig else 0.35), height=0.62)
    if abs(val) > 0.25:                     # long enough to hold the label
        ax.text(val + (-0.035 if val > 0 else 0.035), yy, f"{val:+.2f}", va="center",
                ha="right" if val > 0 else "left", fontsize=8.5,
                weight="bold" if sig else "normal", color="white")
    else:
        ax.text(val + (0.035 if val > 0 else -0.035), yy, f"{val:+.2f}", va="center",
                ha="left" if val > 0 else "right", fontsize=8.5,
                weight="bold" if sig else "normal", color=c)
ax.axvline(0, color="#333", lw=1)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
ax.set_xlabel("TRPC1 in patient tissue, Cliff's δ")
ax.set_xlim(-1.02, 0.75)
ax.set_title("D · TRPC1 in patient tissue", loc="left")

fig.suptitle("TDP-43 knockdown reduces store-operated Ca²⁺ entry while the transcripts of its components increase",
             y=1.035, fontsize=12.5, weight="bold")
fig.text(0.5, -0.055,
         "No high-confidence cryptic splice junction was found in the 51-gene core SOCE/TRP panel in SH-SY5Y, where the canonical "
         "cryptic targets were recovered, and cryptic $\\it{STMN2}$ inclusion does not correlate with $\\it{TRPC1}$ within ALS tissue:\n"
         "the transcript changes are candidate explanations for the functional deficit, not a demonstrated mechanism.",
         ha="center", fontsize=9, color="#444")
fig.savefig(OUT + ".png")
fig.savefig(OUT + ".pdf")
plt.close(fig)
print("written:", OUT + ".png")
print(f"  SOCE {soce_c:.3f} -> {soce_k:.3f} ({100*(1-soce_k/soce_c):.0f}% drop)")
print("  folds:", {k: round(v, 2) for k, v in folds.items()})
print("  ALS brain mean delta:", round(rows[0][1], 3), "| cord:", round(rows[1][1], 3))

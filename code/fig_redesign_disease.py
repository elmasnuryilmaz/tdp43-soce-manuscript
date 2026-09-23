#!/usr/bin/env python3
"""v3 — Figure 8: direction of TRPC1 change across five diseases.

ALS values: NYGC genome-wide per-region analysis (seven brain regions and three
spinal cord levels). MS values: sample and donor level for GSE138614, pooled regions for
GSE123496. All values are read from Table 5 and Supplementary Table S16b."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

O = "/Users/elmas/Desktop/MAKALE/06_MS_ANALIZI"
FIG = "/Users/elmas/Desktop/MAKALE/11_NEUROCHEMISTRY_INTERNATIONAL_FIGURE_REVISION/figures/redesign_work"
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.titlesize": 9, "axes.titleweight": "bold",
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "savefig.facecolor": "white", "figure.facecolor": "white",
})
UP, DN, GREY = "#D55E00", "#0072B2", "#999999"

# Every value is read from the table that reports it: Table 5 for the per-region
# comparisons (significant: Benjamini-Hochberg q < 0.05) and S16b for the donor-level
# normal-appearing white matter and myelin-adjusted tests (significant: p < 0.05).
PKG = "/Users/elmas/Desktop/MAKALE/10_NEUROCHEMISTRY_INTERNATIONAL"
T5 = pd.read_csv(f"{PKG}/tables/Table5_cross_disease_comparison.csv")
T5 = T5[T5.gene == "TRPC1"]
S16 = pd.read_csv(f"{PKG}/supplementary/S16b_multiple_sclerosis_donor_level.csv")
S16 = S16[S16.gene == "TRPC1"]


def t5(label, group, region, cohort=None):
    r = T5[(T5.group == group) & (T5.region == region)]
    if cohort:
        r = r[r.cohort == cohort]
    assert len(r) == 1, (group, region)
    return label, float(r.cliffs_delta.iloc[0]), "*" if r.q_value.iloc[0] < 0.05 else ""


def s16(label, comparison):
    r = S16[(S16.comparison == comparison) & (S16.unit == "donor")]
    assert len(r) == 1, comparison
    return label, float(r.cliffs_delta.iloc[0]), "†" if r.p.iloc[0] < 0.05 else ""


GAP = (None, None, None)
D = [t5(f"ALS · {reg}", "ALS", reg) for reg in
     ["Hippocampus", "Cerebellum", "Motor cortex (lateral)", "Frontal cortex",
      "Motor cortex (medial)", "Temporal cortex", "Occipital cortex",
      "Spinal cord (lumbar)", "Spinal cord (cervical)", "Spinal cord (thoracic)"]]
D += [GAP,
      t5("Other neurological · Cerebellum", "ONd", "Cerebellum"),
      t5("Other neurological · Frontal", "ONd", "Frontal cortex"),
      t5("Other neurological · Temporal", "ONd", "Temporal cortex"),
      GAP,
      t5("Alzheimer · fusiform gyrus", "AD", "Fusiform gyrus"),
      t5("Parkinson · BA9", "PD", "BA9"),
      GAP,
      t5("MS · donor level", "MS", "Donor level (10 MS vs 5 control donors)", "GSE138614"),
      t5("MS · lesions, sample level", "MS", "MS lesions vs control white matter", "GSE138614"),
      t5("MS · NAWM, sample level", "MS",
         "Normal-appearing white matter vs control white matter", "GSE138614"),
      s16("MS · NAWM, donor level", "NAWM vs control WM"),
      s16("MS · myelin-adjusted, donor level", "MS vs control, myelin+glia-adjusted TRPC1"),
      t5("MS · five regions pooled", "MS", "All five regions pooled (centred within region)",
         "GSE123496")]

fig, ax = plt.subplots(figsize=(7.1, 8.0), layout="constrained")
ys, labels = [], []
y = 0
for lab, d, symbol in D:
    if lab is None:
        y += 0.55
        continue
    c = UP if d > 0 else DN
    ax.barh(y, d, color=c, height=0.66, alpha=(1.0 if symbol else 0.42),
            edgecolor="none", zorder=2)
    ax.text(d + (0.035 if d > 0 else -0.035), y,
            f"{d:+.3f}".replace("-", "\u2212") + symbol,
            va="center", ha="left" if d > 0 else "right", fontsize=8,
            color="black", weight="bold" if symbol else "normal")
    ys.append(y); labels.append(lab)
    y += 1

ax.axvline(0, color="#333", lw=1.1, zorder=3)
ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=8.5)
ax.invert_yaxis()
ax.set_xlim(-1.08, 0.78)
ax.set_xlabel("Cliff's δ  (case − control)")
ax.set_title("TRPC1 expression across five disease settings",
             loc="left", fontsize=9, weight="bold")
ax.annotate("← decreased", xy=(-0.98, -0.75), fontsize=8.5, color="black", weight="bold")
ax.annotate("increased →", xy=(0.46, -0.75), fontsize=8.5, color="black", weight="bold")
# below the axis, where it cannot cross a bar or the zero line
ax.text(1.0, -0.075, "* Benjamini–Hochberg q < 0.05   ·   † uncorrected p < 0.05"
        "   ·   pale bars: not significant",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#666")


fig.savefig(os.path.join(FIG, "Figure8_TRPC1_disease_direction.png"))
fig.savefig(os.path.join(FIG, "Figure8_TRPC1_disease_direction.pdf"))
fig.savefig(os.path.join(FIG, "Figure8_TRPC1_disease_direction.svg"))
print("Figure 8 ->", FIG)

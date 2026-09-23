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
FIG = "/Users/elmas/Desktop/MAKALE/09_YAYIN_PAKETI/figures"
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight"})
UP, DN, GREY = "#b0322c", "#2c6fa8", "#9aa8a6"

# Every value is read from the table that reports it: Table 5 for the per-region
# comparisons (significant: Benjamini-Hochberg q < 0.05) and S16b for the donor-level
# normal-appearing white matter and myelin-adjusted tests (significant: p < 0.05).
PKG = "/Users/elmas/Desktop/MAKALE/09_YAYIN_PAKETI"
T5 = pd.read_csv(f"{PKG}/tables/Table5_cross_disease_comparison.csv")
T5 = T5[T5.gene == "TRPC1"]
S16 = pd.read_csv(f"{PKG}/supplementary/S16b_multiple_sclerosis_donor_level.csv")
S16 = S16[S16.gene == "TRPC1"]


def t5(label, group, region, cohort=None):
    r = T5[(T5.group == group) & (T5.region == region)]
    if cohort:
        r = r[r.cohort == cohort]
    assert len(r) == 1, (group, region)
    return label, float(r.cliffs_delta.iloc[0]), bool(r.q_value.iloc[0] < 0.05)


def s16(label, comparison):
    r = S16[(S16.comparison == comparison) & (S16.unit == "donor")]
    assert len(r) == 1, comparison
    return label, float(r.cliffs_delta.iloc[0]), bool(r.p.iloc[0] < 0.05)


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

fig, ax = plt.subplots(figsize=(8.8, 7.9))
ys, labels = [], []
y = 0
for lab, d, sig in D:
    if lab is None:
        y += 0.55
        continue
    c = UP if d > 0 else DN
    ax.barh(y, d, color=c, height=0.66, alpha=(1.0 if sig else 0.42),
            edgecolor="none", zorder=2)
    ax.text(d + (0.035 if d > 0 else -0.035), y,
            f"{d:+.3f}".replace("-", "\u2212") + ("*" if sig else ""),
            va="center", ha="left" if d > 0 else "right", fontsize=8,
            color=c, weight="bold" if sig else "normal")
    ys.append(y); labels.append(lab)
    y += 1

ax.axvline(0, color="#333", lw=1.1, zorder=3)
ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=8.5)
ax.invert_yaxis()
ax.set_xlim(-1.08, 0.78)
ax.set_xlabel("Cliff's δ  (case − control)")
ax.set_title("Direction of TRPC1 expression change across five diseases\n"
             "ALS is the only setting in which it increases significantly",
             loc="left", fontsize=11, weight="bold")
ax.annotate("← decreased", xy=(-0.98, -0.75), fontsize=8.5, color=DN, weight="bold")
ax.annotate("increased →", xy=(0.46, -0.75), fontsize=8.5, color=UP, weight="bold")
# below the axis, where it cannot cross a bar or the zero line
ax.text(1.0, -0.075, "* significant (Benjamini–Hochberg q < 0.05; donor-level tests: p < 0.05)"
        "   ·   pale bars: not significant",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#666")
ax.grid(axis="x", color="#e6ecec", lw=0.8, zorder=0)
ax.set_axisbelow(True)
fig.savefig(os.path.join(FIG, "Figure8_TRPC1_disease_direction.png"))
fig.savefig(os.path.join(FIG, "Figure8_TRPC1_disease_direction.pdf"))
print("Figure 8 ->", FIG)

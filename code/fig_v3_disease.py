#!/usr/bin/env python3
"""v3 — Figure 8: direction of TRPC1 change across five diseases.

ALS values: NYGC genome-wide per-region analysis (seven brain regions and three
spinal cord levels). MS values: donor level where the design allows it."""
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

# ALS/ONd/AD/PD degerleri mevcut raporlardan; MS bu calismadan
D = [
    ("ALS · Hippocampus",              +0.699, True,  "NYGC GSE153960"),
    ("ALS · Cerebellum",               +0.538, True,  "NYGC GSE153960"),
    ("ALS · Motor cortex (lateral)",   +0.511, True,  "NYGC GSE153960"),
    ("ALS · Frontal cortex",           +0.486, True,  "NYGC GSE153960"),
    ("ALS · Motor cortex (medial)",    +0.465, True,  "NYGC GSE153960"),
    ("ALS · Temporal cortex",          +0.447, True,  "NYGC GSE153960"),
    ("ALS · Occipital cortex",         +0.176, False, "NYGC GSE153960"),
    ("ALS · Spinal cord (lumbar)",     +0.112, False, "NYGC GSE153960"),
    ("ALS · Spinal cord (cervical)",   +0.058, False, "NYGC GSE153960"),
    ("ALS · Spinal cord (thoracic)",   -0.321, False, "NYGC GSE153960"),
    (None, None, None, None),
    ("Other neurological · Cerebellum", -0.407, True, "NYGC, same controls"),
    ("Other neurological · Frontal",    -0.717, True, "NYGC, same controls"),
    ("Other neurological · Temporal",   -0.571, True, "NYGC, same controls"),
    (None, None, None, None),
    ("Alzheimer · fusiform gyrus",  -0.447, True,  "GSE125583"),
    ("Parkinson · BA9",             -0.677, True,  "GSE68719"),
    (None, None, None, None),
    ("MS · donor level",                    -0.840, True,  "GSE138614 (this study)"),
    ("MS · lesions (sample level)",         -0.594, True,  "GSE138614 (this study)"),
    ("MS · NAWM, donor level",              -0.771, True,  "GSE138614 (this study)"),
    ("MS · myelin-adjusted, donor level",   -0.640, False, "GSE138614 (this study)"),
    ("MS · five regions pooled",            -0.226, False, "GSE123496 (this study)"),
]

fig, ax = plt.subplots(figsize=(8.8, 7.6))
ys, labels = [], []
y = 0
for lab, d, sig, src in D:
    if lab is None:
        y += 0.55
        continue
    c = UP if d > 0 else DN
    ax.barh(y, d, color=c, height=0.66, alpha=(1.0 if sig else 0.42),
            edgecolor="none", zorder=2)
    ax.text(d + (0.035 if d > 0 else -0.035), y, f"{d:+.3f}" + ("*" if sig else ""),
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
ax.text(0.995, 0.015, "* significant (BH-corrected q < 0.05, or p < 0.05 at donor level)"
        "   ·   pale bars: not significant",
        transform=ax.transAxes, ha="right", fontsize=7.5, color="#666")
ax.grid(axis="x", color="#e6ecec", lw=0.8, zorder=0)
ax.set_axisbelow(True)
fig.savefig(os.path.join(FIG, "Figure8_TRPC1_disease_direction.png"))
fig.savefig(os.path.join(FIG, "Figure8_TRPC1_disease_direction.pdf"))
print("Figure 8 ->", FIG)

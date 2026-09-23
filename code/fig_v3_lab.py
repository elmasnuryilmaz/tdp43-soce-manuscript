#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v3 — Figure 6: functional consequences of TDP-43 knockdown in SH-SY5Y.

Source data (all from the thesis laboratory records):
  A  TARDBP RT-qPCR            SHSY5Y_TDP43_qPCR_Ct_Data.xlsx, sheet TARDBP_Knockdown_Ct
  B  target-gene RT-qPCR       10_GRAPHPAD_PRISM_DOSYALARI/01_tez_sekil_kaynaklari/sekil_4.17_4.18_qPCR.pzfx
  C  WST-1 viability           .../sekil_4.21_wst1.pzfx
  D  representative Fura-2 traces, relabelled in English with the confirmed 10 uM CPA
     concentration (09_YAYIN_PAKETI/source_data/representative_Fura2_traces_relabelled.png)
  E  Fura-2 ER release / SOCE  .../sekil_4.20_fura2.pzfx
The same values are written to the supplementary source-data workbook by build_source_data.py.
"""
import os
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

TEZ = "/Users/elmas/Desktop/TEZ"
PZ = os.path.join(TEZ, "10_GRAPHPAD_PRISM_DOSYALARI", "01_tez_sekil_kaynaklari")
PKG = "/Users/elmas/Desktop/MAKALE/09_YAYIN_PAKETI"
FIG = f"{PKG}/figures"
TRACE = f"{PKG}/source_data/representative_Fura2_traces_relabelled.png"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titleweight": "bold",
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
})
KD_C, CT_C, MID = "#c0392b", "#2c6fa8", "#7f8c8d"


def pzfx(path):
    """Return {table title: {column title: [values]}}."""
    out = {}
    for tab in ET.parse(path).getroot().iter("Table"):
        cols = {}
        for col in tab.findall("YColumn"):
            cols[col.findtext("Title")] = [float(d.text) for d in col.iter("d")]
        out[tab.findtext("Title")] = cols
    return out


def sem(v):
    return np.std(v, ddof=1) / np.sqrt(len(v))


def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


# ---------------------------------------------------------------- source data
qp = pzfx(os.path.join(PZ, "sekil_4.17_4.18_qPCR.pzfx"))
fu = pzfx(os.path.join(PZ, "sekil_4.20_fura2.pzfx"))
ws = pzfx(os.path.join(PZ, "sekil_4.21_wst1.pzfx"))

td = pd.read_excel(os.path.join(TEZ, "SHSY5Y_TDP43_qPCR_Ct_Data.xlsx"),
                   sheet_name="TARDBP_Knockdown_Ct", skiprows=4, nrows=12,
                   names=["grup", "tekrar", "tarih", "TARDBP_Ct", "GAPDH_Ct",
                          "dCt", "ddCt", "rel"])
GRP = {"Kontrol": "Untransduced", "Hedef dışı kontrol": "Non-targeting shRNA",
       "shTDP-43": "shTDP-43"}
td["group"] = td.grup.map(GRP)

fig = plt.figure(figsize=(14.6, 7.6))
gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.35, 1.0], height_ratios=[1.0, 1.02],
                      wspace=0.30, hspace=0.42)
w = 0.36

# ------------------------------------------------------------------ A: TARDBP
ax = fig.add_subplot(gs[0, 0])
order = ["Untransduced", "Non-targeting shRNA", "shTDP-43"]
cols = [MID, CT_C, KD_C]
for i, (g, c) in enumerate(zip(order, cols)):
    v = td.loc[td.group == g, "rel"].astype(float).values
    ax.bar(i, v.mean(), 0.6, color=c, alpha=.85)
    ax.errorbar(i, v.mean(), yerr=sem(v), color="#222", capsize=4, lw=1.2)
    ax.scatter(np.full(len(v), i) + np.linspace(-.13, .13, len(v)), v, s=22,
               color="#222", zorder=3)
kd = td.loc[td.group == "shTDP-43", "rel"].astype(float).values
un = td.loc[td.group == "Untransduced", "rel"].astype(float).values
nt = td.loc[td.group == "Non-targeting shRNA", "rel"].astype(float).values
F, p_anova = stats.f_oneway(np.log2(un), np.log2(nt), np.log2(kd))
ax.set_xticks(range(3))
ax.set_xticklabels(["Untransduced", "Non-targeting\nshRNA", "shTDP-43"], fontsize=8)
ax.set_ylabel("TARDBP mRNA (2$^{-\\Delta\\Delta Ct}$)")
ax.set_ylim(0, 1.32)
ax.plot([1, 1, 2, 2], [1.14, 1.19, 1.19, 1.14], lw=1, color="#222")
ax.text(1.5, 1.20, "****", ha="center", fontsize=10)
ax.text(2, 0.45, "−94.4% vs\nnon-targeting\n−94.8% vs\nuntransduced",
        ha="center", va="center", fontsize=7.2, color="#444")
ax.set_title("A · TARDBP knockdown efficiency", loc="left")

# ------------------------------------------------------------ B: target mRNAs
ax = fig.add_subplot(gs[0, 1])
genes = ["TRPC1", "STIM1", "ORAI1", "ATP2A3"]
x = np.arange(len(genes))
for i, g in enumerate(genes):
    c = np.array(qp[g]["Control"]); k = np.array(qp[g]["shTDP-43"])
    ax.bar(i - w/2, c.mean(), w, color=CT_C, alpha=.85)
    ax.bar(i + w/2, k.mean(), w, color=KD_C, alpha=.85)
    ax.errorbar(i - w/2, c.mean(), yerr=sem(c), color="#222", capsize=3, lw=1.1)
    ax.errorbar(i + w/2, k.mean(), yerr=sem(k), color="#222", capsize=3, lw=1.1)
    ax.scatter(np.full(4, i - w/2) + np.linspace(-.09, .09, 4), c, s=16, color="#222", zorder=3)
    ax.scatter(np.full(4, i + w/2) + np.linspace(-.09, .09, 4), k, s=16, color="#222", zorder=3)
    p = stats.ttest_ind(c, k).pvalue
    top = max(k.mean() + sem(k), c.mean() + sem(c))
    ax.text(i, top + 0.22, stars(p), ha="center", fontsize=10)
    ax.text(i, top + 0.62, f"{k.mean()/c.mean():.1f}×", ha="center", fontsize=7.6, color="#444")
ax.set_xticks(x); ax.set_xticklabels(genes, fontsize=9)
ax.set_ylabel("relative mRNA (2$^{-\\Delta\\Delta Ct}$)")
ax.set_ylim(0, 5.0)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=CT_C, alpha=.85),
                   plt.Rectangle((0, 0), 1, 1, color=KD_C, alpha=.85)],
          labels=["Non-targeting shRNA", "shTDP-43"], frameon=False, fontsize=8, loc="upper left")
ax.set_title("B · SOCE-associated target mRNAs", loc="left")

# ------------------------------------------------------------------ C: WST-1
ax = fig.add_subplot(gs[0, 2])
for i, (lab, tab) in enumerate([("24 h", ws["WST_1_24h"]), ("48 h", ws["WST_1_48h"])]):
    c = np.array(tab["Control"]); k = np.array(tab["shTDP-43"])
    ax.bar(i - w/2, c.mean(), w, color=CT_C, alpha=.85)
    ax.bar(i + w/2, k.mean(), w, color=KD_C, alpha=.85)
    ax.errorbar(i - w/2, c.mean(), yerr=sem(c), color="#222", capsize=3, lw=1.1)
    ax.errorbar(i + w/2, k.mean(), yerr=sem(k), color="#222", capsize=3, lw=1.1)
    ax.scatter(np.full(4, i - w/2) + np.linspace(-.09, .09, 4), c, s=16, color="#222", zorder=3)
    ax.scatter(np.full(4, i + w/2) + np.linspace(-.09, .09, 4), k, s=16, color="#222", zorder=3)
    p = stats.ttest_ind(c, k).pvalue
    ax.text(i, max(c.mean(), k.mean()) + 7, stars(p), ha="center", fontsize=10)
ax.axhline(100, color="#999", lw=.8, ls="--")
ax.set_xticks([0, 1]); ax.set_xticklabels(["24 h", "48 h"], fontsize=9)
ax.set_ylabel("viability (% of non-targeting control)")
ax.set_ylim(0, 145)
ax.set_title("C · WST-1 viability (n = 4 wells)", loc="left")

# ------------------------------------------------- D: representative traces
ax = fig.add_subplot(gs[1, 0:2])
ax.imshow(mpimg.imread(TRACE))
ax.set_axis_off()
ax.set_title("D · Representative Fura-2 traces", loc="left")

# --------------------------------------------------------- E: Fura-2 group data
ax = fig.add_subplot(gs[1, 2])
pairs = [("ER Ca²⁺ release", fu["ER_Ca2_release"]), ("SOCE", fu["SOCE"])]
for i, (lab, tab) in enumerate(pairs):
    c = np.array(tab["Control"]); k = np.array(tab["TDP-43 KD"])
    ax.bar(i - w/2, c.mean(), w, color=CT_C, alpha=.85)
    ax.bar(i + w/2, k.mean(), w, color=KD_C, alpha=.85)
    ax.errorbar(i - w/2, c.mean(), yerr=sem(c), color="#222", capsize=3, lw=1.1)
    ax.errorbar(i + w/2, k.mean(), yerr=sem(k), color="#222", capsize=3, lw=1.1)
    ax.scatter(np.full(3, i - w/2) + np.linspace(-.08, .08, 3), c, s=20, color="#222", zorder=3)
    ax.scatter(np.full(3, i + w/2) + np.linspace(-.08, .08, 3), k, s=20, color="#222", zorder=3)
    p = stats.ttest_ind(c, k).pvalue
    ax.text(i, max(c.max(), k.max()) + 0.13,
            f"p = {p:.3f}" if p >= 0.05 else f"p = {p:.4f}", ha="center", fontsize=8)
ax.set_xticks([0, 1])
ax.set_xticklabels(["ER Ca²⁺ release\n(CPA, Ca²⁺-free)",
                    "SOCE\n(1.5 mM Ca²⁺)"], fontsize=8)
ax.set_ylabel("Δ (F340/F380)")
ax.set_ylim(0, 2.25)
ax.set_title("E · Cytosolic Ca²⁺ responses (n = 3)", loc="left")

fig.suptitle("Functional consequences of TDP-43 knockdown in SH-SY5Y cells",
             y=0.98, fontsize=11, weight="bold")
fig.savefig(os.path.join(FIG, "Figure6_functional_consequences.png"))
fig.savefig(os.path.join(FIG, "Figure6_functional_consequences.pdf"))
plt.close(fig)

# ---------------------------------------------------------------- stat report
print("Figure 6 written ->", FIG)
print(f"A  one-way ANOVA on log2 relative expression: F = {F:.1f}, p = {p_anova:.3g}")
print(f"   TARDBP: untransduced {un.mean():.3f}, non-targeting {nt.mean():.3f}, "
      f"shTDP-43 {kd.mean():.4f}  ->  -{100*(1-kd.mean()/nt.mean()):.1f}% / "
      f"-{100*(1-kd.mean()/un.mean()):.1f}%")
for g in genes:
    c = np.array(qp[g]["Control"]); k = np.array(qp[g]["shTDP-43"])
    print(f"B  {g:7s} {k.mean()/c.mean():.2f}x  p = {stats.ttest_ind(c,k).pvalue:.3g}  "
          f"(mean+-SEM {k.mean():.2f}+-{sem(k):.2f})")
for lab, tab in [("24 h", ws["WST_1_24h"]), ("48 h", ws["WST_1_48h"])]:
    c = np.array(tab["Control"]); k = np.array(tab["shTDP-43"])
    print(f"C  WST-1 {lab}: control {c.mean():.1f}+-{sem(c):.1f}  KD {k.mean():.1f}+-{sem(k):.1f}  "
          f"p = {stats.ttest_ind(c,k).pvalue:.3g}")
for lab, tab in pairs:
    c = np.array(tab["Control"]); k = np.array(tab["TDP-43 KD"])
    print(f"E  {lab:18s} control {c.mean():.3f}+-{sem(c):.3f}  KD {k.mean():.3f}+-{sem(k):.3f}  "
          f"p = {stats.ttest_ind(c,k).pvalue:.4f}")

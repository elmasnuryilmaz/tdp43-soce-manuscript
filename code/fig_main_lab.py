#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 1: functional consequences of TDP-43 knockdown in SH-SY5Y.

Source data (all from the thesis laboratory records):
  A  TARDBP RT-qPCR            SHSY5Y_TDP43_qPCR_Ct_Data.xlsx, sheet TARDBP_Knockdown_Ct
  B  target-gene RT-qPCR       10_GRAPHPAD_PRISM_DOSYALARI/01_tez_sekil_kaynaklari/sekil_4.17_4.18_qPCR.pzfx
  C  WST-1 metabolic activity  .../sekil_4.21_wst1.pzfx
  D/E Fura-2 ER release / SOCE .../sekil_4.20_fura2.pzfx
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

TEZ = "/Users/elmas/Desktop/TEZ"
PZ = os.path.join(TEZ, "10_GRAPHPAD_PRISM_DOSYALARI", "01_tez_sekil_kaynaklari")
PKG = "/Users/elmas/Desktop/MAKALE/11_NEUROCHEMISTRY_INTERNATIONAL_FIGURE_REVISION"
FIG = f"{PKG}/figures/main"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.titlesize": 9, "axes.titleweight": "bold",
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "savefig.facecolor": "white", "figure.facecolor": "white",
})
KD_C, CT_C, MID = "#D55E00", "#0072B2", "#999999"


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

fig = plt.figure(figsize=(7.1, 8.7), layout="constrained")
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.9, 0.95])
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

ax.set_title("A  TARDBP mRNA", loc="left")

# ------------------------------------------------------------ B: target mRNAs
ax = fig.add_subplot(gs[0, 1])
genes = ["TRPC1", "STIM1", "ORAI1", "ATP2A3"]
x = np.arange(len(genes))
raw_p = {g: stats.ttest_ind(np.array(qp[g]["Control"]), np.array(qp[g]["shTDP-43"])).pvalue for g in genes}
_order = sorted(raw_p, key=raw_p.get)
holm = {}
running = 0.0
for rank, g in enumerate(_order):
    running = max(running, raw_p[g] * (len(genes) - rank))
    holm[g] = min(1.0, running)
for i, g in enumerate(genes):
    c = np.array(qp[g]["Control"]); k = np.array(qp[g]["shTDP-43"])
    ax.bar(i - w/2, c.mean(), w, color=CT_C, alpha=.85)
    ax.bar(i + w/2, k.mean(), w, color=KD_C, alpha=.85)
    ax.errorbar(i - w/2, c.mean(), yerr=sem(c), color="#222", capsize=3, lw=1.1)
    ax.errorbar(i + w/2, k.mean(), yerr=sem(k), color="#222", capsize=3, lw=1.1)
    ax.scatter(np.full(4, i - w/2) + np.linspace(-.09, .09, 4), c, s=16, color="#222", zorder=3)
    ax.scatter(np.full(4, i + w/2) + np.linspace(-.09, .09, 4), k, s=16, color="#222", zorder=3)
    p = holm[g]
    top = max(k.mean() + sem(k), c.mean() + sem(c))
    ax.text(i, top + 0.22, stars(p), ha="center", fontsize=10)
    ax.text(i, top + 0.62, f"{k.mean()/c.mean():.1f}×", ha="center", fontsize=7.6, color="#444")
ax.set_xticks(x); ax.set_xticklabels(genes, fontsize=9)
ax.set_ylabel("relative mRNA (2$^{-\\Delta\\Delta Ct}$)")
ax.set_ylim(0, 5.0)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=CT_C, alpha=.85),
                   plt.Rectangle((0, 0), 1, 1, color=KD_C, alpha=.85)],
          labels=["Non-targeting shRNA", "shTDP-43"], frameon=False, fontsize=8, loc="upper left")
ax.set_title("B  Ca$^{2+}$-related transcripts", loc="left")

# ------------------------------------------------------------------ C: WST-1
ax = fig.add_subplot(gs[1, 0])
tab = ws["WST_1_48h"]
c = np.array(tab["Control"]); k = np.array(tab["shTDP-43"])
ax.bar(-w/2, c.mean(), w, color=CT_C, alpha=.85)
ax.bar(w/2, k.mean(), w, color=KD_C, alpha=.85)
ax.errorbar(-w/2, c.mean(), yerr=sem(c), color="#222", capsize=3, lw=1.1)
ax.errorbar(w/2, k.mean(), yerr=sem(k), color="#222", capsize=3, lw=1.1)
ax.scatter(np.full(4, -w/2) + np.linspace(-.09, .09, 4), c, s=16, color="#222", zorder=3)
ax.scatter(np.full(4, w/2) + np.linspace(-.09, .09, 4), k, s=16, color="#222", zorder=3)
ax.axhline(100, color="#999", lw=.8, ls="--")
ax.set_xticks([-w/2, w/2]); ax.set_xticklabels(["Non-targeting\nshRNA", "shTDP-43"], fontsize=8)
ax.set_ylabel("WST-1 signal (% of control)")
ax.set_ylim(0, 125)
ax.text(0, 112, "−38.5%", ha="center", fontsize=8, color="#444")
ax.set_title("C  WST-1 at 48 h", loc="left")

# --------------------------------------------------------- D/E: Fura-2 group data
pairs = [("ER Ca$^{2+}$ release", fu["ER_Ca2_release"]), ("SOCE", fu["SOCE"])]
for label, tab, slot, panel, ylim in [
    ("ER Ca$^{2+}$ release", fu["ER_Ca2_release"], gs[1, 1], "D", 0.52),
    ("SOCE after Ca$^{2+}$ readdition", fu["SOCE"], gs[2, :], "E", 2.25),
]:
    ax = fig.add_subplot(slot)
    c = np.array(tab["Control"]); k = np.array(tab["TDP-43 KD"])
    for x0, v, color in [(0, c, CT_C), (1, k, KD_C)]:
        ax.bar(x0, v.mean(), width=0.55, color=color, alpha=0.9)
        ax.errorbar(x0, v.mean(), yerr=sem(v), color="#222222", capsize=4, lw=1.1)
        ax.scatter(x0 + np.linspace(-0.12, 0.12, len(v)), v,
                   s=30, color="#222222", zorder=3)
    pval = stats.ttest_ind(c, k).pvalue
    ax.set_xticks([0, 1], ["Non-targeting shRNA", "shTDP-43"])
    ax.set_xlim(-0.5, 1.5); ax.set_ylim(0, ylim)
    ax.set_ylabel("Δ(F340/F380)")
    ax.set_title(f"{panel}  {label}", loc="left")
    ax.text(0.5, 0.96, f"p = {pval:.4f}" if pval < 0.05 else f"p = {pval:.3f}",
            transform=ax.transAxes, ha="center", va="top", fontsize=8, color="black")

fig.savefig(os.path.join(FIG, "Figure1_functional_consequences.png"))
fig.savefig(os.path.join(FIG, "Figure1_functional_consequences.pdf"))
plt.close(fig)

# ---------------------------------------------------------------- stat report
print("Figure 6 written ->", FIG)
print(f"A  one-way ANOVA on log2 relative expression: F = {F:.1f}, p = {p_anova:.3g}")
print(f"   TARDBP: untransduced {un.mean():.3f}, non-targeting {nt.mean():.3f}, "
      f"shTDP-43 {kd.mean():.4f}  ->  -{100*(1-kd.mean()/nt.mean()):.1f}% / "
      f"-{100*(1-kd.mean()/un.mean()):.1f}%")
for g in genes:
    c = np.array(qp[g]["Control"]); k = np.array(qp[g]["shTDP-43"])
    print(f"B  {g:7s} {k.mean()/c.mean():.2f}x  raw p = {stats.ttest_ind(c,k).pvalue:.3g}, Holm p = {holm[g]:.3g}  "
          f"(mean+-SEM {k.mean():.2f}+-{sem(k):.2f})")
tab = ws["WST_1_48h"]
c = np.array(tab["Control"]); k = np.array(tab["shTDP-43"])
print(f"C  WST-1 48 h: control {c.mean():.1f}+-{sem(c):.1f}  KD {k.mean():.1f}+-{sem(k):.1f}  (descriptive; wells from one experiment)")
for lab, tab in pairs:
    c = np.array(tab["Control"]); k = np.array(tab["TDP-43 KD"])
    print(f"D  {lab:18s} control {c.mean():.3f}+-{sem(c):.3f}  KD {k.mean():.3f}+-{sem(k):.3f}  "
          f"p = {stats.ttest_ind(c,k).pvalue:.4f}")

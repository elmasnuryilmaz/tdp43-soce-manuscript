#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v3 — junction-level figures for the submission package.

Figure 4 (cryptic discovery + NMD), Figure 5 (specificity + corrected APA),
Figure 9 (NYGC junction-level cryptic STMN2). Numbering follows MANUSCRIPT_v3."""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

D = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI"
OUT = f"{D}/sonuclar"; FIG = "/Users/elmas/Desktop/MAKALE/09_YAYIN_PAKETI/figures"; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titleweight": "bold",
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight"})
KD_C, CT_C, ACC, WARN = "#c0392b", "#2c6fa8", "#0d6259", "#b07d0a"
PANEL = "/Users/elmas/Desktop/TEZ/output/calcium_panel_overlap_analysis_2026-04-25/01_gene_panels"

def kaydet(fig, ad):
    fig.savefig(f"{FIG}/{ad}.png"); fig.savefig(f"{FIG}/{ad}.pdf"); plt.close(fig)
    print("written:", ad)

# ============================================ FIGURE 7: kriptik keşif + NMD
fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.0),
                         gridspec_kw={"width_ratios": [1.1, 1, 1.15], "wspace": 0.34})
kr = pd.read_csv(f"{OUT}/YUKSEK_GUVEN_OWN_SH_SY5Y.tsv", sep="\t")
POS = ["STMN2", "ACTL6B", "UNC13A", "PFKP", "HDGFL2", "AGRN", "ATG4B",
       "ELAVL3", "SETD5", "ARHGAP32", "GPSM2", "RSF1"]

ax = axes[0]
d = kr[kr["gene"].isin(POS)].sort_values("dPSI").drop_duplicates("gene", keep="last")
ax.barh(d["gene"], d["dPSI"], color=ACC)
for y, (v, ok, oc) in enumerate(zip(d["dPSI"], d["okuma_KD"], d["okuma_CTRL"])):
    ax.text(v + .012, y, f"{int(ok)}/{int(oc)}", va="center", fontsize=7, color="#555")
ax.set_xlim(0, 1.32); ax.set_xlabel("ΔPSI (knockdown − control)")
ax.set_title("A · Cryptic positive controls, recovered de novo\n(read counts KD/control)",
             loc="left")

ax = axes[1]
adlar = ["Tier1_SOCE_TRP_51", "Tier2_Channel_Release_Transport_117",
         "Tier3_Curated_Calcium_Handling_258", "Tier4_Expanded_Calcium_Associated_732"]
kisa = ["Tier 1\n(51)", "Tier 2\n(117)", "Tier 3\n(258)", "Tier 4\n(732)"]
n_gen = []
for a in adlar:
    gl = set(pd.read_csv(f"{PANEL}/{a}.csv")["gene_upper"])
    n_gen.append(kr[kr["gene"].isin(gl)]["gene"].nunique())
ax.bar(kisa, n_gen, color=[KD_C if v == 0 else CT_C for v in n_gen])
for i, v in enumerate(n_gen):
    ax.text(i, v + .1, str(v), ha="center", fontsize=10, fontweight="bold")
ax.set_ylabel("genes carrying a cryptic event")
ax.set_ylim(0, max(n_gen) + 1.4)
ax.set_title("B · Cryptic events in the Ca²⁺ panels\n(none in the core SOCE/TRP set)",
             loc="left")

ax = axes[2]
# Use the dependence-aware NMD sensitivity analysis: four intervention
# conditions are the inferential units, rather than eight contrasts that
# reuse the same control libraries.
nmd = pd.read_csv(f"{OUT}/NMD_etkilesim_paylasimli_kontrol_t4.tsv", sep="\t")
old_nmd = pd.read_csv(f"{OUT}/NMD_etkilesim_tum_genler.tsv", sep="\t", usecols=["ens", "sembol"])
nmd = nmd.merge(old_nmd, on="ens", how="left")
nmd["symbol"] = nmd["sembol"].fillna(nmd["ens"])
krg = {g for g in pd.read_csv(f"{OUT}/YUKSEK_GUVEN_OWN_SH_SY5Y.tsv", sep="\t")["gene"].astype(str)
       if not g.startswith("ENSG") and g != "."}
arka = nmd[~nmd["symbol"].isin(krg)]["interaction_log2"]
on = nmd[nmd["symbol"].isin(krg)]["interaction_log2"]
bins = np.linspace(-2, 2, 60)
ax.hist(arka, bins=bins, density=True, color="#ccc", label=f"other genes (n={len(arka)})")
ax.hist(on, bins=bins, density=True, histtype="step", lw=2, color=KD_C,
        label=f"cryptic-junction genes (n={len(on)})")
cb = nmd[nmd["symbol"] == "CBARP"]["interaction_log2"].iloc[0]
ax.axvline(cb, color=ACC, lw=2)
ax.annotate("CBARP\n+1.52", (cb - 0.06, ax.get_ylim()[1] * 0.70), fontsize=8, color=ACC,
            ha="right", fontweight="bold")
ax.set_xlabel("TDP-43-specific NMD interaction (log₂)")
ax.set_ylabel("density"); ax.legend(frameon=False, fontsize=7.5)
ax.set_title("C · Exploratory NMD interaction\n(four condition-level units)", loc="left")
kaydet(fig, "Figure4_cryptic_discovery_and_NMD")

# ============================== FIGURE 8: NYGC birleşim düzeyi kriptik STMN2
S = pd.read_csv(f"{OUT}/NYGC_kriptik_PSI_ornek_duzeyi.tsv", sep="\t", low_memory=False)
S = S[S["STMN2_toplam"] >= 20]
S["kontrol"] = S["grup"].astype(str).str.contains("Non-Neurological")
S["als"] = S["grup"].astype(str).str.strip() == "ALS Spectrum MND"
SIRA = ["Spinal Cord Lumbar", "Spinal Cord Cervical", "Spinal Cord Thoracic",
        "Cortex Motor Unspecified", "Cortex Motor Lateral", "Cortex Motor Medial",
        "Cortex Temporal", "Hippocampus", "Cortex Frontal", "Cortex Occipital", "Cerebellum"]
KISA = {"Spinal Cord Lumbar": "Spinal cord, lumbar",
        "Spinal Cord Cervical": "Spinal cord, cervical",
        "Spinal Cord Thoracic": "Spinal cord, thoracic",
        "Cortex Motor Unspecified": "Motor cortex (unspec.)",
        "Cortex Motor Lateral": "Motor cortex, lateral",
        "Cortex Motor Medial": "Motor cortex, medial",
        "Cortex Temporal": "Temporal cortex", "Hippocampus": "Hippocampus",
        "Cortex Frontal": "Frontal cortex", "Cortex Occipital": "Occipital cortex",
        "Cerebellum": "Cerebellum"}

fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.4),
                         gridspec_kw={"width_ratios": [1.7, 1], "wspace": 0.24})
ax = axes[0]
oran_als, oran_ctl, etiket = [], [], []
for d_ in SIRA:
    a = S[S["als"] & (S["doku"] == d_)]; c = S[S["kontrol"] & (S["doku"] == d_)]
    if len(a) < 20: continue
    oran_als.append(100 * (a["STMN2_kriptik"] > 0).mean())
    oran_ctl.append(100 * (c["STMN2_kriptik"] > 0).mean() if len(c) else 0.0)
    etiket.append(f"{KISA[d_]} (n={len(a)})")
x = np.arange(len(etiket)); w = 0.38
ax.bar(x - w/2, oran_als, w, color=KD_C, label="ALS")
ax.bar(x + w/2, oran_ctl, w, color=CT_C, label="Non-neurological control")
ax.set_xticks(x); ax.set_xticklabels(etiket, fontsize=7.5, rotation=38, ha="right")
ax.set_ylabel("samples carrying the cryptic STMN2 junction (%)")
ax.set_title("A · Cryptic STMN2 by region, NYGC cohort", loc="left")
ax.legend(frameon=False, fontsize=8)

ax = axes[1]
R = pd.read_csv(f"{OUT}/NYGC_kriptikPSI_korelasyon.tsv", sep="\t")
p = R[R["hedef"].isin(["SNAP25", "TRPC1", "SARAF", "ATP2A2", "CBARP"])].pivot_table(
    index=["hedef", "doku"], columns="vekil", values="rho").dropna()
renk = {"SNAP25": "#444", "TRPC1": KD_C, "SARAF": ACC, "ATP2A2": WARN, "CBARP": "#7d3c98"}
for h in p.index.get_level_values(0).unique():
    s_ = p.loc[h]
    ax.scatter(s_["gen_duzeyi_STMN2"], s_["kriptik_STMN2_PSI"], s=42,
               color=renk[h], label=h, alpha=.85, edgecolor="w", linewidth=.6)
ax.axhline(0, color="#999", lw=.8); ax.axvline(0, color="#999", lw=.8)
ax.plot([-1, 1], [-1, 1], "--", color="#bbb", lw=.8)
ax.set_xlim(-0.8, 1.0); ax.set_ylim(-0.45, 0.45)
ax.set_xlabel("Spearman ρ — proxy: gene-level STMN2")
ax.set_ylabel("Spearman ρ — proxy: cryptic STMN2 PSI")
ax.set_title("B · Two proxies for TDP-43 loss,\nsame ALS samples", loc="left")
ax.legend(frameon=False, fontsize=7, loc="upper left")
kaydet(fig, "Figure9_NYGC_cryptic_STMN2")

# ================================ FIGURE 9: RBP özgüllüğü + APA
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3),
                         gridspec_kw={"width_ratios": [1, 1.2, 1.0], "wspace": 0.42})

ax = axes[0]
POSG = {"STMN2","UNC13A","HDGFL2","ACTL6B","AGRN","KALRN","ARHGAP32","PFKP","ATG4B",
        "SETD5","ELAVL3","POLDIP3","CAMK2B","RSF1","GPSM2","SYNJ2"}
SIRA2 = ["SH_SY5Y","SH_SY5Y_DOZ25","iPSC_koloni","K562_mRNA","NSC34","iPSC_MN",
         "K562_totalRNA","C2C12","Fare_striatum","iPSC_MN_FUS","iPSC_MN_TAF15"]
ETI = {"SH_SY5Y":"SH-SY5Y 75 ng/mL","SH_SY5Y_DOZ25":"SH-SY5Y 25 ng/mL",
       "iPSC_koloni":"iPSC colonies","K562_mRNA":"K562 poly(A)+","NSC34":"NSC34",
       "iPSC_MN":"iPSC-MN · TDP-43","K562_totalRNA":"K562 total RNA","C2C12":"C2C12",
       "Fare_striatum":"Mouse striatum","iPSC_MN_FUS":"iPSC-MN · FUS",
       "iPSC_MN_TAF15":"iPSC-MN · TAF15"}
vals, cols, labs = [], [], []
for ds in SIRA2:
    f = f"{OUT}/RT_KRIPTIK_{ds}.tsv"
    if not os.path.exists(f): continue
    gu = set(pd.read_csv(f, sep="\t")["gene"].astype(str).str.upper())
    vals.append(len(gu & POSG)); labs.append(ETI[ds])
    cols.append("#8e44ad" if ds.endswith(("FUS", "TAF15"))
                else (KD_C if ds.startswith("SH_SY5Y") else CT_C))
y = np.arange(len(vals))[::-1]
ax.barh(y, vals, color=cols)
for yy, v in zip(y, vals):
    ax.text(v + 0.25, yy, str(v), va="center", fontsize=8, color="#444")
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=7.5)
ax.set_xlabel("positive-control genes recovered (of 16)")
ax.set_xlim(0, 18.5)
ax.legend(handles=[Patch(color=KD_C, label="TDP-43 KD (primary model)"),
                   Patch(color=CT_C, label="TDP-43 KD (other)"),
                   Patch(color="#8e44ad", label="FUS / TAF15 KD")],
          frameon=False, fontsize=7, loc="lower right")
ax.set_title("A · Positive-control recovery by comparison\n(permissive definition)", loc="left")

ax = axes[1]
A = pd.read_csv(f"{OUT}/APA_corrected_full_core_summary.tsv", sep="\t")
sig = A[A["delta"].abs() >= 0.05].sort_values("delta")
renk2 = np.where(sig["measure"] == "IPA_index", ACC, "#8e44ad")
ax.barh(np.arange(len(sig)), sig["delta"], color=renk2)
ax.set_yticks(np.arange(len(sig)))
ax.set_yticklabels([f"{g} · {b.replace('termexon','terminal exon')}"
                    for g, b in zip(sig["gene"], sig["unit"])], fontsize=6.5)
ax.axvline(0, color="#666", lw=.8)
ax.set_xlabel("Δ coverage index (knockdown − control)")
ax.legend(handles=[Patch(color=ACC, label="intronic polyadenylation index"),
                   Patch(color="#8e44ad", label="distal 3′UTR usage index")],
          frameon=True, framealpha=0.92, edgecolor="none", fontsize=7, loc="lower right")
ax.set_title("B · Alternative polyadenylation in Ca²⁺ genes,\nSH-SY5Y", loc="left")

ax = axes[2]
MANx = pd.read_csv(f"{D}/kod/ornekler.tsv", sep="\t")
sub = MANx[MANx["dataset"] == "SH_SY5Y"]
kd = list(sub[sub.group == "KD"]["sample"]); ct = list(sub[sub.group == "CTRL"]["sample"])
b = pd.read_csv(f"{OUT}/apa_bedcov_SH_SY5Y_corrected_full.tsv", sep="\t", header=None,
                names=["chrom", "start", "end", "isim"] + kd + ct)
st = b[b["isim"].str.startswith("STMN2|intron2|")]
if len(st) == 2:
    i5 = st[st["isim"].str.contains("I5")].iloc[0]
    i3 = st[st["isim"].str.contains("I3")].iloc[0]
    d5k = i5[kd].values / (i5["end"] - i5["start"]); d3k = i3[kd].values / (i3["end"] - i3["start"])
    d5c = i5[ct].values / (i5["end"] - i5["start"]); d3c = i3[ct].values / (i3["end"] - i3["start"])
    xs = np.array([0, 1])
    for j in range(len(kd)):
        # v3 fix: a measured coverage of zero is data, not missing data.
        # It is plotted at zero and marked with an open symbol.
        ax.plot(xs, [d5c[j], d3c[j]], "-o", color=CT_C, alpha=.8, ms=5,
                label="Control" if j == 0 else None)
        ax.plot(xs + 0.08, [d5k[j], d3k[j]], "-o", color=KD_C, alpha=.8, ms=5,
                label="TDP-43 KD" if j == 0 else None)
        if d3k[j] == 0:
            ax.plot([xs[1] + 0.08], [0.0], "o", ms=10, mfc="white", mec=KD_C,
                    mew=1.8, zorder=5)
            ax.annotate("measured zero", xy=(xs[1] + 0.08, 0.0),
                        xytext=(-6, 14), textcoords="offset points",
                        fontsize=7, color=KD_C, ha="right")
    ax.set_xticks([0.04, 1.04]); ax.set_xticklabels(["intron 5′ end", "intron 3′ end"])
    ax.set_ylabel("mean coverage depth")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("C · Corrected coverage across STMN2 intron 2\n(one knockdown 3′ window has a measured coverage of zero)",
                 loc="left")
kaydet(fig, "Figure5_specificity_and_APA")

#!/usr/bin/env python3
"""A5 — NMD inhibisyonu etkileşim testi (GSE307054, PRJNA1235234).

Tasarım: i3Neuron, iki faktör — TDP-43 (WT/KD) ve NMD (WT/inhibe).
Yığın (batch) TDP durumuyla karışık olduğu için etkileşim, her yığının
KENDİ kontrolüne göre hesaplanan farkların farkı olarak kurulur:

  d_KD(c) = log2CPM(T{c}) - log2CPM(TDPKD)     [yığın 2 içinde]
  d_WT(c) = log2CPM(C{c}) - log2CPM(CON)       [yığın 1 içinde]
  etkilesim(c) = d_KD(c) - d_WT(c)

c ∈ {X=XRN1, XS=XRN1+SMG6, XU=XRN1+UPF1, US=UPF1+SMG6}. Yığın etkisi
her farkın içinde sadeleşir. Pozitif etkileşim = transkript TDP-43 kaybında
üretilip NMD ile yıkılıyor (kriptik PTC imzası).
"""
import os, sys, gzip
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_junc as L

BASE = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI"
OUT = f"{BASE}/sonuclar"; os.makedirs(OUT, exist_ok=True)

cnt = pd.read_csv(f"{BASE}/nmd_GSE307054/GSE307054_counts.csv.gz", index_col=0)
meta = pd.read_csv(f"{BASE}/nmd_GSE307054/GSE307054_meta.csv.gz", index_col=0)
print("sayım matrisi:", cnt.shape)

# --- normalizasyon: yazarların DESeq2 boyut faktörleri
sf = meta.set_index("sample")["sizeFactor"]
norm = cnt.div(sf[cnt.columns].values, axis=1)
lg = np.log2(norm + 1)

# --- gen sembolü
g2n = pd.read_csv("/Volumes/10TBElmas/thesis_addendum_2026/ref/gene_id2name_v47.tsv",
                  sep="\t", header=None, names=["gid", "sym"])
g2n["ens"] = g2n["gid"].str.split(".").str[0]
sym = g2n.drop_duplicates("ens").set_index("ens")["sym"]

PAIRS = [("X", "TX", "CX"), ("XS", "TXS", "CXS"), ("XU", "TXU", "CXU"), ("US", "TUS", "CUS")]
REPS = ["1", "2"]

inter = {}
for c, tpref, cpref in PAIRS:
    for r in REPS:
        t, tc = f"{tpref}_{r}", f"TDPKD_{r}"
        w, wc = f"{cpref}_{r}", f"CON_{r}"
        inter[f"{c}_{r}"] = (lg[t] - lg[tc]) - (lg[w] - lg[wc])
I = pd.DataFrame(inter)

# yeterli ifade filtresi: NMDi altında TDP-KD örneklerinde ortalama sayım
tdp_nmdi = [f"{p[1]}_{r}" for p in PAIRS for r in REPS]
expressed = (norm[tdp_nmdi].mean(axis=1) >= 10) & (norm[["TDPKD_1", "TDPKD_2"]].mean(axis=1) >= 5)
I = I[expressed]
print("test edilen gen:", len(I))

res = pd.DataFrame({
    "ens": I.index,
    "etkilesim_log2": I.mean(axis=1).values,
    "sd": I.std(axis=1).values,
    "n_pozitif": (I > 0).sum(axis=1).values,
})
tt = stats.ttest_1samp(I.values, 0.0, axis=1)
res["p"] = tt.pvalue
res["q"] = L.bh(res["p"].values)
res["sembol"] = res["ens"].map(sym).fillna(res["ens"])
res = res.sort_values("etkilesim_log2", ascending=False).reset_index(drop=True)
res["sira"] = np.arange(1, len(res) + 1)
res["yuzdelik"] = 100 * (1 - (res["sira"] - 1) / len(res))
res.to_csv(f"{OUT}/NMD_etkilesim_tum_genler.tsv", sep="\t", index=False)

# --- pozitif kontroller ve SOCE paneli
POS = ["STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "ARHGAP32", "PFKP", "ATG4B",
       "KALRN", "CAMK2B", "RSF1", "SETD5", "ATG4B"]
SOCE = ["STIM1", "STIM2", "ORAI1", "ORAI2", "ORAI3", "TRPC1", "TRPC3", "TRPC4", "TRPC5",
        "TRPC6", "SARAF", "STIMATE", "CBARP", "ATP2A1", "ATP2A2", "ATP2A3", "ATP2B1",
        "ATP2B4", "MCU", "MCUB", "MICU1", "MICU2", "MICU3", "SELENOK", "SELENON",
        "CRACR2A", "CRACR2B", "SEPTIN4", "ITPR1", "ITPR2", "ITPR3", "RYR1", "RYR2", "RYR3"]
def rapor(adlar, etiket):
    d = res[res["sembol"].isin(adlar)][
        ["sembol", "etkilesim_log2", "n_pozitif", "p", "q", "sira", "yuzdelik"]]
    d = d.sort_values("etkilesim_log2", ascending=False)
    d.to_csv(f"{OUT}/NMD_{etiket}.tsv", sep="\t", index=False)
    print(f"\n=== {etiket} ===")
    print(d.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    return d
pos = rapor(POS, "pozitif_kontroller")
soce = rapor(SOCE, "SOCE_paneli")

# pozitif kontrollerin genel dağılıma göre konumu
allv = res["etkilesim_log2"].values
pv = res[res["sembol"].isin(POS)]["etkilesim_log2"].values
u = stats.mannwhitneyu(pv, allv, alternative="greater")
print(f"\nPozitif kontrol paneli genel dağılıma karşı: Mann-Whitney p = {u.pvalue:.3g} "
      f"(medyan {np.median(pv):+.3f} vs {np.median(allv):+.3f})")

# Ca paneli toptan
P = "/Users/elmas/Desktop/TEZ/output/calcium_panel_overlap_analysis_2026-04-25/01_gene_panels"
for tier in ["Tier1_SOCE_TRP_51", "Tier2_Channel_Release_Transport_117",
             "Tier3_Curated_Calcium_Handling_258", "Tier4_Expanded_Calcium_Associated_732"]:
    gl = set(pd.read_csv(f"{P}/{tier}.csv")["gene_upper"])
    v = res[res["sembol"].isin(gl)]["etkilesim_log2"].values
    if len(v) < 5: continue
    uu = stats.mannwhitneyu(v, allv, alternative="greater")
    print(f"{tier}: n={len(v)} medyan {np.median(v):+.4f} vs {np.median(allv):+.4f}, p={uu.pvalue:.3g}")

#!/usr/bin/env python3
"""A8 — Kendi kriptik gen kümemizin NMD duyarlılığı.

A1'de SH-SY5Y'de de novo bulunan kriptik birleşimleri taşıyan genler, bağımsız
bir NMD inhibisyon deneyinde (GSE307054) TDP-43'e ÖZGÜ NMD etkileşimi
göstermeli. Bu, kriptik çağrılarımızın dış geçerliliğini sınar ve hangi
SOCE geninin kriptik-NMD imzası taşıdığını belirler.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
import lib_junc as L

OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
nmd = pd.read_csv(f"{OUT}/NMD_etkilesim_tum_genler.tsv", sep="\t")
kr = pd.read_csv(f"{OUT}/KRIPTIK_olaylar_SH_SY5Y.tsv", sep="\t")

# kriptik gen kümesi: etkiye göre katmanlandır
kr_gen = kr.groupby("gene").agg(max_dPSI=("dPSI", "max"),
                                n_olay=("dPSI", "size"),
                                min_q=("q", "min"),
                                okuma=("okuma_KD", "max")).reset_index()
kr_gen = kr_gen[~kr_gen["gene"].str.startswith("ENSG")]
kr_gen.to_csv(f"{OUT}/KRIPTIK_genler_SH_SY5Y.tsv", sep="\t", index=False)
print(f"SH-SY5Y kriptik gen: {len(kr_gen)}")

m = nmd.merge(kr_gen, left_on="sembol", right_on="gene", how="left")
m["kriptik"] = m["gene"].notna()
arka = m[~m["kriptik"]]["etkilesim_log2"].values
on = m[m["kriptik"]]["etkilesim_log2"].values
u = stats.mannwhitneyu(on, arka, alternative="greater")
print(f"\nKriptik gen taşıyanlar (n={len(on)}) vs diğerleri (n={len(arka)}):")
print(f"  medyan NMD etkileşimi {np.median(on):+.4f} vs {np.median(arka):+.4f}")
print(f"  Mann-Whitney tek yönlü p = {u.pvalue:.3g}")

# etki büyüklüğüne göre katman
for lo, hi, ad in [(0.05, 0.2, "zayıf"), (0.2, 0.5, "orta"), (0.5, 1.01, "güçlü")]:
    v = m[m["kriptik"] & (m["max_dPSI"] >= lo) & (m["max_dPSI"] < hi)]["etkilesim_log2"].values
    if len(v) < 5: continue
    uu = stats.mannwhitneyu(v, arka, alternative="greater")
    print(f"  ΔPSI {lo}-{hi} ({ad}): n={len(v)}, medyan {np.median(v):+.4f}, p={uu.pvalue:.3g}")

# kriptik + NMD-duyarlı kesişimi
kesisim = m[m["kriptik"] & (m["etkilesim_log2"] > 0) & (m["q"] < 0.10)]
kesisim = kesisim.sort_values("etkilesim_log2", ascending=False)
kesisim[["sembol", "max_dPSI", "n_olay", "okuma", "etkilesim_log2", "n_pozitif", "p", "q", "yuzdelik"]] \
    .to_csv(f"{OUT}/KRIPTIK_ve_NMD_duyarli_genler.tsv", sep="\t", index=False)
print(f"\nHem kriptik hem NMD-duyarlı (etkileşim>0, q<0.10): {len(kesisim)} gen")
print(kesisim[["sembol", "max_dPSI", "okuma", "etkilesim_log2", "q", "yuzdelik"]]
      .head(25).to_string(index=False, float_format=lambda x: f"{x:.4g}"))

#!/usr/bin/env python3
"""A11 — Veri setleri arası kriptik özet, RBP özgüllüğü ve STIM2 SOAR ekzonu."""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import lib_junc as L
OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
TAB = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/tablolar"; os.makedirs(TAB, exist_ok=True)

POS = ["STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
       "ATG4B", "SETD5", "ELAVL3", "POLDIP3", "CAMK2B", "RSF1", "GPSM2", "SYNJ2"]
POSU = {p.upper() for p in POS}
SIRA = ["SH_SY5Y", "SH_SY5Y_DOZ25", "iPSC_koloni", "iPSC_MN", "K562_mRNA",
        "K562_totalRNA", "C2C12", "NSC34", "Fare_striatum", "iPSC_MN_FUS", "iPSC_MN_TAF15"]

# ---------------------------------------------------- 1. kriptik özet + RBP özgüllüğü
rows = []
for f in glob.glob(f"{OUT}/RT_KRIPTIK_*.tsv"):
    ds = os.path.basename(f)[11:-4]
    d = pd.read_csv(f, sep="\t")
    d["gen_u"] = d["gene"].astype(str).str.upper()
    rows.append(dict(veri_seti=ds, kriptik_olay=len(d), kriptik_gen=d["gene"].nunique(),
                     pozitif_kontrol_gen=d.loc[d["gen_u"].isin(POSU), "gene"].nunique(),
                     bulunan=",".join(sorted(set(d.loc[d["gen_u"].isin(POSU), "gen_u"])))))
S = pd.DataFrame(rows)
S["sira"] = S["veri_seti"].apply(lambda x: SIRA.index(x) if x in SIRA else 99)
S = S.sort_values("sira").drop(columns="sira")
S.to_csv(f"{TAB}/S13_kriptik_ozet_ve_RBP_ozgullugu.tsv", sep="\t", index=False)
print("=== Kriptik olay özeti ve RBP özgüllüğü ===")
print(S.to_string(index=False))

# ---------------------------------------------------- 2. STIM2 SOAR ekzonu (24 nt)
# İnsan: chr4:27.007.982-27.008.006 ; fare: chr5:54.110.115-54.110.139
HED = {"human": ("4", 27007982, 27008006, "STIM2"), "mouse": ("5", 54110115, 54110139, "Stim2")}
rows = []
for f in glob.glob(f"{OUT}/RT_LSV_*.tsv.gz"):
    ds = os.path.basename(f)[7:-7]
    d = pd.read_csv(f, sep="\t", low_memory=False)
    sp = "mouse" if ds in ("C2C12", "NSC34", "Fare_striatum") else "human"
    c, e_s, e_e, gen = HED[sp]
    sub = d[(d["gene"] == gen) & (d["chrom"].astype(str) == c)]
    # ekzonun iki yanındaki birleşimler: end == e_s-1 (yukarı intron) veya start == e_e+1 (aşağı intron)
    inc = sub[(sub["end"] == e_s - 1) | (sub["start"] == e_e + 1)]
    for _, r in inc.iterrows():
        rows.append(dict(veri_seti=ds, yan="yukari" if r["end"] == e_s - 1 else "asagi",
                         chrom=r["chrom"], start=r["start"], end=r["end"], sinif=r["sinif"],
                         lsv_tipi=r["lsv_tipi"], psi_KD=r["psi_KD"], psi_CTRL=r["psi_CTRL"],
                         dPSI=r["dPSI"], GA_alt=r["GA_alt"], GA_ust=r["GA_ust"], q=r["q"],
                         okuma_KD=r["okuma_KD"], okuma_CTRL=r["okuma_CTRL"],
                         toplam_KD=r["toplam_KD"], toplam_CTRL=r["toplam_CTRL"]))
T = pd.DataFrame(rows)
if len(T):
    T.to_csv(f"{TAB}/S14_STIM2_SOAR_ekzonu_birlesim_duzeyi.tsv", sep="\t", index=False)
    print("\n=== STIM2 SOAR ekzonu (24 nt) — birleşim düzeyi ===")
    print(T.sort_values(["veri_seti", "yan"]).to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    # havuzlanmış işaret testi (yalnızca TDP-43 karşılaştırmaları)
    tdp = T[~T["veri_seti"].isin(["iPSC_MN_FUS", "iPSC_MN_TAF15"])]
    ic = tdp[tdp["lsv_tipi"].isin(["verici", "alici"])]
    print(f"\nTDP-43 karşılaştırmalarında ekzon-yanı birleşim sayısı: {len(ic)}, "
          f"pozitif ΔPSI: {(ic['dPSI'] > 0).sum()}, medyan ΔPSI: {ic['dPSI'].median():+.4f}")

#!/usr/bin/env python3
"""A10 — Makale için ek tabloların üretimi."""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
TAB = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/tablolar"; os.makedirs(TAB, exist_ok=True)

def yaz(df, ad, **kw):
    df.to_csv(f"{TAB}/{ad}.tsv", sep="\t", index=False, **kw)
    print(f"{ad}: {len(df)} satır")

# --- S5: kriptik pozitif kontroller (tüm veri setleri)
rows = []
for f in sorted(glob.glob(f"{OUT}/YUKSEK_GUVEN_*.tsv")):
    ds = os.path.basename(f)[len("YUKSEK_GUVEN_"):-len(".tsv")]
    d = pd.read_csv(f, sep="\t")
    d["veri_seti"] = ds
    rows.append(d)
if rows:
    K = pd.concat(rows, ignore_index=True)
    K = K[["veri_seti", "gene", "chrom", "start", "end", "strand", "sinif", "lsv_tipi",
           "psi_KD", "psi_CTRL", "dPSI", "GA_alt", "GA_ust", "q",
           "okuma_KD", "okuma_CTRL", "toplam_KD", "toplam_CTRL"]]
    yaz(K, "S5_yuksek_guven_kriptik_olaylar")
    # the sixteen literature cryptic targets, as defined in A2_lsv_sorgu.py (POS_H),
    # with the mouse spellings of the same genes
    POS = ["STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
           "ATG4B", "SETD5", "CAMK2B", "ELAVL3", "POLDIP3", "RSF1", "GPSM2", "SYNJ2"]
    POS = POS + [g.capitalize() for g in POS]
    yaz(K[K["gene"].isin(POS)], "S6_kriptik_pozitif_kontroller")
    # veri seti x gen matrisi
    K["gen_u"] = K["gene"].str.upper()
    piv = (K[K["gen_u"].isin([p.upper() for p in POS])]
           .pivot_table(index="gen_u", columns="veri_seti", values="dPSI", aggfunc="max"))
    piv.round(3).to_csv(f"{TAB}/S6b_kriptik_pozitif_kontrol_matrisi.tsv", sep="\t")
    print("S6b matrisi yazıldı")

# --- S7: SOCE panelinde anotasyondan bağımsız olaylar
rows = []
for f in sorted(glob.glob(f"{OUT}/RT_LSV_*.tsv.gz")):
    ds = os.path.basename(f)[7:-7]
    d = pd.read_csv(f, sep="\t", low_memory=False)
    SOCE = ["STIM1", "STIM2", "ORAI1", "ORAI2", "ORAI3", "TRPC1", "SARAF", "STIMATE",
            "CBARP", "ATP2A1", "ATP2A2", "ATP2A3", "MCU", "MCUB", "MICU1", "MICU2",
            "MICU3", "CRACR2A", "CRACR2B", "SELENOK", "SELENON"]
    SOCE = SOCE + [g.capitalize() for g in SOCE]
    d = d[d["gene"].isin(SOCE) & (d["q"] < 0.05) & (d["dPSI"].abs() >= 0.10)]
    if len(d):
        d = d.copy(); d["veri_seti"] = ds; rows.append(d)
if rows:
    S = pd.concat(rows, ignore_index=True)
    S = S[["veri_seti", "gene", "chrom", "start", "end", "strand", "sinif", "lsv_tipi",
           "psi_KD", "psi_CTRL", "dPSI", "GA_alt", "GA_ust", "q",
           "okuma_KD", "okuma_CTRL", "toplam_KD", "toplam_CTRL"]]
    yaz(S.sort_values(["gene", "veri_seti"]), "S7_SOCE_anotasyonsuz_analiz")

# --- S8: NYGC kriptik PSI
for src, ad in ((f"{OUT}/NYGC_kriptik_PSI_ALS_vs_kontrol.tsv", "S8_NYGC_kriptikPSI_ALS_vs_kontrol"),
                (f"{OUT}/NYGC_kriptikPSI_korelasyon.tsv", "S9_NYGC_kriptikPSI_korelasyon")):
    if os.path.exists(src):
        yaz(pd.read_csv(src, sep="\t"), ad)

# --- S10: NMD etkileşimi (dependence-aware, four condition-level units)
if os.path.exists(f"{OUT}/NMD_etkilesim_paylasimli_kontrol_SOCE_pozitif.tsv"):
    yaz(pd.read_csv(f"{OUT}/NMD_etkilesim_paylasimli_kontrol_SOCE_pozitif.tsv", sep="\t"),
        "S10_NMD_etkilesimi_SOCE")
if os.path.exists(f"{OUT}/KRIPTIK_ve_NMD_duyarli_genler.tsv"):
    yaz(pd.read_csv(f"{OUT}/KRIPTIK_ve_NMD_duyarli_genler.tsv", sep="\t"),
        "S11_kriptik_ve_NMD_duyarli")

# --- S12: corrected SH-SY5Y APA candidate gradients
if os.path.exists(f"{OUT}/APA_corrected_full_core_summary.tsv"):
    A = pd.read_csv(f"{OUT}/APA_corrected_full_core_summary.tsv", sep="\t")
    yaz(A[A["delta"].abs() >= 0.05].sort_values(["gene", "unit"]),
        "S12_APA_anlamli_olaylar")

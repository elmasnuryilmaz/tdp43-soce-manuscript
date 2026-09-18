#!/usr/bin/env python3
"""A2 — LSV sonuçlarının hedefli sorgusu: pozitif kontroller, SOCE paneli, kriptik olaylar."""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import lib_junc as L

OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
P = "/Users/elmas/Desktop/TEZ/output/calcium_panel_overlap_analysis_2026-04-25/01_gene_panels"

POS_H = ["STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
         "ATG4B", "SETD5", "CAMK2B", "ELAVL3", "POLDIP3", "RSF1", "GPSM2", "SYNJ2"]
SOCE_H = ["STIM1", "STIM2", "ORAI1", "ORAI2", "ORAI3", "TRPC1", "SARAF", "STIMATE",
          "CBARP", "ATP2A1", "ATP2A2", "ATP2A3", "MCU", "MCUB", "MICU1", "MICU2",
          "MICU3", "CRACR2A", "CRACR2B", "SELENOK", "SELENON", "ITPR1", "ITPR2", "ITPR3"]
def fare(x): return [g.capitalize() for g in x]

KOL = ["gene", "chrom", "start", "end", "strand", "sinif", "lsv_tipi", "n_junction_lsv",
       "psi_KD", "psi_CTRL", "dPSI", "GA_alt", "GA_ust", "p", "q",
       "okuma_KD", "okuma_CTRL", "toplam_KD", "toplam_CTRL",
       "n_KD_pozitif", "n_CTRL_pozitif"]

def yukle(ds):
    f = f"{OUT}/LSV_{ds}.tsv.gz"
    return pd.read_csv(f, sep="\t") if os.path.exists(f) else None

def rapor(ds, genler, etiket):
    D = yukle(ds)
    if D is None: return None
    d = D[D["gene"].isin(genler)][KOL].copy()
    d = d.sort_values(["gene", "q"])
    d.to_csv(f"{OUT}/LSV_{etiket}_{ds}.tsv", sep="\t", index=False)
    return d

def kriptik_cagri(ds, min_dpsi=0.05):
    """Kriptik olay: yeni birleşim, KD'de daha yüksek PSI, kontrolde yok denecek kadar az."""
    D = yukle(ds)
    if D is None: return None
    k = D[(D["sinif"] != "anotasyonlu") & (D["dPSI"] >= min_dpsi) & (D["q"] < 0.05)
          & (D["psi_CTRL"] <= 0.05) & (D["GA_alt"] > 0)].copy()
    k = k.sort_values("dPSI", ascending=False)
    k[KOL].to_csv(f"{OUT}/KRIPTIK_olaylar_{ds}.tsv", sep="\t", index=False)
    return k

def main():
    datasets = sorted({os.path.basename(f)[4:-7] for f in glob.glob(f"{OUT}/LSV_*.tsv.gz")})
    print("veri setleri:", datasets)
    ozet = []
    for ds in datasets:
        insan = ds in ("SH_SY5Y", "iPSC_koloni", "iPSC_MN", "K562_mRNA", "K562_totalRNA")
        pos = POS_H if insan else fare(POS_H)
        soce = SOCE_H if insan else fare(SOCE_H)
        k = kriptik_cagri(ds)
        dp = rapor(ds, pos, "pozitif_kontrol")
        dsoce = rapor(ds, soce, "SOCE")
        n_kriptik_gen = k["gene"].nunique() if k is not None else 0
        ozet.append(dict(veri_seti=ds, kriptik_olay=len(k) if k is not None else 0,
                         kriptik_gen=n_kriptik_gen,
                         pozitif_kontrol_geni=int(k["gene"].isin(pos).sum()) if k is not None else 0))
        print(f"\n{'='*70}\n### {ds}: kriptik olay {len(k)}, gen {n_kriptik_gen}")
        if k is not None and len(k):
            print(k[["gene","chrom","start","end","sinif","psi_KD","psi_CTRL","dPSI","q",
                     "okuma_KD","okuma_CTRL"]].head(15).to_string(index=False,
                     float_format=lambda x: f"{x:.4g}"))
        if dp is not None and len(dp):
            sig = dp[(dp["q"] < 0.05) & (dp["dPSI"].abs() >= 0.05)]
            print(f"\n-- pozitif kontrol genlerinde anlamlı olay: {len(sig)} / {len(dp)} test")
            if len(sig): print(sig.head(12).to_string(index=False, float_format=lambda x: f"{x:.4g}"))
        if dsoce is not None and len(dsoce):
            sig = dsoce[(dsoce["q"] < 0.05) & (dsoce["dPSI"].abs() >= 0.10)]
            print(f"\n-- SOCE panelinde anlamlı olay: {len(sig)} / {len(dsoce)} test")
            if len(sig): print(sig.head(20).to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    pd.DataFrame(ozet).to_csv(f"{OUT}/LSV_ozet.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""A13 — Kriptik çağrı ölçütünün kalibrasyonu ve yanlış pozitif oranı.

Kontrol-karşı-kontrol boş testi, gevşek ölçütle yanlış pozitif oranının bazı veri
setlerinde gerçek çağrı sayısını aştığını gösterdi. Ölçüt burada ETKİ TEMELLİ olarak
yeniden tanımlanır (kontrolde sıfır olma koşulu, STMN2 gibi bazal sızıntısı olan
gerçek kriptik ekzonları elediği için kaldırılmıştır):

  YÜKSEK GÜVEN = anotasyonsuz kırpılma bölgesi (yeni_bolge / tamamen_yeni)
                 + ΔPSI ≥ 0,20 + q < 0,05
                 + bootstrap %95 GA alt sınırı > 0,05
                 + KD'de ≥ 20 okuma
                 + tüm KD replikalarında sıfırdan büyük

Her kademede boş test ile yanlış pozitif oranı ve pozitif kontrol geri kazanımı
birlikte raporlanır; doğru eşik, yanlış pozitifi düşürürken pozitif kontrolleri
korumayı sürdüren eşiktir.
"""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd

D = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI"
OUT, TAB = f"{D}/sonuclar", f"{D}/tablolar"
POS = {"STMN2","UNC13A","HDGFL2","ACTL6B","AGRN","KALRN","ARHGAP32","PFKP","ATG4B",
       "SETD5","ELAVL3","POLDIP3","CAMK2B","RSF1","GPSM2","SYNJ2"}
CEK = {"STIM1","STIM2","ORAI1","ORAI2","ORAI3","TRPC1","SARAF","STIMATE","CBARP","ATP2A1",
       "ATP2A2","ATP2A3","MCU","MCUB","MICU1","MICU2","CRACR2A","CRACR2B","SELENOK"}
P = "/Users/elmas/Desktop/TEZ/output/calcium_panel_overlap_analysis_2026-04-25/01_gene_panels"
T1 = set(pd.read_csv(f"{P}/Tier1_SOCE_TRP_51.csv")["gene_upper"])
MAN = pd.read_csv(f"{D}/kod/ornekler_regtools.tsv", sep="\t")
KADEME = [("K1", 0.10, 0.00, 10, False), ("K2", 0.20, 0.05, 20, True), ("K3", 0.30, 0.10, 30, True)]

def sup(d, dp, ga, ok, tum, n_kd):
    m = (d["sinif"].isin(["yeni_bolge", "tamamen_yeni"])) & (d["dPSI"] >= dp) \
        & (d["q"] < 0.05) & (d["GA_alt"] > ga) & (d["okuma_KD"] >= ok)
    if tum:
        m &= d["n_KD_pozitif"] >= n_kd
    return d[m]

# ------------------------------------------------- 1. kalibrasyon
rows = []
for ds in ["iPSC_koloni", "K562_totalRNA", "Fare_striatum"]:
    fn = f"{OUT}/NULL_KRIPTIK_{ds}.tsv"
    if not os.path.exists(fn): continue
    G = pd.read_csv(f"{OUT}/RT_LSV_{ds}.tsv.gz", sep="\t", low_memory=False)
    N = pd.read_csv(fn, sep="\t")
    n_kd = MAN[(MAN.dataset == ds) & (MAN.group == "KD")]["replika"].nunique()
    for ad, dp, ga, ok, tum in KADEME:
        g = sup(G, dp, ga, ok, tum, n_kd); n = sup(N, dp, ga, ok, tum, max(n_kd // 2, 1))
        rows.append(dict(veri_seti=ds, kademe=ad, olcut=f"ΔPSI≥{dp}, GA>{ga}, ≥{ok} okuma",
                         gercek_olay=len(g), bos_olay=len(n),
                         yanlis_pozitif_orani=round(len(n) / max(len(g), 1), 2),
                         pozitif_kontrol=len(set(g.gene.astype(str).str.upper()) & POS)))
K = pd.DataFrame(rows)
K.to_csv(f"{TAB}/S15_kriptik_esik_kalibrasyonu.tsv", sep="\t", index=False)
print("=== Eşik kalibrasyonu (boş test: kontrol karşı kontrol) ===")
print(K.to_string(index=False))

# ------------------------------------------------- 2. seçilen ölçüt tüm veri setlerine
rows2 = []
for f in sorted(glob.glob(f"{OUT}/RT_LSV_*.tsv.gz")):
    ds = os.path.basename(f)[7:-7]
    n_kd = MAN[(MAN.dataset == ds) & (MAN.group == "KD")]["replika"].nunique()
    d = sup(pd.read_csv(f, sep="\t", low_memory=False), 0.20, 0.05, 20, True, n_kd)
    d.to_csv(f"{OUT}/YUKSEK_GUVEN_{ds}.tsv", sep="\t", index=False)
    gu = set(d.gene.astype(str).str.upper())
    rows2.append(dict(veri_seti=ds, olay=len(d), gen=d.gene.nunique(),
                      pozitif_kontrol=len(gu & POS),
                      bulunan=",".join(sorted(gu & POS)) or "—",
                      Tier1=",".join(sorted(gu & T1)) or "—",
                      cekirdek_SOCE=",".join(sorted(gu & CEK)) or "—"))
S = pd.DataFrame(rows2)
S.to_csv(f"{TAB}/S13_yuksek_guven_kriptik_ozet.tsv", sep="\t", index=False)
print("\n=== Yüksek güven ölçütüyle (K2) ===")
print(S.to_string(index=False))

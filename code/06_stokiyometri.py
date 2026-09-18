#!/usr/bin/env python3
"""
Oncelik 3.3 — SOCE bilesenlerinin MUTLAK stokiyometrisi.

Tez, log2FC uzerinden "bilesenler artti ama SOCE azaldi" diyor. Asil soru,
hangi izoformun BASKIN oldugu: baskin izoform azalirken minor izoform artiyorsa
toplam fonksiyon dusebilir. Bunu mutlak ifade duzeyinde gosteriyoruz.
"""
import os
import numpy as np
import pandas as pd

BASE = "/Users/elmas/Desktop/TEZ"
OUT = "/Users/elmas/Desktop/MAKALE"
TAB = os.path.join(OUT, "03_TABLOLAR")

CTRL = ["SRR33374996", "SRR33375001", "SRR33374995"]   # 0 ng/mL
KD   = ["SRR33374999", "SRR33374997", "SRR33375000"]   # 75 ng/mL

cnt = pd.read_csv(os.path.join(BASE, "output/GSE296712_doz_DEG_22.08.26/gene_counts_salmon.csv"),
                  index_col=0)
cnt.index = cnt.index.astype(str).str.upper()
cnt = cnt[CTRL + KD]

# CPM benzeri normalizasyon (kutuphane buyuklugu)
lib = cnt.sum(0)
cpm = cnt.divide(lib, axis=1) * 1e6

deg = pd.read_csv(os.path.join(BASE, "output/GSE296712_doz_DEG_22.08.26/DESeq2_ctrl_vs_75.csv"))
deg["gene"] = deg["gene"].astype(str).str.upper()
deg = deg.drop_duplicates("gene").set_index("gene")

FAMILIES = {
    "STIM (ER Ca sensoru)":        ["STIM1", "STIM2"],
    "ORAI (CRAC kanali)":          ["ORAI1", "ORAI2", "ORAI3"],
    "SERCA (ER'ye Ca geri alimi)": ["ATP2A1", "ATP2A2", "ATP2A3"],
    "TRPC":                        ["TRPC1", "TRPC3", "TRPC4", "TRPC5", "TRPC6"],
    "SOCE duzenleyicileri":        ["SARAF", "STIMATE", "CRACR2A", "CRACR2B", "CBARP"],
    "Mitokondriyal Ca alimi":      ["MCU", "MICU1", "MICU2", "MICU3", "MCUR1", "MCUB"],
    "PMCA (Ca disari atim)":       ["ATP2B1", "ATP2B2", "ATP2B3", "ATP2B4"],
}

rows = []
for fam, genes in FAMILIES.items():
    present = [g for g in genes if g in cpm.index]
    fam_ctrl_total = cpm.loc[present, CTRL].mean(1).sum() if present else np.nan
    for g in genes:
        if g not in cpm.index:
            rows.append(dict(aile=fam, gen=g, durum="veride yok"))
            continue
        c = cpm.loc[g, CTRL].mean()
        k = cpm.loc[g, KD].mean()
        l2 = deg.loc[g, "log2FoldChange"] if g in deg.index else np.nan
        pa = deg.loc[g, "padj"] if g in deg.index else np.nan
        rows.append(dict(
            aile=fam, gen=g, durum="",
            CPM_kontrol=round(c, 2), CPM_KD=round(k, 2),
            ailedeki_pay_kontrol_pct=round(100 * c / fam_ctrl_total, 1) if fam_ctrl_total else np.nan,
            baskin_mi=("BASKIN" if (fam_ctrl_total and c / fam_ctrl_total > 0.5) else ""),
            log2FC=round(l2, 3) if pd.notna(l2) else np.nan,
            padj=f"{pa:.2e}" if pd.notna(pa) else "NA",
            yon=("artis" if pd.notna(l2) and l2 > 0 else ("azalis" if pd.notna(l2) else "NA")),
            mutlak_degisim_CPM=round(k - c, 2),
        ))

T = pd.DataFrame(rows)
T.to_csv(os.path.join(TAB, "SOCE_mutlak_stokiyometri.tsv"), sep="\t", index=False)

pd.set_option("display.width", 220)
print("=" * 110)
print("SOCE MAKINESININ MUTLAK STOKIYOMETRISI — SH-SY5Y (GSE296712, 0 vs 75 ng/mL)")
print("CPM = kutuphane buyuklugune gore normalize edilmis Salmon gen sayimlari")
print("=" * 110)
for fam in FAMILIES:
    sub = T[T.aile == fam]
    print(f"\n### {fam}")
    show = sub[sub.durum == ""].drop(columns=["aile", "durum"])
    if show.empty:
        print("   (veride yok)")
        continue
    print(show.to_string(index=False))

# ---------------------------------------------------------------- net etki
print("\n" + "=" * 110)
print("AILE DUZEYINDE NET DEGISIM (mutlak CPM toplamlari)")
print("=" * 110)
net = []
for fam, genes in FAMILIES.items():
    present = [g for g in genes if g in cpm.index]
    if not present:
        continue
    c = cpm.loc[present, CTRL].mean(1).sum()
    k = cpm.loc[present, KD].mean(1).sum()
    net.append(dict(aile=fam, toplam_CPM_kontrol=round(c, 1), toplam_CPM_KD=round(k, 1),
                    net_degisim=round(k - c, 1),
                    yuzde_degisim=round(100 * (k - c) / c, 1) if c else np.nan))
N = pd.DataFrame(net)
print(N.to_string(index=False))
N.to_csv(os.path.join(TAB, "SOCE_aile_duzeyi_net_degisim.tsv"), sep="\t", index=False)

print("\n" + "=" * 110)
print("YORUM ICIN KRITIK ORANLAR")
print("=" * 110)
def ratio(a, b, label):
    if a in cpm.index and b in cpm.index:
        rc = cpm.loc[a, CTRL].mean() / max(cpm.loc[b, CTRL].mean(), 1e-9)
        rk = cpm.loc[a, KD].mean() / max(cpm.loc[b, KD].mean(), 1e-9)
        print(f"  {label:28s} kontrol = {rc:8.2f}   KD = {rk:8.2f}   degisim = {rk-rc:+.2f}")
ratio("ORAI2", "ORAI1", "ORAI2 / ORAI1")
ratio("ATP2A2", "ATP2A3", "ATP2A2 / ATP2A3 (SERCA2/3)")
ratio("STIM1", "ORAI1", "STIM1 / ORAI1")
ratio("STIM2", "STIM1", "STIM2 / STIM1")
ratio("TRPC1", "ORAI1", "TRPC1 / ORAI1")
print("\nTablolar ->", TAB)

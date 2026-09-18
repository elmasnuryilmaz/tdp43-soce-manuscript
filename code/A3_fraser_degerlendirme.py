#!/usr/bin/env python3
"""A3 — FRASER çıktısının gen anotasyonlu yeniden değerlendirmesi.

Temmuz 2026'da SH-SY5Y doksisiklin serisinde (9 örnek) çalıştırılan FRASER
sonuçları gen sembolü olmadan bırakılmış, anlamlı olay bulunamamıştı.
Burada: (1) her aykırı olay gene atanır, (2) örnek başına aykırı yükü
TDP-43 tükenmiş ve kontrol örnekleri arasında karşılaştırılır, (3) SOCE
paneli ve kriptik pozitif kontrol lokusları özel olarak sorgulanır.

FRASER bir AYKIRI DEĞER yöntemidir: kohortun beklenen kırpılma örüntüsünü
öğrenip sapan ÖRNEKLERİ işaretler. 9 örnekli ve KD'nin çoğunluk olduğu bir
kohortta grup düzeyi güç düşüktür; sonuç bu sınır içinde raporlanır.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_junc as L

F = "/Volumes/10TBElmas/thesis_addendum_2026/03_fraser/FRASER_SHSY5Y_all_results.tsv"
OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
GRUP = {"SRR33374996": "CTRL", "SRR33375001": "CTRL", "SRR33374995": "CTRL",
        "SRR33374939": "KD25", "SRR33375002": "KD25", "SRR33374998": "KD25",
        "SRR33374999": "KD75", "SRR33374997": "KD75", "SRR33375000": "KD75"}

use = ["seqnames", "start", "end", "strand", "sampleID", "type",
       "pValue", "padjust", "zScore", "psiValue", "deltaPsi",
       "counts", "totalCounts"]
print("FRASER tablosu okunuyor (427 MB)...", flush=True)
it = pd.read_csv(F, sep="\t", usecols=use, dtype={"seqnames": str}, chunksize=2_000_000)
keep = []
n_tot = 0
for ch in it:
    n_tot += len(ch)
    keep.append(ch[ch["pValue"] < 1e-3])
df = pd.concat(keep, ignore_index=True)
print(f"toplam satır {n_tot:,}; p<1e-3 olan {len(df):,}", flush=True)
df["grup"] = df["sampleID"].map(GRUP)
df["grup2"] = np.where(df["grup"] == "CTRL", "CTRL", "KD")

# --- örnek başına aykırı yükü
for thr in (1e-3, 1e-4, 1e-5, 1e-6):
    sub = df[(df["pValue"] < thr) & (df["deltaPsi"].abs() >= 0.10)]
    burden = sub.groupby(["sampleID"]).size().reindex(GRUP.keys()).fillna(0).astype(int)
    b = pd.DataFrame({"n_aykiri": burden, "grup": [GRUP[s] for s in burden.index]})
    kd = b[b.grup != "CTRL"]["n_aykiri"].values; ct = b[b.grup == "CTRL"]["n_aykiri"].values
    u = stats.mannwhitneyu(kd, ct, alternative="greater") if len(ct) and len(kd) else None
    print(f"\np<{thr:g}, |dPsi|>=0.10 — örnek başına aykırı olay")
    print(b.to_string())
    if u: print(f"KD ortalama {kd.mean():.1f} vs CTRL {ct.mean():.1f}, Mann-Whitney p = {u.pvalue:.3f}")
    if thr == 1e-5:
        b.to_csv(f"{OUT}/FRASER_ornek_basina_aykiri_yuk.tsv", sep="\t")

# --- gen ataması
ai, genes = L.load_annotation("human")
sig = df[(df["pValue"] < 1e-5) & (df["deltaPsi"].abs() >= 0.10)].copy()
gidx = {}
for c, s, e, g in zip(genes["chrom"], genes["start"], genes["end"], genes["gene"]):
    gidx.setdefault(c, []).append((s, e, g))
for c in gidx: gidx[c].sort()
starts = {c: np.array([x[0] for x in v]) for c, v in gidx.items()}
def gene_of(c, s):
    if c not in gidx: return "."
    i = np.searchsorted(starts[c], s, side="right")
    for j in range(max(0, i - 400), min(len(gidx[c]), i + 1)):
        gs, ge, gn = gidx[c][j]
        if gs <= s <= ge: return gn
    return "."
sig["gene"] = [gene_of(c, s) for c, s in zip(sig["seqnames"], sig["start"])]
sig.to_csv(f"{OUT}/FRASER_aykiri_olaylar_genli.tsv", sep="\t", index=False)
print(f"\np<1e-5 & |dPsi|>=0.10 olay: {len(sig)}, gen atanan: {(sig['gene']!='.').sum()}")

# --- KD'ye özgü aykırı genler
kd_only = sig[sig["grup2"] == "KD"]["gene"].value_counts()
ct_gen = set(sig[sig["grup2"] == "CTRL"]["gene"])
ozel = kd_only[[g not in ct_gen and g != "." for g in kd_only.index]]
print("\nYalnızca TDP-43 tükenmiş örneklerde aykırı olan ilk 25 gen:")
print(ozel.head(25).to_string())
ozel.to_frame("n_KD_ornegi").to_csv(f"{OUT}/FRASER_KD_ozgu_genler.tsv", sep="\t")

# --- ilgi panelleri
POS = ["STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "ARHGAP32", "PFKP", "ATG4B", "KALRN"]
SOCE = ["STIM1", "STIM2", "ORAI1", "ORAI2", "ORAI3", "TRPC1", "SARAF", "STIMATE", "CBARP",
        "ATP2A2", "ATP2A3", "MCU", "MCUB", "MICU1", "MICU2", "CRACR2A", "CRACR2B", "SELENOK"]
for etiket, panel in (("kriptik_pozitif_kontrol", POS), ("SOCE", SOCE)):
    d = sig[sig["gene"].isin(panel)][
        ["gene", "seqnames", "start", "end", "type", "sampleID", "grup",
         "pValue", "padjust", "deltaPsi", "counts", "totalCounts"]]
    d = d.sort_values(["gene", "pValue"])
    d.to_csv(f"{OUT}/FRASER_{etiket}.tsv", sep="\t", index=False)
    print(f"\n=== FRASER aykırı olayları — {etiket} ({len(d)} satır) ===")
    if len(d): print(d.head(30).to_string(index=False, float_format=lambda x: f"{x:.3g}"))

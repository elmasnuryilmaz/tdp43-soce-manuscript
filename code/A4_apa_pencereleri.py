#!/usr/bin/env python3
"""A4a — Alternatif poliadenilasyon (APA) için kapsam pencereleri kurar.

İki ölçüm:
  IPA (intronik poliadenilasyon): her intronun 5' ve 3' uçlarındaki kapsam oranı.
      Erken sonlanma varsa intronun 5' ucu 3' ucuna göre zenginleşir.
  3'UTR kullanımı (DaPars mantığı): terminal ekzonun proksimal / distal yarısı.

Pencereler transkripsiyon yönüne göre tanımlanır (şerit dikkate alınır).
"""
import os, sys, gzip
import numpy as np, pandas as pd

REF = "/Users/elmas/Desktop/TEZ/tez_duzeltmeler/reference"
P = "/Users/elmas/Desktop/TEZ/output/calcium_panel_overlap_analysis_2026-04-25/01_gene_panels"
ND = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI"
W = 500      # pencere uzunluğu (bp)
PAD = 50     # intron uçlarından kaçınılan tampon
MIN_INTRON = 1200
MIN_UTR = 400

POS = ["STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
       "ATG4B", "SETD5", "CAMK2B", "ELAVL3", "POLDIP3", "TARDBP", "GPSM2", "SYNJ2"]
tier4 = list(pd.read_csv(f"{P}/Tier4_Expanded_Calcium_Associated_732.csv")["gene_upper"])
HEDEF = set(tier4) | set(POS)

ex = pd.read_csv(f"{REF}/gencode.v47.basic.annotation.exons.txt.gz", sep="\t")
ex["gene_name"] = ex["gene_name"].astype(str)
ex = ex[ex["gene_name"].isin(HEDEF)]
ex["chrom"] = ex["chr"].str.replace("^chr", "", regex=True)
print(f"hedef gen: {len(HEDEF)}, ekzon kaydı: {len(ex)}, bulunan gen: {ex.gene_name.nunique()}")

rows = []
for (gene, chrom, strand), g in ex.groupby(["gene_name", "chrom", "strand"]):
    iv = sorted(zip(g["start"], g["end"]))
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    if len(merged) < 2:
        continue
    # --- intronlar
    for i in range(len(merged) - 1):
        istart = merged[i][1] + 1
        iend = merged[i + 1][0] - 1
        if iend - istart + 1 < MIN_INTRON:
            continue
        a5s, a5e = istart + PAD, istart + PAD + W - 1          # genomik sol
        a3s, a3e = iend - PAD - W + 1, iend - PAD               # genomik sağ
        if strand == "+":
            up, dn = (a5s, a5e), (a3s, a3e)
        else:
            up, dn = (a3s, a3e), (a5s, a5e)                     # - şeritte ters
        idx = i + 1 if strand == "+" else len(merged) - 1 - i
        rows.append((chrom, up[0], up[1], f"{gene}|intron{idx}|I5|{strand}"))
        rows.append((chrom, dn[0], dn[1], f"{gene}|intron{idx}|I3|{strand}"))
    # --- terminal ekzon 3'UTR yarıları
    term = merged[-1] if strand == "+" else merged[0]
    ln = term[1] - term[0] + 1
    if ln >= MIN_UTR:
        mid = term[0] + ln // 2
        if strand == "+":
            prox, dist = (term[0], mid - 1), (mid, term[1])
        else:
            prox, dist = (mid, term[1]), (term[0], mid - 1)
        rows.append((chrom, prox[0], prox[1], f"{gene}|termexon|Uprox|{strand}"))
        rows.append((chrom, dist[0], dist[1], f"{gene}|termexon|Udist|{strand}"))

bed = pd.DataFrame(rows, columns=["chrom", "start", "end", "isim"])
bed = bed[bed["end"] > bed["start"]]
bed["start0"] = bed["start"] - 1
bed = bed[["chrom", "start0", "end", "isim"]].sort_values(["chrom", "start0"])
bed.to_csv(f"{ND}/kod/apa_pencereleri_human.bed", sep="\t", header=False, index=False)
print(f"pencere: {len(bed)} ({bed.isim.str.split('|').str[0].nunique()} gen)")
print(bed.head(6).to_string(index=False, header=False))

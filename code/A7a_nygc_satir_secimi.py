#!/usr/bin/env python3
"""A7a — NYGC (recount3 SRP270799) birleşim matrisinden ilgi genlerinin satırlarını seç."""
import os, sys, gzip
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ND = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/nygc_junction"
P = "/Users/elmas/Desktop/TEZ/output/calcium_panel_overlap_analysis_2026-04-25/01_gene_panels"

# TDP-43 kriptik pozitif kontrolleri (literatür) + SOCE/Ca hedefleri
POS = ["STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
       "ATG4B", "SETD5", "CAMK2B", "RSF1", "ELAVL3", "POLDIP3", "GPSM2", "SYNJ2",
       "MYO18A", "CEP72", "SLC24A3", "TARDBP"]
SOCE = ["TRPC1", "SARAF", "CBARP", "STIMATE", "STIM1", "STIM2", "ORAI1", "ORAI2",
        "ORAI3", "ATP2A1", "ATP2A2", "ATP2A3", "MCU", "MCUB", "MICU1", "MICU2",
        "MICU3", "CRACR2A", "CRACR2B", "SELENOK", "SELENON", "ITPR1", "ITPR2", "ITPR3"]
MARKER = ["SNAP25", "RBFOX3", "SYT1", "SYN1", "NEFL", "NEFM", "ENO2", "MAP2", "TUBB3",
          "GAP43", "SYP", "NRGN", "CAMK2A", "SLC17A7", "GAD1", "GAD2", "GFAP", "AIF1",
          "MBP", "PLP1", "AQP4", "CX3CR1", "OLIG2", "SOX10"]
tier2 = list(pd.read_csv(f"{P}/Tier2_Channel_Release_Transport_117.csv")["gene_upper"])
HEDEF = sorted(set(POS + SOCE + MARKER + tier2))
print(f"hedef gen: {len(HEDEF)}")

genes = pd.read_csv("/Volumes/10TBElmas/thesis_addendum_2026/ref/genes_human.bed",
                    sep="\t", header=None,
                    names=["chrom", "start0", "end", "info", "n", "strand"], dtype={"chrom": str})
genes["gene"] = genes["info"].str.split("|").str[0]
g = genes[genes["gene"].isin(HEDEF)].copy()
# gen başına en geniş aralık
win = g.groupby("gene").agg(chrom=("chrom", "first"), start=("start0", "min"),
                            end=("end", "max"), strand=("strand", "first")).reset_index()
win["chrom"] = "chr" + win["chrom"].astype(str)
print(f"eşleşen gen aralığı: {len(win)}  (bulunamayan: {sorted(set(HEDEF)-set(win['gene']))[:15]})")
win.to_csv(f"{ND}/hedef_gen_araliklari.tsv", sep="\t", index=False)

# kromozoma göre aralık listesi
by_chr = {}
for _, r in win.iterrows():
    by_chr.setdefault(r["chrom"], []).append((r["start"] - 5000, r["end"] + 5000, r["gene"]))
for c in by_chr: by_chr[c].sort()

rows = []
with gzip.open(f"{ND}/SRP270799.ALL.RR.gz", "rt") as fh:
    hdr = fh.readline()
    for i, line in enumerate(fh, start=1):        # 1-tabanlı satır = MM satır indeksi
        f = line.rstrip("\n").split("\t")
        c = f[0]
        iv = by_chr.get(c)
        if not iv: continue
        s = int(f[1]); e = int(f[2])
        for a, b, gn in iv:
            if a <= s and e <= b:
                rows.append((i, c, s, e, f[4], f[5], f[6], f[7], f[8], f[9], gn))
                break
R = pd.DataFrame(rows, columns=["satir", "chrom", "start", "end", "strand", "annotated",
                                "left_motif", "right_motif", "left_annotated",
                                "right_annotated", "gene"])
R.to_csv(f"{ND}/hedef_satirlar.tsv", sep="\t", index=False)
print(f"seçilen birleşim satırı: {len(R)}")
print(R.groupby("gene").size().sort_values(ascending=False).head(12).to_string())
print("\nSTMN2 birleşimleri (ilk 15, koordinata göre):")
print(R[R.gene == "STMN2"].head(15).to_string(index=False))

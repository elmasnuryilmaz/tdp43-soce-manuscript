#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Transcript-family abundance of SOCE-related genes in SH-SY5Y (Table 4, Figure 7).

TPM sums to the same total in every library, so a large gain by a few abundant transcripts
lowers the TPM of every other gene. In the TDP-43 knockdown libraries CHGA, ribosomal-protein
and mitochondrially encoded transcripts take a larger share of the total, and the median
expressed gene has about a quarter less TPM than in the controls although DESeq2, which
normalises for library composition, finds it unchanged. Between-condition changes are
therefore computed after rescaling each library by a median-of-ratios factor over the genes
with TPM > 1 in all six libraries (the normalisation principle of DESeq2). Shares within a
family are computed from the same adjusted values; within a library the scaling cancels.

Reads   03_TABLOLAR/v3/gene_TPM_full.csv, 03_TABLOLAR/v3/DESeq2_ctrl_vs_75_fullmap.csv
Writes  03_TABLOLAR/v3/Table4_family_TPM.csv            (Table 4, Figure 7, graphical abstract)
        03_TABLOLAR/v3/Table4_composition_scaling.csv   (scaling factors and the global shift)
        and the same two files in 09_YAYIN_PAKETI/source_data/
"""
import os
import numpy as np
import pandas as pd

M = "/Users/elmas/Desktop/MAKALE"
V3 = f"{M}/03_TABLOLAR/v3"
SRC = f"{M}/09_YAYIN_PAKETI/source_data"
CTRL = ["SRR33374996", "SRR33375001", "SRR33374995"]      # 0 ng/mL doxycycline
KD = ["SRR33374999", "SRR33374997", "SRR33375000"]        # 75 ng/mL doxycycline

FAMILIES = [
    ("STIM (ER Ca2+ sensor)", ["STIM1", "STIM2"]),
    ("ORAI (CRAC channel)", ["ORAI1", "ORAI2", "ORAI3"]),
    ("SERCA (Ca2+ re-uptake into ER)", ["ATP2A1", "ATP2A2", "ATP2A3"]),
    ("TRPC", ["TRPC1", "TRPC3", "TRPC4", "TRPC5", "TRPC6"]),
    ("SOCE regulators", ["SARAF", "STIMATE", "CRACR2A", "CRACR2B", "CBARP"]),
    ("Mitochondrial Ca2+ uptake", ["MCU", "MICU1", "MICU2", "MICU3", "MCUR1", "MCUB"]),
    ("PMCA (Ca2+ extrusion)", ["ATP2B1", "ATP2B2", "ATP2B3", "ATP2B4"]),
]

tpm = pd.read_csv(f"{V3}/gene_TPM_full.csv", index_col=0)[CTRL + KD]
de = pd.read_csv(f"{V3}/DESeq2_ctrl_vs_75_fullmap.csv", index_col=0)

# ---- per-library median-of-ratios factors on TPM
expr_all = tpm[(tpm > 1).all(axis=1)]
geo = np.exp(np.log(expr_all).mean(axis=1))
sf = expr_all.div(geo, axis=0).median()
adj = tpm.div(sf, axis=1)

# ---- the global shift that the scaling removes
both = (tpm[CTRL].mean(1) > 5) & (tpm[KD].mean(1) > 5)
raw_ratio = float(np.median(tpm.loc[both, KD].mean(1) / tpm.loc[both, CTRL].mean(1)))
adj_ratio = float(np.median(adj.loc[both, KD].mean(1) / adj.loc[both, CTRL].mean(1)))
g_de = de.index.intersection(tpm.index[both])
l2_adj = np.log2(adj.loc[g_de, KD].mean(1) / adj.loc[g_de, CTRL].mean(1))
diff_de = float(np.median(l2_adj - de.loc[g_de, "log2FoldChange"]))
gain = (tpm[KD].mean(1) - tpm[CTRL].mean(1)).sort_values(ascending=False)
scal = pd.DataFrame(
    [{"quantity": f"scaling factor {s}", "value": round(float(sf[s]), 4)} for s in CTRL + KD] +
    [{"quantity": "genes used for the scaling factors (TPM > 1 in all six libraries)", "value": len(expr_all)},
     {"quantity": "genes with mean TPM > 5 in both groups", "value": int(both.sum())},
     {"quantity": "median KD/control TPM ratio, unadjusted", "value": round(raw_ratio, 4)},
     {"quantity": "median KD/control TPM ratio, adjusted", "value": round(adj_ratio, 4)},
     {"quantity": "median (adjusted log2 ratio - DESeq2 log2FC)", "value": round(diff_de, 4)},
     {"quantity": "TPM gained in KD by the 100 largest gainers", "value": round(float(gain[:100].sum()), 0)},
     {"quantity": "CHGA mean TPM, control", "value": round(float(tpm.loc["CHGA", CTRL].mean()), 1)},
     {"quantity": "CHGA mean TPM, knockdown", "value": round(float(tpm.loc["CHGA", KD].mean()), 1)}])

# ---- family table
rows = []
for fam, genes in FAMILIES:
    c_adj = adj.loc[genes, CTRL].mean(1)
    k_adj = adj.loc[genes, KD].mean(1)
    tot_c = c_adj.sum()
    share = 100 * c_adj / tot_c
    dom = share.idxmax() if share.max() > 50 else None
    for g in genes:
        lfc = de.loc[g, "log2FoldChange"] if g in de.index else np.nan
        padj = de.loc[g, "padj"] if g in de.index else np.nan
        rows.append(dict(
            Family=fam, Gene=g,
            TPM_control=round(float(tpm.loc[g, CTRL].mean()), 2),
            TPM_KD=round(float(tpm.loc[g, KD].mean()), 2),
            TPM_control_adjusted=round(float(c_adj[g]), 2),
            TPM_KD_adjusted=round(float(k_adj[g]), 2),
            share_of_family_control_pct=round(float(share[g]), 1),
            dominant=("dominant" if g == dom else ""),
            change_pct_adjusted=(round(100 * (k_adj[g] / c_adj[g] - 1), 1) if c_adj[g] >= 0.05 else np.nan),
            delta_TPM_adjusted=round(float(k_adj[g] - c_adj[g]), 2),
            log2FC=(round(float(lfc), 3) if pd.notna(lfc) else np.nan),
            padj=(float(f"{padj:.3g}") if pd.notna(padj) else np.nan)))
    rows.append(dict(
        Family=fam, Gene="FAMILY TOTAL",
        TPM_control=round(float(tpm.loc[genes, CTRL].mean(1).sum()), 2),
        TPM_KD=round(float(tpm.loc[genes, KD].mean(1).sum()), 2),
        TPM_control_adjusted=round(float(tot_c), 2), TPM_KD_adjusted=round(float(k_adj.sum()), 2),
        share_of_family_control_pct=100.0, dominant="",
        change_pct_adjusted=round(100 * (k_adj.sum() / tot_c - 1), 1),
        delta_TPM_adjusted=round(float(k_adj.sum() - tot_c), 2), log2FC=np.nan, padj=np.nan))
t4 = pd.DataFrame(rows)

for d in (V3, SRC):
    os.makedirs(d, exist_ok=True)
    t4.to_csv(f"{d}/Table4_family_TPM.csv", index=False)
    scal.to_csv(f"{d}/Table4_composition_scaling.csv", index=False)

print("scaling factors:", {k: round(float(v), 3) for k, v in sf.items()})
print(f"median KD/control TPM ratio: unadjusted {raw_ratio:.3f}, adjusted {adj_ratio:.3f}; "
      f"median adjusted-minus-DESeq2 log2 difference {diff_de:+.3f}")
fam = t4[t4.Gene == "FAMILY TOTAL"][["Family", "change_pct_adjusted"]]
print(fam.to_string(index=False))

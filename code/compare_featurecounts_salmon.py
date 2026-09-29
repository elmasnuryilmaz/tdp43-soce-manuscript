#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Concordance of the two DESeq2 runs on the six primary SH-SY5Y libraries (release v1.0.8).

Both runs use GSE296712, three 0 ng/mL control libraries against three 75 ng/mL knockdown libraries
(SRR33374996, SRR33375001, SRR33374995 versus SRR33374999, SRR33374997, SRR33375000):

  * Salmon  : estimated counts summed to genes with the complete GENCODE v47 transcript-to-gene map
              (code/deseq_full.R -> source_data/DESeq2_ctrl_vs_75_fullmap.csv). Every DESeq2 statistic
              quoted in the manuscript comes from this run.
  * featureCounts : gene_name counts against the GENCODE v47 basic annotation, thesis-era run
              (source_data/DESeq2_ctrl_vs_75_featureCounts.csv; copied from the July 2026 reanalysis).

Writes source_data/DESeq2_featureCounts_vs_Salmon_calcium_genes.csv and prints the genome-wide
agreement quoted in Methods 2.2.

Run from the repository root:  /usr/bin/python3 code/compare_featurecounts_salmon.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
sa = pd.read_csv(R / "source_data" / "DESeq2_ctrl_vs_75_fullmap.csv", index_col=0)
fc = pd.read_csv(R / "source_data" / "DESeq2_ctrl_vs_75_featureCounts.csv", index_col=0)

common = fc.index.intersection(sa.index)
a, b = fc.loc[common, "log2FoldChange"], sa.loc[common, "log2FoldChange"]
ok = a.notna() & b.notna()
print(f"genes in featureCounts table: {len(fc)}; in Salmon table: {len(sa)}; counted in both: {int(ok.sum())}")
print(f"log2 fold change, Pearson r = {np.corrcoef(a[ok], b[ok])[0, 1]:.3f}; Spearman r = {a[ok].rank().corr(b[ok].rank()):.3f}")

genes = ["TARDBP", "STMN2", "STIM1", "STIM2", "ORAI1", "ORAI2", "ORAI3", "TRPC1", "SARAF", "STIMATE",
         "CBARP", "ATP2A2", "ATP2A3", "GAPDH"]
rows = []
for g in genes:
    rows.append({"gene": g,
                 "log2FC_featureCounts": fc.loc[g, "log2FoldChange"], "padj_featureCounts": fc.loc[g, "padj"],
                 "log2FC_Salmon": sa.loc[g, "log2FoldChange"], "padj_Salmon": sa.loc[g, "padj"]})
out = pd.DataFrame(rows)
out.to_csv(R / "source_data" / "DESeq2_featureCounts_vs_Salmon_calcium_genes.csv", index=False)
print(out.round(4).to_string(index=False))
print("\nSTMN2 is the one gene of this list on which the two count sets differ markedly; the manuscript quotes no "
      "gene-level DESeq2 value for STMN2 (its cryptic splicing is measured at junction level).")

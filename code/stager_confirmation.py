#!/usr/bin/env python3
"""Read the stage-wise confirmation values quoted for STIM1 in Section 3.4.

Screening stage: DRIMSeq gene-level adjusted p. Confirmation stage: stageR transcript-level
adjusted p. The isoform-level q values and the premature-termination-codon flags come from the
IsoformSwitchAnalyzeR object of the same run.
"""
import pandas as pd

D = ("/Users/elmas/Desktop/TEZ/output/reanalysis_corrected_full_2026-07-22/"
     "isoform_0_vs_75/drimseq_stager")
st = pd.read_csv(f"{D}/stageR_adjusted_pvalues.csv")
s = st[st.geneID == "STIM1"]
print(f"STIM1 transcripts in the stage-wise procedure: {len(s)}")
print(f"  screening-stage q (DRIMSeq gene level): {s.gene.iloc[0]:.3g}")
print(f"  confirmation-stage adjusted p: all {s.transcript.unique()}")
print(f"  transcripts confirmed at 0.05: {(s.transcript < 0.05).sum()}")

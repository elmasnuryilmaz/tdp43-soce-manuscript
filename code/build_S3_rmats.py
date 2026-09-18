#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v3 — S3: rMATS events that meet the manuscript thresholds, six datasets, JC and JCEC.

The complete rMATS output is 23 GB and is deposited separately; this file contains
every event with FDR < 0.05 and |dPSI| >= 0.10 together with the raw junction counts,
so the coverage pre-filter of Section 2.3 can be reproduced from it.
"""
import os
import numpy as np, pandas as pd
from importlib.machinery import SourceFileLoader
core = SourceFileLoader("core", "/Users/elmas/Desktop/MAKALE/04_KOD/01_rmats_core.py").load_module()
SUP = "/Users/elmas/Desktop/MAKALE/09_YAYIN_PAKETI/supplementary"
DS = {"GSE296712_SHSY5Y": "SH-SY5Y (GSE296712)", "GSE230647_iPSC_koloni": "iPSC colonies (GSE230647)",
      "GSE77702_iPSC_MN": "iPSC-derived motor neurons (GSE77702)", "GSE27394_mouse_SE": "Mouse striatum (GSE27394)",
      "Mouse_PE_C2C12": "C2C12 (GSE171714)", "Mouse_PE_NSC34": "NSC34 (GSE171714)"}
out = []
for kind in ("JC", "JCEC"):
    d = core.load_all(kind=kind)
    d = d[(d.FDR < 0.05) & (d.IncLevelDifference.abs() >= 0.10)].copy()
    d["model"] = kind
    out.append(d)
    print(kind, len(d))
A = pd.concat(out, ignore_index=True)
A["dataset"] = A.dataset.map(lambda x: DS.get(x, x))
A["passes_coverage_filter"] = (A.mean_reads_per_sample >= 10) & (A.min_informative_reads >= 5)
cols = ["dataset", "model", "eventType", "geneSymbol", "GeneID", "chr", "strand",
        "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
        "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
        "IncLevel1", "IncLevel2", "IncLevelDifference", "PValue", "FDR",
        "mean_reads_per_sample", "min_informative_reads", "passes_coverage_filter"]
cols = [c for c in cols if c in A.columns]
A = A[cols].rename(columns={"eventType": "event_class", "geneSymbol": "gene",
                            "IncLevelDifference": "delta_PSI",
                            "IncLevel1": "PSI_knockdown_per_replicate",
                            "IncLevel2": "PSI_control_per_replicate"})
p = f"{SUP}/S3_rMATS_significant_events.csv.gz"
A.to_csv(p, index=False, compression="gzip")
print("written:", p, len(A), "rows", round(os.path.getsize(p)/1e6, 1), "MB")
print(A.groupby(["dataset", "model"]).size())

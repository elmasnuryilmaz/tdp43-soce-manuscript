#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Automated consistency check between the reviewed manuscript and the package files.

Every check re-reads the number from the file that produced it and compares it with the
string in the manuscript. Writes 09_YAYIN_PAKETI/logs/consistency_check_reviewed.txt
"""
import io, os, re
import numpy as np
import pandas as pd

M = "/Users/elmas/Desktop/MAKALE"
P = f"{M}/09_YAYIN_PAKETI"
TXT = io.open(f"{P}/manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.md", encoding="utf-8").read()
SEARCH_TXT = re.sub(r"\s+", " ", TXT.replace("*", "").replace("\\", "").replace("--", "–").replace("'", "’")).strip()
out, ok, bad = [], 0, 0


def check(label, snippet, expected=True):
    global ok, bad
    needle = re.sub(r"\s+", " ", snippet.replace("*", "").replace("\\", "").replace("--", "–").replace("'", "’")).strip()
    present = needle in SEARCH_TXT
    good = (present == expected)
    out.append(f"{'PASS' if good else 'FAIL'}  {label}: {snippet[:90]}")
    if good:
        ok += 1
    else:
        bad += 1


def close(label, a, b, tol=0.005):
    global ok, bad
    good = abs(float(a) - float(b)) <= tol
    out.append(f"{'PASS' if good else 'FAIL'}  {label}: manuscript {a} vs data {b}")
    if good:
        ok += 1
    else:
        bad += 1


out.append("=== Table 2: coverage pre-filter ===")
t1 = pd.read_csv(f"{P}/tables/Table2_coverage_prefilter.csv")
close("events removed range low", 18, round(t1.events_removed_pct.min()), tol=0.6)
close("events removed range high", 78, round(t1.events_removed_pct.max()), tol=0.6)
close("significant lost low", 33, round(t1.significant_lost_pct.min()), tol=0.6)
close("significant lost high", 76, round(t1.significant_lost_pct.max()), tol=0.6)
close("SH-SY5Y significant before filter", 7854,
      int(t1.loc[t1.dataset.str.startswith("SH-SY5Y"), "significant_before_filter"].iloc[0]), tol=0)

out.append("\n=== Table 3: robust SOCE events ===")
t2 = pd.read_csv(f"{P}/tables/Table3_robust_SOCE_splicing_events.csv").set_index("gene")
for g, d, lo, hi in [("STIMATE", 0.244, 0.095, 0.368), ("ORAI3", -0.269, -0.404, -0.107),
                     ("STIM2", -0.120, -0.165, -0.064), ("STIM1", 0.145, -0.002, 0.293)]:
    close(f"{g} delta_PSI", d, t2.loc[g, "delta_PSI"])
    close(f"{g} CI low", lo, t2.loc[g, "CI95_low"])
    close(f"{g} CI high", hi, t2.loc[g, "CI95_high"])
close("STIM1 exon start 1-based", 4088702, int(t2.loc["STIM1", "start_1based"]), tol=0)
close("STIM1 exon length", 37, int(t2.loc["STIM1", "exon_bp"]), tol=0)

out.append("\n=== Table 4 / S14: cryptic calls and null test ===")
t3 = pd.read_csv(f"{P}/tables/Table4_cryptic_events_eleven_comparisons.csv").set_index("comparison")
close("SH-SY5Y stringent-filter events (regtools)", 165,
      int(t3.loc["SH-SY5Y 75 ng/mL", "stringent_filter_events"]), tol=0)
close("SH-SY5Y genes", 113, int(t3.loc["SH-SY5Y 75 ng/mL", "genes"]), tol=0)
close("iPSC-MN TDP-43 stringent-filter events", 18, int(t3.loc["iPSC-MN, TDP-43 KD", "stringent_filter_events"]), tol=0)
close("iPSC-MN FUS stringent-filter events", 26, int(t3.loc["iPSC-MN, FUS KD", "stringent_filter_events"]), tol=0)
close("iPSC-MN TAF15 stringent-filter events", 23, int(t3.loc["iPSC-MN, TAF15 KD", "stringent_filter_events"]), tol=0)
close("iPSC-MN TDP-43 permissive positive controls", 2,
      int(t3.loc["iPSC-MN, TDP-43 KD", "positive_controls_permissive"]), tol=0)
close("iPSC-MN FUS permissive positive controls", 0,
      int(t3.loc["iPSC-MN, FUS KD", "positive_controls_permissive"]), tol=0)
close("SH-SY5Y permissive positive controls", 13,
      int(t3.loc["SH-SY5Y 75 ng/mL", "positive_controls_permissive"]), tol=0)
close("iPSC colonies permissive positive controls", 15,
      int(t3.loc["iPSC colonies", "positive_controls_permissive"]), tol=0)
close("permissive call counts, TDP-43 / FUS / TAF15", 141,
      int(t3.loc["iPSC-MN, TDP-43 KD", "permissive_events"]), tol=0)
close("burden minimum", 12, int(t3.stringent_filter_events.min()), tol=0)
_mouse = t3.loc[["C2C12", "NSC34", "Mouse striatum"]]
close("Table 4 mouse rows: positive controls not assessed", 3,
      int((_mouse.positive_controls_stringent_filter == "n/a (mouse)").sum()), tol=0)
close("burden maximum", 477, int(t3.stringent_filter_events.max()), tol=0)
s14 = pd.read_csv(f"{P}/supplementary/S14_control_vs_control_null_test.csv")
_s14p = s14[s14.tier == "permissive"].set_index("dataset")
close("permissive null ratio iPSC colonies (Methods 2.5: 0.98 per real call)", 0.98,
      _s14p.loc["iPSC colonies", "null_to_real_ratio"])
close("permissive null ratio K562 total RNA (Methods 2.5: 2.01)", 2.01,
      _s14p.loc["K562 total RNA", "null_to_real_ratio"])
s14 = s14[s14.tier == "tier 2 (stringent filter)"].set_index("dataset")
close("null ratio iPSC colonies", 0.64, s14.loc["iPSC colonies", "null_to_real_ratio"])
close("null ratio K562 total RNA", 2.17, s14.loc["K562 total RNA", "null_to_real_ratio"])
close("null ratio mouse striatum", 0.83, s14.loc["Mouse striatum", "null_to_real_ratio"])

out.append("\n=== Table 1: transcript-family abundance (TPM) ===")
t4 = pd.read_csv(f"{P}/tables/Table1_transcript_family_abundance.csv")
g = t4[t4.Gene != "FAMILY TOTAL"].set_index("Gene")
close("ORAI2 share", 77, round(g.loc["ORAI2", "share_of_family_control_pct"]), tol=0.6)
close("ATP2A2 share", 98, round(g.loc["ATP2A2", "share_of_family_control_pct"]), tol=0.6)
close("ORAI1 share", 15, round(g.loc["ORAI1", "share_of_family_control_pct"]), tol=0.6)
close("SARAF share", 82, round(g.loc["SARAF", "share_of_family_control_pct"]), tol=0.6)
close("CBARP log2FC", -1.254, g.loc["CBARP", "log2FC"])
close("STIM1 log2FC", 0.929, g.loc["STIM1", "log2FC"])
close("ORAI1 log2FC", 0.433, g.loc["ORAI1", "log2FC"])
close("ATP2A3 log2FC", 1.306, g.loc["ATP2A3", "log2FC"])
close("TRPC1 log2FC", 0.958, g.loc["TRPC1", "log2FC"])
close("STIM1 share", 64, round(g.loc["STIM1", "share_of_family_control_pct"]), tol=0.6)
close("TRPC1 share", 98, round(g.loc["TRPC1", "share_of_family_control_pct"]), tol=0.6)
close("ATP2A3 share", 1.5, g.loc["ATP2A3", "share_of_family_control_pct"], tol=0.05)
close("ORAI2 log2FC (-0.17)", -0.17, g.loc["ORAI2", "log2FC"], tol=0.005)
close("ATP2A2 log2FC (-0.25)", -0.25, g.loc["ATP2A2", "log2FC"], tol=0.005)
close("ORAI3 log2FC (2.06)", 2.06, g.loc["ORAI3", "log2FC"], tol=0.005)
close("SARAF log2FC (0.55)", 0.55, g.loc["SARAF", "log2FC"], tol=0.005)
fam = t4[t4.Gene == "FAMILY TOTAL"].set_index("Family")
# composition-adjusted TPM (build_family_abundance.py); unadjusted TPM understates every
# change by about a quarter because a few abundant transcripts gain share in knockdown
for f, v in [("SERCA (Ca2+ re-uptake into ER)", -12), ("Mitochondrial Ca2+ uptake", -25),
             ("STIM (ER Ca2+ sensor)", 48), ("ORAI (CRAC channel)", 27), ("Ca2+-entry regulators", 30)]:
    close(f"family net, adjusted, {f}", v, round(fam.loc[f, "change_pct_adjusted"]), tol=0.6)
_o = t4[t4.Family == "ORAI (CRAC channel)"].set_index("Gene")
_okd = _o.loc["FAMILY TOTAL", "TPM_KD_adjusted"]
close("ORAI3 share of ORAI pool, control (8%)", 8, round(_o.loc["ORAI3", "share_of_family_control_pct"]), tol=0.6)
close("ORAI3 share of ORAI pool, knockdown (30%)", 30, round(100 * _o.loc["ORAI3", "TPM_KD_adjusted"] / _okd), tol=0.6)
close("ORAI2 share of ORAI pool, knockdown (53%)", 53, round(100 * _o.loc["ORAI2", "TPM_KD_adjusted"] / _okd), tol=0.6)
_sc = pd.read_csv(f"{P}/source_data/Table4_composition_scaling.csv").set_index("quantity").value
close("median KD/control TPM ratio, unadjusted (24% lower)", 24,
      round(100 * (1 - _sc["median KD/control TPM ratio, unadjusted"])), tol=0.6)
close("genes used for the scaling factors", 13907,
      _sc["genes used for the scaling factors (TPM > 1 in all six libraries)"], tol=0)
close("median adjusted-minus-DESeq2 log2 difference (-0.02)", -0.02,
      _sc["median (adjusted log2 ratio - DESeq2 log2FC)"], tol=0.005)
close("CHGA control TPM (1,857)", 1857, _sc["CHGA mean TPM, control"], tol=0.6)
close("CHGA knockdown TPM (10,675)", 10675, _sc["CHGA mean TPM, knockdown"], tol=0.6)

out.append("\n=== Table 5 / S16: patient tissue ===")
t5 = pd.read_csv(f"{P}/tables/Table5_cross_disease_comparison.csv")
als = t5[(t5.group == "ALS") & (t5.gene == "TRPC1")]
brain = als[~als.region.str.startswith("Spinal")]
sig = brain[brain.q_value < 0.05]
close("TRPC1 significant brain regions", 6, len(sig), tol=0)
close("TRPC1 brain regions tested", 7, len(brain), tol=0)
close("TRPC1 delta min (significant)", 0.447, sig.cliffs_delta.min())
close("TRPC1 delta max (significant)", 0.699, sig.cliffs_delta.max())
cord = als[als.region.str.startswith("Spinal")]
close("spinal cord levels", 3, len(cord), tol=0)
close("thoracic cord delta", -0.321, float(cord.loc[cord.region.str.contains("thoracic"), "cliffs_delta"].iloc[0]))
close("significant cord levels", 0, int((cord.q_value < 0.05).sum()), tol=0)
a3 = t5[(t5.group == "ALS") & (t5.gene == "TRPC1")]
don = pd.read_csv(f"{P}/supplementary/S16b_multiple_sclerosis_donor_level.csv")
nawm = don[(don.comparison == "NAWM vs control WM") & (don.gene == "TRPC1") & (don.unit == "donor")]
close("MS NAWM donor-level delta", -0.771, float(nawm.cliffs_delta.iloc[0]))
close("MS NAWM donor-level p", 0.030, float(nawm.p.iloc[0]), tol=0.001)
adj = don[(don.comparison.str.contains("myelin")) & (don.unit == "donor")]
close("MS myelin-adjusted donor-level delta", -0.640, float(adj.cliffs_delta.iloc[0]))
close("MS myelin-adjusted donor-level p", 0.055, float(adj.p.iloc[0]), tol=0.001)

out.append("\n=== S8 / S9: cryptic STMN2 in ALS tissue ===")
s8 = pd.read_csv(f"{P}/supplementary/S8_cryptic_STMN2_ALS_vs_control.csv")
s8 = s8[s8.gene == "STMN2"].set_index("region")
close("lumbar cord delta", 0.705, s8.loc["Spinal Cord Lumbar", "cliffs_delta"])
close("cervical cord delta", 0.619, s8.loc["Spinal Cord Cervical", "cliffs_delta"])
close("medial motor cortex delta", 0.337, s8.loc["Cortex Motor Medial", "cliffs_delta"])
close("temporal cortex delta", 0.283, s8.loc["Cortex Temporal", "cliffs_delta"])
close("cerebellum delta", -0.040, s8.loc["Cerebellum", "cliffs_delta"])
s9 = pd.read_csv(f"{P}/supplementary/S9_cryptic_PSI_correlations_within_ALS.csv")
k = s9[s9.proxy == "cryptic STMN2 PSI"]
close("correlations against junction marker", 110, len(k), tol=0)
close("correlations in total", 231, len(s9), tol=0)
close("surviving correction", 8, int((k.q_value < 0.05).sum()), tol=0)
tr = k[k.target_gene == "TRPC1"]
close("TRPC1 rho min", -0.19, tr.spearman_rho.min(), tol=0.006)
close("TRPC1 rho max", 0.18, tr.spearman_rho.max(), tol=0.006)
close("TRPC1 significant correlations", 0, int((tr.q_value < 0.05).sum()), tol=0)

out.append("\n=== S11: corrected APA, both models ===")
s11 = pd.read_csv(f"{P}/supplementary/S11_APA_candidate_gradients.csv")
sh = s11[s11.dataset.str.startswith("SH-SY5Y")].set_index(["gene", "unit"])
close("SH-SY5Y STIM1 intron17 delta", -0.173, sh.loc[("STIM1", "intron17"), "delta"])
close("SH-SY5Y STIM2 intron13 delta", -0.156, sh.loc[("STIM2", "intron13"), "delta"])
close("SH-SY5Y STIM2 intron13 control index (0.96)", 0.96, sh.loc[("STIM2", "intron13"), "index_control"], tol=0.005)
_w = [int(v) for v in sh.loc[("STIM2", "intron13"), "window_5prime_or_proximal"].split("-")]
_t2 = pd.read_csv(f"{P}/tables/Table3_robust_SOCE_splicing_events.csv").set_index("gene")
close("STIM2 intron13 5' window contains the Table 3 STIM2 exon", 1,
      int(_w[0] <= int(_t2.loc["STIM2", "start_1based"]) and int(_t2.loc["STIM2", "end"]) <= _w[1]), tol=0)
close("S11: every unit has a genomic window", 0, int(s11.window_5prime_or_proximal.isna().sum()), tol=0)
mn = s11[s11.dataset.str.startswith("iPSC-MN")]
close("iPSC-MN units passing the depth filter", 59, len(mn), tol=0)
close("iPSC-MN genes", 29, mn.gene.nunique(), tol=0)
mni = mn.set_index(["gene", "unit"])
close("iPSC-MN STMN2 intron2 delta", 0.249, mni.loc[("STMN2", "intron2"), "delta"], tol=0.001)
close("iPSC-MN STMN2 intron2 index KD", 0.819, mni.loc[("STMN2", "intron2"), "index_knockdown"], tol=0.001)
close("iPSC-MN STMN2 intron2 index control", 0.570, mni.loc[("STMN2", "intron2"), "index_control"], tol=0.001)
close("iPSC-MN ATP2A2 intron3 delta", -0.164, mni.loc[("ATP2A2", "intron3"), "delta"], tol=0.001)
close("iPSC-MN SARAF intron5 delta", 0.097, mni.loc[("SARAF", "intron5"), "delta"], tol=0.001)
close("iPSC-MN TRPC1 intron1 delta", 0.083, mni.loc[("TRPC1", "intron1"), "delta"], tol=0.001)
close("iPSC-MN STIM2 terminal exon delta", -0.074, mni.loc[("STIM2", "termexon"), "delta"], tol=0.001)

c2 = s11[s11.dataset == "C2C12"]; ns = s11[s11.dataset == "NSC34"]
close("C2C12 units", 74, len(c2), tol=0)
close("NSC34 units", 131, len(ns), tol=0)
close("NSC34 Saraf intron5 delta", 0.204,
      ns.set_index(["gene", "unit"]).loc[("Saraf", "intron5"), "delta"], tol=0.001)
close("C2C12 Atp2a2 intron6 delta", 0.285,
      c2.set_index(["gene", "unit"]).loc[("Atp2a2", "intron6"), "delta"], tol=0.001)
close("C2C12 Trpc1 intron7 delta", 0.260,
      c2.set_index(["gene", "unit"]).loc[("Trpc1", "intron7"), "delta"], tol=0.001)

out.append("\n=== isoform-level tests quoted for CBARP ===")
_iso = pd.read_csv("/Users/elmas/Desktop/TEZ/output/reanalysis_corrected_full_2026-07-22/tables/"
                   "corrected_isoform_nmd_summary.csv").set_index("gene")
close("CBARP IsoformSwitchAnalyzeR gene-level q (0.12)", 0.12, _iso.loc["CBARP", "gene_switch_q"], tol=0.005)
close("CBARP DRIMSeq q (0.34)", 0.34, _iso.loc["CBARP", "DRIMSeq_q"], tol=0.005)
close("CBARP significant switching isoforms", 0, int(_iso.loc["CBARP", "n_significant_switch_isoforms"]), tol=0)
_sw = pd.read_csv("/Users/elmas/Desktop/TEZ/output/reanalysis_corrected_full_2026-07-22/isoform_0_vs_75/"
                  "top_switching_genes.csv").set_index("gene_name")
close("CBARP-DT IsoformSwitchAnalyzeR gene-level q (0.008)", 0.008, _sw.loc["CBARP-DT", "gene_switch_q_value"], tol=0.0005)

out.append("\n=== S1: laboratory source data ===")
x = pd.read_excel(f"{P}/supplementary/S1_laboratory_source_data.xlsx", sheet_name="Summary_stats")
x = x.set_index(x.measurement + " | " + x.group)
close("SOCE control mean", 1.542, x.loc["SOCE delta F340/F380 | Non-targeting shRNA control", "mean"])
close("SOCE knockdown mean", 0.245, x.loc["SOCE delta F340/F380 | shTDP-43", "mean"])
close("ER release control mean", 0.268, x.loc["ER Ca2+ release delta F340/F380 | Non-targeting shRNA control", "mean"])
close("WST-1 48 h knockdown", 61.5, x.loc["WST-1 signal 48 h | shTDP-43", "mean"], tol=0.05)

import subprocess as _sp
close("manuscript and DATA_AVAILABILITY name the same release", 1,
      int(len(set(re.findall(r"releases/tag/(v\d+\.\d+\.\d+)",
                            TXT + io.open(f"{P}/DATA_AVAILABILITY.md", encoding="utf-8").read()))) == 1), tol=0)
out.append("\n=== comparison group and cell composition (Tables S18, S18b, S18c) ===")
import gzip as _gz
_B = "/Users/elmas/Desktop/TEZ/output/ek_analizler_2026-07-31"
_meta = pd.read_csv(f"{_B}/06_WGCNA_ALS_kohort/T60_NYGC_ornek_ustverisi.csv", low_memory=False)
with _gz.open(f"{_B}/00_ham_veri_onbellek/gse153960/GSE153960_counts.txt.gz", "rt") as _fh:
    _cols = {c for c in _fh.readline().rstrip("\n").split("\t") if c.startswith("CGND")}
_meta = _meta.dropna(subset=["ornek_id"]).drop_duplicates("ornek_id")
_meta = _meta[_meta.ornek_id.isin(_cols)]
_g = _meta.grup.astype(str).str.strip().str.replace("DIsorders", "Disorders").value_counts()
close("NYGC samples with metadata in the count matrix", 1640, len(_meta), tol=0)
close("samples labelled both ALS and comparison group", 266, int(_g.get("ALS Spectrum MND, Other Neurological Disorders", 0)), tol=0)
close("samples labelled both Pre-fALS and comparison group", 2, int(_g.get("Pre-fALS, Other Neurological Disorders", 0)), tol=0)

_s18 = pd.read_csv(f"{P}/supplementary/S18_TRPC1_cell_composition_adjustment.csv")
_t = _s18[_s18.variable == "TRPC1"]
_six = ["Cerebellum", "Cortex Frontal", "Cortex Motor Lateral", "Cortex Motor Medial", "Cortex Temporal", "Hippocampus"]
_a = _t[(_t.group == "ALS") & _t.region.isin(_six) & (_t.adjustment != "none")]
close("ALS adjusted TRPC1 delta, lowest", 0.26, _a.cliffs_delta.min(), tol=0.005)
close("ALS adjusted TRPC1 delta, highest", 0.66, _a.cliffs_delta.max(), tol=0.005)
close("ALS adjusted TRPC1 delta positive in every region and model", 1, int((_a.cliffs_delta > 0).all()), tol=0)
_rob = sorted(r for r in _six if (_a[_a.region == r].q_value < 0.05).all())
close("significant under every model: cerebellum, frontal, medial motor", 1,
      int(_rob == ["Cerebellum", "Cortex Frontal", "Cortex Motor Medial"]), tol=0)
close("significant in all six after RBFOX3", 6, int((_a[_a.adjustment == "RBFOX3"].q_value < 0.05).sum()), tol=0)
_o = _t[_t.group == "Other neurological disorders"].set_index(["region", "adjustment"]).cliffs_delta
close("comparison group, temporal TRPC1 before adjustment", -0.571, _o[("Cortex Temporal", "none")], tol=0.0005)
close("comparison group, temporal TRPC1 after SNAP25", -0.169, _o[("Cortex Temporal", "SNAP25")], tol=0.0005)
close("comparison group, temporal TRPC1 after RBFOX3", -0.390, _o[("Cortex Temporal", "RBFOX3")], tol=0.0005)
_fr = _o.loc["Cortex Frontal"].drop("none"); _ce = _o.loc["Cerebellum"].drop("none")
close("comparison group, frontal adjusted range low", -0.68, _fr.min(), tol=0.005)
close("comparison group, frontal adjusted range high", -0.42, _fr.max(), tol=0.005)
close("comparison group, cerebellum adjusted range low", -0.42, _ce.min(), tol=0.005)
close("comparison group, cerebellum adjusted range high", -0.39, _ce.max(), tol=0.005)
_mk = _s18[_s18.group == "Other neurological disorders"].set_index(["variable", "region"]).cliffs_delta
for _v, _r, _x in (("SNAP25", "Cortex Frontal", -0.57), ("SNAP25", "Cortex Temporal", -0.55),
                   ("GFAP", "Cortex Frontal", 0.59), ("GFAP", "Cortex Temporal", 0.64)):
    close(f"comparison group {_v} delta, {_r}", _x, _mk[(_v, _r)], tol=0.005)

_b = pd.read_csv(f"{P}/supplementary/S18b_cryptic_STMN2_by_group_and_region.csv").set_index(["group", "region"])
_OND, _ALS, _CT = "Other neurological disorders", "ALS", "Non-neurological control"
for _grp, _r, _k, _n in ((_OND, "Cortex Frontal", 21, 42), (_OND, "Cortex Temporal", 22, 35),
                         (_ALS, "Cortex Frontal", 2, 154), (_ALS, "Cortex Temporal", 2, 25),
                         (_CT, "Cortex Frontal", 0, 55), (_CT, "Cortex Temporal", 0, 24),
                         (_OND, "Cerebellum", 0, 49), (_ALS, "Cerebellum", 0, 157)):
    close(f"S18b {_grp[:12]} {_r}: samples above 1%", _k, int(_b.loc[(_grp, _r), "n_PSI_above_0_01"]), tol=0)
    close(f"S18b {_grp[:12]} {_r}: samples", _n, int(_b.loc[(_grp, _r), "n_samples"]), tol=0)
close("comparison group vs control delta, frontal", 0.73, _b.loc[(_OND, "Cortex Frontal"), "cliffs_delta_vs_control"], tol=0.005)
close("comparison group vs control delta, temporal", 0.71, _b.loc[(_OND, "Cortex Temporal"), "cliffs_delta_vs_control"], tol=0.005)
close("comparison group vs control q below 1e-6 in cortex", 1,
      int(max(_b.loc[(_OND, "Cortex Frontal"), "q_value"], _b.loc[(_OND, "Cortex Temporal"), "q_value"]) < 1e-6), tol=0)
_bb = _b.reset_index()
_alsb = _bb[(_bb.group == _ALS) & ~_bb.region.str.startswith("Spinal")]
_alsc = _bb[(_bb.group == _ALS) & _bb.region.str.startswith("Spinal")]
_ondc = _bb[(_bb.group == _OND) & _bb.region.str.startswith("Cortex")]
close("ALS brain: at most 14% above 1%", 1, int(_alsb.fraction_PSI_above_0_01.max() <= 0.14), tol=0)
close("ALS cord: lowest fraction 41%", 0.41, _alsc.fraction_PSI_above_0_01.min(), tol=0.005)
close("ALS cord: highest fraction 67%", 0.67, _alsc.fraction_PSI_above_0_01.max(), tol=0.005)
close("comparison-group cortex: lowest fraction 50%", 0.50, _ondc.fraction_PSI_above_0_01.min(), tol=0.005)
close("comparison-group cortex: highest fraction 63%", 0.63, _ondc.fraction_PSI_above_0_01.max(), tol=0.005)

_c = pd.read_csv(f"{P}/supplementary/S18c_cryptic_STMN2_within_comparison_group.csv")
_w = _c[_c.family.str.startswith("cryptic")].set_index(["region", "y"])
for _r, _y, _x in (("Cortex Frontal", "TRPC1", -0.50), ("Cortex Temporal", "TRPC1", -0.40),
                   ("Cortex Frontal", "SNAP25", -0.49), ("Cortex Temporal", "SNAP25", -0.47)):
    close(f"within comparison group rho {_y} {_r}", _x, _w.loc[(_r, _y), "spearman_rho"], tol=0.005)
close("partial rho TRPC1 frontal", -0.24, _w.loc[("Cortex Frontal", "TRPC1"), "partial_rho_given_SNAP25"], tol=0.005)
close("partial p TRPC1 frontal", 0.12, _w.loc[("Cortex Frontal", "TRPC1"), "partial_p_value"], tol=0.005)
close("partial rho TRPC1 temporal", 0.00, _w.loc[("Cortex Temporal", "TRPC1"), "partial_rho_given_SNAP25"], tol=0.005)
close("partial p TRPC1 temporal", 0.98, _w.loc[("Cortex Temporal", "TRPC1"), "partial_p_value"], tol=0.005)
_k = _c[_c.family.str.startswith("TRPC1 with")]
close("TRPC1-SNAP25 rho, lowest", 0.34, _k.spearman_rho.min(), tol=0.005)
close("TRPC1-SNAP25 rho, highest", 0.86, _k.spearman_rho.max(), tol=0.005)
close("TRPC1-SNAP25 rho positive everywhere", 1, int((_k.spearman_rho > 0).all()), tol=0)

for _i in ("S18", "S18b", "S18c"):
    close(f"Supplementary Table {_i} cited in the text", 1,
          int(re.search(rf"Table {_i}(?![0-9a-z])", SEARCH_TXT.split("## References")[0]) is not None), tol=0)
for _s in ["1,640 samples with metadata after filtering",
           "Groups were used as single labels: 266 samples annotated with both ALS Spectrum MND and Other Neurological Disorders",
           "the adjusted δ remained positive under every marker combination (+0.26 to +0.66)",
           "Adjustment for SNAP25 removed most of the temporal-cortex decrease",
           "(δ = −0.571 before and −0.169 after; −0.390 after adjustment for RBFOX3)",
           "21 of 42 frontal and 22 of 35 temporal cortex samples, against 2 of 154 and 2 of 25 ALS samples",
           "whereas in cerebellum no comparison-group or ALS sample exceeded 1%",
           "although the diagnoses cannot be checked",
           "(partial ρ = −0.24, p = 0.12, and 0.00, p = 0.98",
           "exceeded 1% of reads in at most 14% of samples, showed no statistically significant difference in ALS spinal cord (41–67%)",
           "fell in the cortex of the comparison group (50–63%)",
           "the direction it shares with the cellular model does not by itself link the two",
           "which regression on marker genes adjusts for only partially",
           "the cryptic STMN2 junction indicates TDP-43 loss of function rather than a diagnosis"]:
    check("present", _s, True)
for _s in ["1,641 samples", "make a simple neurodegeneration explanation less compelling",
           "these observations identify disease-associated candidates rather than a direct TDP-43-driven mechanism"]:
    check("absent", _s, False)

out.append("\n=== stage-wise confirmation for STIM1 ===")
_st = pd.read_csv("/Users/elmas/Desktop/TEZ/output/reanalysis_corrected_full_2026-07-22/"
                  "isoform_0_vs_75/drimseq_stager/stageR_adjusted_pvalues.csv")
_st1 = _st[_st.geneID == "STIM1"]
close("STIM1 transcripts in the stage-wise procedure", 7, len(_st1), tol=0)
close("STIM1 screening-stage q", 5.49e-7, float(_st1.gene.iloc[0]), tol=1e-8)
close("STIM1 transcripts confirmed at 0.05", 0, int((_st1.transcript < 0.05).sum()), tol=0)
close("every evaluable STIM1 transcript has a confirmation-stage p of 1.0", 1,
      int((_st1.transcript == 1.0).all()), tol=0)
for _s in ["none of seven evaluable transcripts passed stageR confirmation",
           "all confirmation-stage adjusted p values = 1.0",
           "the two predicted PTC isoforms did not change (q = 0.33 and 0.77)",
           "stageR supplies transcript-level confirmation"]:
    check("present", _s, True)
check("absent", "Isoform-level testing supported STIM1", False)

out.append("\n=== svaseq sensitivity, recomputed 28 September 2026 ===")
_sva = pd.read_csv(f"{P}/source_data/svaseq_sensitivity_SHSY5Y.csv")
close("genes with an adjusted p value in both fits", 14012, len(_sva), tol=0)
_de0 = _sva[(_sva.padj_published < 0.05) & (_sva.log2FC_published.abs() >= 1)]
_de1 = _sva[(_sva.padj_with_SV < 0.05) & (_sva.log2FC_with_SV.abs() >= 1)]
_shared = set(_de0.gene) & set(_de1.gene)
close("differentially expressed in the published model", 1694, len(_de0), tol=0)
close("differentially expressed with surrogate variables", 1067, len(_de1), tol=0)
close("shared between the two fits", 850, len(_shared), tol=0)
close("median change in log2 fold change", 0.20,
      (_sva.log2FC_published - _sva.log2FC_with_SV).abs().median(), tol=0.005)
_sh = _sva[_sva.gene.isin(_shared)]
close("direction agrees for every shared gene", 1,
      int((np.sign(_sh.log2FC_published) == np.sign(_sh.log2FC_with_SV)).all()), tol=0)
_g = _sva.set_index("gene")
for _n in ["STIM1", "TRPC1", "ORAI3", "SARAF", "CBARP"]:
    close(f"{_n} retains q < 0.05 in both fits (not necessarily the FC threshold)", 1,
          int(_g.loc[_n, "padj_with_SV"] < 0.05 and _g.loc[_n, "padj_published"] < 0.05), tol=0)
close("only ORAI3 and CBARP meet both DE cutoffs in both fits among these five", 1,
      int(set([n for n in ["STIM1", "TRPC1", "ORAI3", "SARAF", "CBARP"] if n in _shared]) == {"ORAI3", "CBARP"}), tol=0)
for _n in ["ORAI1", "ATP2A3", "ATP2A2", "STIM2"]:
    close(f"{_n} does not stay significant with surrogate variables", 0,
          int(_g.loc[_n, "padj_with_SV"] < 0.05), tol=0)
for _s in ["which estimated two surrogate variables in the primary SH-SY5Y comparison",
           "reduced the number meeting them from 1,694 to 1,067",
           "kept the direction of all 850 genes",
           "source_data/svaseq_sensitivity_SHSY5Y.csv"]:
    check("present", _s, True)
check("absent", "the principal findings were unchanged", False)

out.append("\n=== analyses moved to the supplementary material ===")
_supp2 = _sp.run(["pandoc", "-t", "plain", "--wrap=none",
                  f"{P}/supplementary/SUPPLEMENTARY_MATERIAL.docx"], capture_output=True, text=True).stdout
for _v in ["The median absolute difference in log2 fold change was 0.20", "14,012 genes"]:
    close(f"svaseq statistic kept in Supplementary Results 5: {_v[:40]}", 1, int(_v in _supp2), tol=0)
for _n, _title in ((1, "Aberrant splicing outlier detection (FRASER)"), (2, "Nonsense-mediated decay interaction"),
                   (3, "Coverage-based polyadenylation screen: details"), (4, "Multiple sclerosis in detail"),
                   (5, "Surrogate-variable sensitivity analysis (svaseq)"), (6, "STIM1 isoform-level testing")):
    close(f"Supplementary Results {_n} present in the supplementary document", 1,
          int(f"Supplementary Results {_n}. {_title}" in _supp2), tol=0)
    close(f"Supplementary Results {_n} cited in the manuscript", 1,
          int(f"Supplementary Results {_n}" in SEARCH_TXT), tol=0)

for _v in ["human SARAF locus in iPSC-derived motor neurons (+0.097)",
           "human SH-SY5Y estimate was 0.000 (interval −0.264 to +0.241)",
           "(δ = −0.482, p = 0.005), although not after Benjamini–Hochberg correction (q = 0.10)",
           "from δ = −0.562 to Cliff’s δ of −0.418 on the residuals at sample level (p = 0.002)"]:
    close(f"moved text kept in the supplementary: {_v[:40]}", 1, int(_v in _supp2), tol=0)
for _v in ["without t tests, sign tests, panel-enrichment tests or FDR claims", "median interaction of −0.11 log₂", "not validated NMD rescue in this experiment",
           "4.8 versus 2.3 events", "Atp2a2 intron 6 (+0.285)", "Trpc1 intron 7 (+0.260)",
           "δ = −0.562 to Cliff’s δ of −0.418", "donor-level δ = −1.000 for both",
           "range +1.19 to +2.23", "STIM1 had a mean interaction of −0.325 log₂"]:
    close(f"value kept in the supplementary text: {_v[:40]}", 1, int(_v in _supp2), tol=0)
for _s in ["FRASER, run on the nine-sample doxycycline series, returned no genome-wide",
           "Because the four interventions reuse the same control and knockdown libraries",
           "RYR2 and MCU also had positive interactions",
           "the myelin markers MBP (δ = +0.051), PLP1 (−0.074), MOG (−0.257) and MAG (−0.299)"]:
    check("absent", _s, False)

out.append("\n=== 27 September 2026 structural scan ===")
for _s in ["Two of the 24 dataset–panel combinations had nominal permutation p < 0.05",
           "Neither survived adjustment for the number of comparisons",
           "3.4 Splicing changes in core SOCE-pathway genes that pass the robustness checks",
           "alternative polyadenylation, which generates the truncated STMN2 transcript",
           "coupling to nonsense-mediated decay",
           "Two of the 24 dataset–panel combinations had nominal permutation p < 0.05 after matching"]:
    check("present", _s, True)
for _s in ["This is the central observation of the study",
           "tentative, motor-neuron-associated observation",
           "Robust splicing changes concentrate in SOCE regulators",
           "they removed isolated events and apparent panel enrichment"]:
    check("absent", _s, False)

out.append("\n=== 26 September 2026 referee round ===")
from scipy import stats as _st2
close("smallest two-sided rank-test p at 3 vs 3", 0.10,
      _st2.mannwhitneyu([3, 2, 1], [6, 5, 4], alternative="two-sided").pvalue, tol=0.001)
_ga = io.open(f"{P}/figures/graphical_abstract.svg", encoding="utf-8").read()
close("graphical abstract names the one-plate design", 1, int("one plate: 3 wells/group" in _ga), tol=0)
close("graphical abstract omits the well-level p value", 0, int("p = 0.035" in _ga), tol=0)
close("graphical abstract avoids an established SOCE decrease claim", 0,
      int("reduces store-operated Ca" in _ga), tol=0)
close("graphical abstract names the observed readdition response", 1,
      int("-readdition observation" in _ga), tol=0)
close("graphical abstract no longer quotes the Student p value", 0, int("0.0115" in _ga), tol=0)
for _s in ["prepared for measurement 72 h after transduction, while puromycin selection was still in progress",
           "onto disinfected glass coverslips in 24-well plates one day before the measurement",
           "measured while still attached to the coverslip, which was mounted in the cuvette",
           "three wells on the same plate (n = 3 technical replicates)",
           "Amplification efficiencies were not determined",
           "was not re-measured in the June and July set",
           "GSE307054, the dataset of a preprint",
           "SARAF and CBARP changed in the same direction in cervical and lumbar cord",
           "cannot distinguish altered ER store depletion from changes in STIM activation",
           "neither the cell number nor the dye loading of each cuvette was recorded",
           "BTP2, Synta66 or Gd³⁺",
           "the novel-site junction at this locus carried six",
           "The original recordings, with the CPA and CaCl₂ additions, are shown at their own axis ranges in Supplementary Figure S9",
           "no value is smoothed or rescaled",
           "which requires independent biological replication",
           "it does not estimate between-experiment variability",
           "WST-1 signal at 48 h, four wells from one experiment",
           "The APA and NMD analyses are hypothesis-generating"]:
    check("present", _s, True)
close("no independent-culture Fura-2 claims remain", 0,
      len(re.findall(r"three independent cultures", SEARCH_TXT)), tol=0)

out.append("\n=== 25 September 2026 audit corrections ===")
_s9 = pd.read_csv(f"{P}/supplementary/S9_cryptic_PSI_correlations_within_ALS.csv")
_s9 = _s9[(_s9.proxy == "cryptic STMN2 PSI") & _s9.target_gene.isin(["TRPC1", "SARAF", "CBARP"])
          & (_s9.in_correction_family == "yes")]
close("proxy correlations, number of tests", 30, len(_s9), tol=0)
close("proxy correlations, lowest rho", -0.29, _s9.spearman_rho.min(), tol=0.005)
close("proxy correlations, highest rho", 0.18, _s9.spearman_rho.max(), tol=0.005)
close("proxy correlations passing q < 0.05", 1, int((_s9.q_value < 0.05).sum()), tol=0)
_sig = _s9[_s9.q_value < 0.05].iloc[0]
close("the single significant proxy correlation is negative", 1, int(_sig.spearman_rho < 0), tol=0)
check("present", "ranged from ρ = −0.29 to +0.18, and only one of these thirty tests reached significance", True)
check("present", "(SARAF in lumbar cord, ρ = −0.194, q = 0.049)", True)
_apa11 = pd.read_csv(f"{P}/supplementary/S11_APA_candidate_gradients.csv")
_apa11 = _apa11[_apa11.dataset == "SH-SY5Y (GSE296712)"].set_index(["gene", "unit"])
for _g, _u, _v in [("SARAF", "termexon", -0.036), ("ORAI2", "termexon", 0.012),
                   ("SARAF", "intron5", 0.000)]:
    close(f"S11 holds the {_g} {_u} value quoted in Section 3.7", _v,
          _apa11.loc[(_g, _u), "delta"], tol=0.0005)
close("S11 units below the candidate threshold are flagged", 1,
      int((_apa11.candidate_gradient == "no").sum() > 0), tol=0)
for _s in ["the twelve core SOCE-pathway genes", "The only stringent-filter unannotated change anywhere in the nineteen-gene SOCE-associated set",
           "three terminal-exon units outside this set also excluded zero",
           "this model provides no evidence of enrichment beyond the length and expression properties",
           "Fisher’s exact test p = 0.48", "four independent datasets covering Alzheimer’s disease, Parkinson’s disease and multiple sclerosis",
           "In the two human models compared at junction level", "Figure 1A", "Figure 1B", "Figure 1C"]:
    check("present", _s, True)
for _s in ["twelve SOCE-machinery genes", "we did not detect high-confidence unannotated splicing changes in the store-operated entry machinery",
           "largely attributable to the length and expression properties",
           "the genes it identifies are.", "Three independent cohorts extended this"]:
    check("absent", _s, False)
_supp = _sp.run(["pandoc", "-t", "plain", "--wrap=none",
                 f"{P}/supplementary/SUPPLEMENTARY_MATERIAL.docx"], capture_output=True, text=True).stdout
close("supplementary document calls S12 the machine-readable version of Table 4", 1,
      int("S12. Machine-readable version of Table 4" in _supp), tol=0)
close("supplementary document describes the S11 candidate_gradient column", 1,
      int("candidate_gradient column" in _supp), tol=0)

out.append("\n=== 28 September 2026 author correction: Fura-2 wells on one plate ===")
from docx import Document as _DocxDocument
_main_doc = _DocxDocument(f"{P}/manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx")
_supp_doc = _DocxDocument(f"{P}/supplementary/SUPPLEMENTARY_MATERIAL.docx")
_highlights_doc = _DocxDocument(f"{P}/highlights_Neurochemistry_International.docx")
close("supplementary and main titles match", 1,
      int(_supp_doc.paragraphs[1].text == _main_doc.paragraphs[0].text), tol=0)
close("third highlight names the recurrent CBARP finding", 1,
      int(_highlights_doc.paragraphs[3].text ==
          "CBARP splicing changed in five of six TDP-43-depletion RNA-seq datasets."), tol=0)
close("highlight names the one-plate design", 1,
      int("One-plate Fura-2 measurements" in _highlights_doc.paragraphs[1].text), tol=0)
close("S1 description names the wells and plate", 1,
      int("three wells per group on one plate" in _supp_doc.paragraphs[12].text), tol=0)
_fp = pd.read_excel(f"{P}/supplementary/S1_laboratory_source_data.xlsx", sheet_name="Fura2_per_culture")
_c = _fp[_fp.group == "Non-targeting shRNA control"]; _k = _fp[_fp.group == "shTDP-43"]
_d = _k.readdition.mean() - _c.readdition.mean()
close("readdition difference (-1.30)", -1.30, _d, tol=0.005)
close("readdition control minimum", 1.013, _c.readdition.min(), tol=0.001)
close("readdition control maximum", 1.975, _c.readdition.max(), tol=0.001)
close("readdition knockdown minimum", 0.080, _k.readdition.min(), tol=0.001)
close("readdition knockdown maximum", 0.338, _k.readdition.max(), tol=0.001)
_d = _k.ER_release.mean() - _c.ER_release.mean()
close("ER release difference (-0.088)", -0.088, _d, tol=0.001)
close("ER release decrease (33%)", 33, 100 * (1 - _k.ER_release.mean() / _c.ER_release.mean()), tol=0.5)
close("readdition decrease (84%)", 84, 100 * (1 - _k.readdition.mean() / _c.readdition.mean()), tol=0.5)
_rc, _rk = _c.readdition_to_release_ratio, _k.readdition_to_release_ratio
close("ratio control mean (5.88)", 5.88, _rc.mean(), tol=0.005)
close("ratio control SEM (1.10)", 1.10, _rc.std(ddof=1) / np.sqrt(3), tol=0.005)
close("ratio knockdown mean (1.36)", 1.36, _rk.mean(), tol=0.005)
close("ratio knockdown SEM (0.01)", 0.01, _rk.std(ddof=1) / np.sqrt(3), tol=0.005)
close("ratio decrease (77%)", 77, 100 * (1 - _rk.mean() / _rc.mean()), tol=0.5)
close("every knockdown ratio below every control ratio", 1, int(_rk.max() < _rc.min()), tol=0)
for _s in ["three wells per group on one plate; difference −1.30",
           "The observed ranges were 1.013–1.975 and 0.080–0.338",
           "ER Ca²⁺ release was 0.268 ± 0.042 and 0.180 ± 0.061",
           "This ratio was 5.88 ± 1.10 in controls and 1.36 ± 0.01 after knockdown, 77% lower",
           "Each group comprised three wells on the same plate"]:
    check("present", _s, True)
for _s in ["Welch’s t-test p = 0.035", "Hedges’ g = −2.9", "n = 3 independent cultures",
           "ER release fell less and not significantly"]:
    check("absent", _s, False)
_j = pd.read_csv(f"{P}/source_data/CBARP_locus/CBARP_junction_counts_per_library.tsv", sep="\t")
_share = []
for (_ds, _g, _smp), _sub in _j.groupby(["dataset", "group", "sample"]):
    _tot = _sub.loc[_sub.end == 1235500, "count"].sum()
    _b = _sub.loc[(_sub.start == 1235343) & (_sub.end == 1235500), "count"].sum()
    _share.append((_ds, _g, _b / _tot))
_share = pd.DataFrame(_share, columns=["dataset", "group", "share"]).groupby(["dataset", "group"]).share.mean()
close("CBARP junction b share, SH-SY5Y control (17%)", 17, 100 * _share[("SH-SY5Y", "Control")], tol=0.5)
close("CBARP junction b share, SH-SY5Y knockdown (86%)", 86, 100 * _share[("SH-SY5Y", "TDP-43 KD")], tol=0.5)
close("CBARP junction b share, iPSC control (1%)", 1, 100 * _share[("iPSC colonies", "Control")], tol=0.5)
close("CBARP junction b share, iPSC knockdown (78%)", 78, 100 * _share[("iPSC colonies", "TDP-43 KD")], tol=0.5)
_sh = _j[_j.dataset == "SH-SY5Y"].groupby(["group", "start", "end"])["count"].sum()
close("SH-SY5Y junction b reads, knockdown (LSV test: 33)", 33, _sh[("TDP-43 KD", 1235343, 1235500)], tol=0)
close("SH-SY5Y junction b reads, control (LSV test: 26)", 26, _sh[("Control", 1235343, 1235500)], tol=0)
close("alternative 3' site to exon 5 acceptor distance (197 nt)", 197, 1235342 - 1235146 + 1, tol=0)
check("present", "197 nucleotides upstream of the exon 5 acceptor", True)
check("present", "carried 17% of exon-4 donor reads in SH-SY5Y controls and 86% after knockdown, and 1% and 78% in iPSC colonies", True)
_first = {}
for _m in re.finditer(r"Figures? (\d)", SEARCH_TXT):
    _first.setdefault(int(_m.group(1)), _m.start())
close("main figures first cited in numerical order", 1,
      int([k for k, _ in sorted(_first.items(), key=lambda kv: kv[1])] == sorted(_first)), tol=0)
for _s in ["All laboratory comparisons are descriptive", "These Fura-2 and WST-1 comparisons are descriptive",
           "What reproduces is the involvement of the locus", "Figure5_CBARP_splicing",
           "UNC13A programs", "cryptic-splicing program."]:
    check("absent", _s, False)

out.append("\n=== cross-file agreement and the corrected correction family ===")
# Verify the accession inventory against the sample count stated in the manuscript.
_s17 = pd.read_csv(f"{P}/supplementary/S17_dataset_accessions.csv")
_ms_design = _s17.loc[_s17.accession == "GSE138614", "design"]
close("S17 GSE138614 row count", 1, len(_ms_design), tol=0)
close("S17 GSE138614 sample count", 98,
      int(re.search(r"(\d+) samples", _ms_design.iloc[0]).group(1)) if len(_ms_design) else -1, tol=0)
_s3 = pd.read_csv(f"{P}/supplementary/S3_rMATS_significant_events.csv.gz", compression="gzip")
_event_coords = {
    "SE": ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "A5SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "A3SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "MXE": ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd",
            "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "RI": ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
}
for _ev, _coords in _event_coords.items():
    _rows = _s3[_s3.event_class == _ev]
    close(f"S3 {_ev} coordinate completeness", 0, int(_rows[_coords].isna().sum().sum()), tol=0)
    close(f"S3 {_ev} unique event IDs", 0,
          int(_rows.duplicated(["dataset", "model", "ID"]).sum()), tol=0)
close("S3 inclusion and skipping form lengths", 0,
      int((_s3.IncFormLen.isna() | _s3.SkipFormLen.isna() |
           (_s3.IncFormLen <= 0) | (_s3.SkipFormLen <= 0)).sum()), tol=0)
# the same GSE138614 donor-level test appears in Table 5 and in S16b; the two pipelines
# must share the low-expression filter, or the per-sample median centring shifts the delta
_t5 = pd.read_csv(f"{P}/tables/Table5_cross_disease_comparison.csv")
_t5d = _t5[(_t5.group == "MS") & (_t5.region.str.startswith("Donor level")) & (_t5.gene == "TRPC1")]
_s16 = pd.read_csv(f"{P}/supplementary/S16b_multiple_sclerosis_donor_level.csv")
_s16d = _s16[(_s16.comparison == "MS vs control, raw TRPC1") & (_s16.unit == "donor")]
close("MS donor-level delta: Table 5 vs S16b", float(_t5d.cliffs_delta.iloc[0]),
      float(_s16d.cliffs_delta.iloc[0]), tol=0.001)
# S9: the eleven self-correlations of gene-level STMN2 carry no information and are
# excluded from the Benjamini-Hochberg family
_s9 = pd.read_csv(f"{P}/supplementary/S9_cryptic_PSI_correlations_within_ALS.csv")
close("S9 tests in the correction family", 220, int((_s9.in_correction_family == "yes").sum()), tol=0)
close("S9 trivial tests excluded", 11, int((_s9.in_correction_family == "no").sum()), tol=0)
close("S9 excluded tests have no q value", 0, int(_s9.loc[_s9.in_correction_family == "no", "q_value"].notna().sum()), tol=0)
_sig = _s9[(_s9.proxy == "cryptic STMN2 PSI") & (_s9.q_value < 0.05)].copy()
_sig["tag"] = _sig.target_gene + " " + _sig.region
for _gene, _reg, _q in [("STIM1", "Spinal Cord Lumbar", 0.0021), ("ORAI1", "Spinal Cord Cervical", 0.0031),
                        ("STMN2", "Spinal Cord Cervical", 0.011), ("ORAI1", "Spinal Cord Lumbar", 0.013),
                        ("ATP2A2", "Cortex Motor Medial", 0.015), ("GFAP", "Cortex Motor Medial", 0.023),
                        ("SNAP25", "Cortex Motor Medial", 0.028), ("SARAF", "Spinal Cord Lumbar", 0.049)]:
    _r = _sig[(_sig.target_gene == _gene) & (_sig.region == _reg)]
    close(f"S9 q, {_gene} {_reg}", _q, float(_r.q_value.iloc[0]) if len(_r) else -9, tol=0.0006)

out.append("\n=== published tables carry English dataset labels and the full control set ===")
_POS16 = {"STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
          "ATG4B", "SETD5", "CAMK2B", "ELAVL3", "POLDIP3", "RSF1", "GPSM2", "SYNJ2"}
for _f in ["S5_stringent_filter_unannotated_splicing_candidates", "S6_cryptic_positive_controls",
           "S7_SOCE_genes_annotation_free"]:
    _d = pd.read_csv(f"{P}/supplementary/{_f}.csv")
    _bad = [c for c in _d.comparison.unique()
            if re.search(r"koloni|DOZ|Fare|^[A-Z0-9]*_", str(c))]
    close(f"{_f}: untranslated dataset labels", 0, len(_bad), tol=0)
_s6 = pd.read_csv(f"{P}/supplementary/S6_cryptic_positive_controls.csv")
close("S6 genes all belong to the sixteen literature controls", 0,
      len({g for g in _s6.gene if g.upper() not in _POS16}), tol=0)
close("S6 contains GPSM2, named in Section 3.3", 1, int("GPSM2" in set(_s6.gene)), tol=0)
_s4 = pd.read_csv(f"{P}/supplementary/S4_matched_permutation_enrichment.csv")
close("S4 decision values are English", 0,
      len(set(_s4.decision) - {"enrichment", "no enrichment"}), tol=0)
close("S4 panels with surviving enrichment", 2, int((_s4.decision == "enrichment").sum()), tol=0)
_t3 = pd.read_csv(f"{P}/tables/Table4_cryptic_events_eleven_comparisons.csv")
close("Table 4 rows = eleven comparisons", 11, len(_t3), tol=0)
_mq = _s6[_s6.comparison == "SH-SY5Y 75 ng/mL (MAPQ-filtered set)"]
close("SH-SY5Y MAPQ-filtered positive controls", 12, _mq.gene.nunique(), tol=0)

out.append("\n=== Section 3.2 CBARP figures and Section 3.5 completeness ===")
_f = pd.read_pickle(f"{M}/03_TABLOLAR/_filtreli_olaylar.pkl")
_cb = _f[_f.geneSymbol.astype(str).str.upper() == "CBARP"]
_cb = _cb[(_cb.FDR < 0.05) & (_cb.IncLevelDifference.abs() >= 0.10)].copy()
_cb["reads"] = _cb.apply(lambda r: sum(int(x) for c in
    ["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2"]
    for x in str(r[c]).split(",")), axis=1)
close("CBARP coverage-qualified events", 32, len(_cb), tol=0)
close("CBARP datasets", 5, _cb.dataset.nunique(), tol=0)
close("CBARP |dPSI| minimum", 0.11, _cb.IncLevelDifference.abs().min(), tol=0.005)
close("CBARP |dPSI| maximum", 0.74, _cb.IncLevelDifference.abs().max(), tol=0.005)
close("CBARP |dPSI| median", 0.37, _cb.IncLevelDifference.abs().median(), tol=0.005)
close("CBARP reads minimum", 61, _cb.reads.min(), tol=0)
close("CBARP reads maximum", 4559, _cb.reads.max(), tol=0)
close("CBARP reads median", 246, _cb.reads.median(), tol=0.5)
close("CBARP events above 100 reads", 26, int((_cb.reads >= 100).sum()), tol=0)
_shsy = _cb[_cb.dataset == "GSE296712_SHSY5Y"]
close("CBARP SH-SY5Y events, all negative", 3, int((_shsy.IncLevelDifference < 0).sum()), tol=0)
# every SH-SY5Y unit whose bootstrap interval excludes zero must be named in Section 3.5
_apa = pd.read_csv(f"{M}/07_DISK_ANALIZLERI/sonuclar/APA_corrected_full_core_summary.tsv", sep="\t")
_POS16b = {"STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
           "ATG4B", "SETD5", "CAMK2B", "ELAVL3", "POLDIP3", "RSF1", "GPSM2", "SYNJ2"}
_ex = _apa[((_apa.boot_low > 0) | (_apa.boot_high < 0)) & ~_apa.gene.str.upper().isin(_POS16b)]
_unnamed = sorted({g for g in _ex.gene if f"*{g}*" not in TXT.split("## **References**")[0]})
close("Ca2+ units excluding zero that Section 3.5 does not name", 0, len(_unnamed), tol=0)
out.append("    " + ("all named: " + ", ".join(sorted(set(_ex.gene)))) if not _unnamed
           else "    unnamed: " + ", ".join(_unnamed))
_s10b = pd.read_csv(f"{P}/supplementary/S10b_NMD_panel_descriptive_summary.csv").set_index("panel")
close("S10b cryptic reference panel descriptive median", -0.1143427734,
      float(_s10b.loc["Cryptic_splicing_reference_genes_16", "median_interaction_log2"]), tol=1e-8)
close("S10b omits inferential tests", 0,
      len([c for c in _s10b if c.startswith(("p_", "q_"))]), tol=0)

out.append("\n=== manuscript strings that must be present ===")
for s in ["10,926 versus 176 reads", "three spinal cord levels", "six brain regions",
          "0.64 calls in iPSC colonies, 2.17 in K562 total RNA and 0.83 in mouse striatum",
           "chr11:4,088,702–4,088,738",
          "chr4:27,007,983–27,008,006", "n = 3", "10 µM cyclopiazonic acid",
          "Albarran L, Lopez JJ, Woodard GE, Salido GM, Rosado JA",
          "cryptic junction was most frequent in spinal cord",
          "Cutadapt v5.2",
          "59 units in 29 genes passed the depth filter",
          "from 0.570 in controls to 0.819 in knockdown (Δ = +0.249",
          "`-p --countReadPairs` for paired-end libraries",
          "74 qualifying units in C2C12 and 131 in NSC34",
          
          "null-to-real call ratios were 0.98 in iPSC colonies and 2.01 in K562 total RNA",
          "this check is available only for the human comparisons",
          "isoform-level testing did not detect a *CBARP* isoform switch",
          "the three mouse comparisons have no conserved control",
          "*ORAI3* went from 8% to 30% of ORAI transcripts and *ORAI2* from 77% to 53%",
          "The STIM (+48%) and Ca²⁺-entry-regulator (+30%) pools also rose, whereas the SERCA and mitochondrial-uptake pools fell (−12% and −25%)",
          "the median gene expressed in both groups (mean TPM > 5) had 24% lower TPM",
          "Neither is independent of the splicing results",
          "Yoast RE, Emrich SM, Zhang X, et al.",
          "3\u2076 = 729 combinations for the three-versus-three comparisons and 2\u2074 = 16",
          "Total *STMN2* expression was therefore not used as a specific indicator",
          "Cellular metabolic activity was assessed 48 h after seeding",
          "10 MS/5 control donors, 98 samples",
          
          "decreased at donor level across all sampled lesion types (\u03b4 = \u22120.840; q = 0.038)"]:
    check("present", s, True)
out.append("\n=== strings that must be gone ===")
for s in ["Sah P, et al.", "10,930", "four spinal cord regions", "p = 0.13)", "log2FC = −1.746",
          "prioritized candidate", "full GENCODE v47 index", "Cutadapt v4.6",
          "the mouse datasets were not analysed at all",
          "all 3\u2076 = 729 replicate combinations", "their p values are therefore optimistic",
          "applied across all 231 correlations", "is not measurable in SH-SY5Y",
          "used in the first version of this analysis", "The original analysis treated eight contrasts",
          "the correction was applied only to the primary SH-SY5Y model",
          # 22 September 2026 audit
          "three independent cultures per group", "Ca²⁺ imaging",
          "isoform-level testing, to converge", "each of which is a minor member",
          "lie outside the core panel", "about two-thirds as many calls",
          "a second control moved in the expected direction", "silent in SH-SY5Y",
          "the dominant members fall", "in the second junction set",
          "in ten of the eleven comparisons", "Length-corrected family summaries",
          "(−32% and −43%", "(−3% and −1%)", "informative in the two neuronal models",
          "the current transcript estimates", "the positive control was informative only in the motor-neuron-like line",
          # revision-history wording removed on 22 September 2026
          "no longer", "before the correction", "The corrected analysis", "the corrected analyses",
          "corrected full-panel", "Corrected APA", "After correcting the APA", "Our first attempt",
          "we withdraw that inference", "our own prior candidate", "the original call",
          "remained the largest", "also retained positive", "is not retained", "Our own *TRPC1* candidate",
          # the laboratory control group is named, and the apoptosis markers are not part of the study
          "with the control set to 1.0", "normalised to the mean of the control group", "BCL2", "BAK1"]:
    check("absent", s, False)

out.append("\n=== every figure and supplementary table is cited in the text ===")
_body = SEARCH_TXT.split("## References")[0]
_first = {}
for _i in range(1, 5):
    _m = re.search(rf"Figure {_i}(?![0-9])", _body)
    close(f"Figure {_i} cited in the text", 1, int(_m is not None), tol=0)
    _first[_i] = _m.start() if _m else -1
close("figures first cited in numerical order", 1,
      int(all(_first[i] < _first[i + 1] for i in range(1, 4))), tol=0)
_results = _body.split("## 3. Results")[1].split("## Declarations")[0]
_results = re.sub(r"<figcaption>.*?</figcaption>", " ", _results, flags=re.S)
for _n in range(1, 7):
    close(f"Figure {_n} cited in the Results or Discussion text", 1,
          int(re.search(rf"Figure {_n}[A-E]?(?![0-9])", _results) is not None), tol=0)
for _i in range(1, 18):
    close(f"Supplementary Table S{_i} cited in the text", 1,
          int(re.search(rf"Table S{_i}(?![0-9])", _body) is not None), tol=0)

out.append("\n=== S1 carries the four reported targets and names the control group ===")
_s1 = pd.ExcelFile(f"{P}/supplementary/S1_laboratory_source_data.xlsx")
_tg = pd.read_excel(_s1, "Target_qPCR_Ct")
close("S1 target genes are the four reported targets", 1,
      int(sorted(_tg.gene.unique()) == ["ATP2A3", "ORAI1", "STIM1", "TRPC1"]), tol=0)
for _sh in ["Target_qPCR_Ct", "Target_qPCR_rel", "Fura2", "WST1"]:
    close(f"S1 {_sh}: control group is the non-targeting shRNA control", 1,
          int(set(pd.read_excel(_s1, _sh).group) == {"Non-targeting shRNA control", "shTDP-43"}), tol=0)
check("Methods name the comparison group", "compared shTDP-43 cells with the non-targeting shRNA control", True)
# Author decisions: Y. Kaymaz is not an author; only the available 48-h WST-1
# 48-h WST-1 dataset is reported as n = 4.
check("author list", "**Elmasnur Yılmazᵃ, Yasemin Eraçᵃ,\\***", True)
# confirmed by the author on 23 September 2026: the affiliation is Pharmacology (the thesis
# belongs to the Biotechnology PhD programme, as the acknowledgement says), and the SH-SY5Y
# medium contained no antibiotic ("penicillin" is among the forbidden strings above)
check("affiliation", "ᵃ Department of Pharmacology, Faculty of Pharmacy, Ege University", True)
check("thesis programme", "Graduate School of Natural and Applied Sciences, Department of Biotechnology, 2026", True)
for _s in ["Kaymaz", "Y.K.", "Bioengineering", "performed three times"]:
    check("absent", _s, False)
for _s in ["n = 4", "n = 3"]:
    check("present", _s, True)
_rd = pd.read_excel(_s1, "README")
close("S1 README omits unsupported WST-1 experiments", 0,
      int(_rd.astype(str).apply(lambda c: c.str.contains("performed three times")).any().any()), tol=0)

out.append("\n=== 22 September 2026, round 3: values quoted from the supplements ===")
_d16 = pd.read_csv(f"{P}/supplementary/S16b_multiple_sclerosis_donor_level.csv")
_r = lambda comp, unit: _d16[(_d16.comparison == comp) & (_d16.gene == "TRPC1") & (_d16.unit == unit)].iloc[0]
close("MS NAWM sample-level delta (-0.482)", -0.482, _r("NAWM vs control WM", "sample").cliffs_delta)
close("MS NAWM sample-level p (0.005)", 0.005, _r("NAWM vs control WM", "sample").p, tol=0.0006)
_t5n = t5[(t5.gene == "TRPC1") & t5.region.str.startswith("Normal-appearing")]
close("MS NAWM sample-level q (0.10)", 0.10, float(_t5n.q_value.iloc[0]), tol=0.005)
close("MS raw sample-level delta (-0.562)", -0.562, _r("MS vs control, raw TRPC1", "sample").cliffs_delta)
close("MS adjusted sample-level delta (-0.418)", -0.418,
      _r("MS vs control, myelin+glia-adjusted TRPC1", "sample").cliffs_delta)
close("MS adjusted sample-level p (0.002)", 0.002,
      _r("MS vs control, myelin+glia-adjusted TRPC1", "sample").p, tol=0.0006)
# STIM2.1 exon: fixed-effect inverse-variance meta-analysis of the six rMATS estimates (S15)
_s15 = pd.read_csv(f"{P}/supplementary/S15_STIM2.1_exon_six_datasets.csv")
_w = 1 / _s15.variance
_est = float((_w * _s15.delta_PSI).sum() / _w.sum())
_se = float(np.sqrt(1 / _w.sum()))
from math import erf, sqrt
_pm = 2 * (1 - 0.5 * (1 + erf(abs(_est / _se) / sqrt(2))))
close("STIM2.1 pooled dPSI (+0.0013)", 0.0013, _est, tol=0.00006)
close("STIM2.1 pooled CI low (-0.022)", -0.022, _est - 1.96 * _se, tol=0.0006)
close("STIM2.1 pooled CI high (+0.024)", 0.024, _est + 1.96 * _se, tol=0.0006)
close("STIM2.1 pooled p (0.914)", 0.914, _pm, tol=0.0006)
close("STIM2.1 exon positive in five of six datasets", 5, int((_s15.delta_PSI > 0).sum()), tol=0)
from scipy import stats as _st
_Q = float((_w * (_s15.delta_PSI - _est) ** 2).sum())
close("STIM2.1 Cochran's Q (4.84)", 4.84, _Q, tol=0.006)
close("STIM2.1 heterogeneity p (0.44)", 0.44, float(_st.chi2.sf(_Q, len(_s15) - 1)), tol=0.006)
close("STIM2.1 I2 (0%)", 0, max(0.0, (_Q - (len(_s15) - 1)) / _Q), tol=0)
close("STIM2.1 weight carried by the mouse datasets (86%)", 0.86,
      float(_w[_s15.species == "mouse"].sum() / _w.sum()), tol=0.006)
_h = _s15[_s15.species == "human"]
_wh = 1 / _h.variance
_eh = float((_wh * _h.delta_PSI).sum() / _wh.sum())
_sh = float(np.sqrt(1 / _wh.sum()))
close("STIM2.1 human-only pooled dPSI (+0.031)", 0.031, _eh, tol=0.0006)
close("STIM2.1 human-only CI low (-0.030)", -0.030, _eh - 1.96 * _sh, tol=0.0006)
close("STIM2.1 human-only CI high (+0.091)", 0.091, _eh + 1.96 * _sh, tol=0.0006)
close("STIM2.1 human-only p (0.32)", 0.32, float(2 * _st.norm.sf(abs(_eh / _sh))), tol=0.006)
close("enrichment tests: dataset-panel combinations (24)", 24, len(_s4), tol=0)
_mn4 = _s4[_s4.dataset == "GSE77702_iPSC_MN"].set_index("panel").p_matched_permutation
close("iPSC-MN channel/transport matched p (0.046)", 0.046, _mn4["Tier2_Channel_Release_Transport_117"], tol=0.0006)
close("iPSC-MN curated Ca2+ matched p (0.032)", 0.032, _mn4["Tier3_Curated_Calcium_Handling_258"], tol=0.0006)
_abs = TXT.split("## **Abstract**")[1].split("**Keywords:**")[0]
close("abstract length, words including headings (at most 350)", 1,
      int(len(re.sub(r"[*]", "", _abs).split()) <= 350), tol=0)
for s_ in ["fixed-effect inverse-variance meta-analysis",
           "the mean per-base depths of its two windows summed to at least 3",
           "derived from the doctoral thesis of Elmasnur Yılmaz",
           
           
           "375,000 per well", "Dharmacon TRC Lentiviral shRNA, cat. no. RHS3979",
           "1 mM EGTA", "Premix WST-1, Takara Bio, cat. no. MK400"]:
    check("present", s_, True)
for s_ in ["penicillin", "dominant-negative", "the donor-level values are the ones reported",
           "is expressed in oligodendrocytes", "log2(TDP-43 KD", "the reduction held",
           "was still reduced", "summed depth of its two windows",
           "while leaving it intact in iPSC-derived motor neurons",
           "the *STIM2.1*/STIM2β question"]:
    check("absent", s_, False)

out.append("\n=== every in-text citation has a reference and every reference is cited ===")
import unicodedata as _ud
_refs_txt = TXT.split("## **References**")[1].split("**S1.**")[0]
_refs = [r.strip() for r in _refs_txt.strip().split("\n\n") if r.strip()]
_key = lambda t: "".join(c for c in _ud.normalize("NFKD", t) if not _ud.combining(c)).lower()
close("reference list in alphabetical order", 1, int([_key(r) for r in _refs] == sorted(_key(r) for r in _refs)), tol=0)
_ref_keys = set()
for _r_ in _refs:
    # the first author ends where the initials begin (initials may be non-ASCII, e.g. "Selli Ç")
    _m = re.match(r"(.+?) [A-Z\u00c0-\u017d][A-Z\u00c0-\u017d\-]*[,.]", _r_)
    _y = re.search(r"\b((?:19|20)\d\d)[;.]", _r_)
    if _m and _y:
        _ref_keys.add((_m.group(1).split(";")[0].strip(), _y.group(1)))
if any("National Center for Biotechnology Information (NCBI)" in r and "2025" in r for r in _refs):
    _ref_keys.add(("NCBI", "2025"))
_cites = set()
_NAME = r"(?:Van den |Van |De )?[A-ZÀ-Ž][\w’'\-]+"
for _g in re.findall(r"\(([^()]*\d{4}[^()]*)\)", _body):
    for _piece in _g.split(";"):
        _m = re.match(rf"^\s*(?:.*, )?({_NAME})(?: et al\.| and {_NAME})?, ((?:19|20)\d\d)\s*$", _piece)
        if _m:
            _cites.add((_m.group(1), _m.group(2)))
for _m in re.finditer(rf"({_NAME})(?: et al\.| and {_NAME})? \(((?:19|20)\d\d)\)", _body):
    _cites.add((_m.group(1), _m.group(2)))
_missing = sorted(c for c in _cites if c not in _ref_keys)
_uncited = sorted(r for r in _ref_keys if r not in _cites)
close("in-text citations without a reference entry", 0, len(_missing), tol=0)
close("reference entries never cited in the text", 0, len(_uncited), tol=0)
if _missing:
    out.append("    missing: " + "; ".join(f"{a} {b}" for a, b in _missing))
if _uncited:
    out.append("    uncited: " + "; ".join(f"{a} {b}" for a, b in _uncited))
out.append(f"    {len(_cites)} distinct citations, {len(_refs)} references")

out.append("\n=== Tables 1-5 in the manuscript match tables/*.csv ===")
# The tables are regenerated from tables/*.csv with the formatting functions of the baseline builder
# and compared with the tables of the CURRENT manuscript by their numbers, row by row. (Until v1.0.7 this
# check regenerated the historical v4 Markdown instead of the current manuscript.)
import importlib.util as _iu
_spec = _iu.spec_from_file_location("baseline_tables", f"{P}/code/build_manuscript_docx.py")
_bt = _iu.module_from_spec(_spec)
_spec.loader.exec_module(_bt)
_NUM = re.compile(r"[+−\-]?\d[\d,]*(?:\.\d+)?(?:\s*×\s*10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)?")


def _table_rows(tbl):
    rows = []
    for line in tbl.splitlines():
        if line.startswith("|") and not set(line.replace("|", "").strip()) <= set(":- "):
            rows.append([c.strip() for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))])
    return rows


def _numbers(cells):
    txt = " ".join(cells).replace("*", "").replace("\\", "")
    return [re.sub(r"\s+", "", x).replace("−", "-") for x in _NUM.findall(txt)]


for _n, _fn in {1: _bt.table4, 2: _bt.table1, 3: _bt.table2, 4: _bt.table3, 5: _bt.table5}.items():
    _cap = re.search(rf"^\*\*Table {_n}\.", TXT, re.M)
    _rest = TXT[_cap.end():]
    _i = _rest.index("|")
    _act = _table_rows(_rest[_i:_rest.index("\n\n", _i)])[1:]
    _exp = _table_rows(_fn()[0])[1:]
    _bad = sum(_numbers(e) != _numbers(a) for e, a in zip(_exp, _act)) + abs(len(_exp) - len(_act))
    close(f"Table {_n} of the manuscript equals the numbers regenerated from tables/*.csv ({len(_exp)} rows)", 0, _bad, tol=0)
for _n in range(1, 6):
    check(f"table {_n} present", f"Table {_n}.", True)


out.append("\n=== release v1.0.8: evaluation of v1.0.7 (claim, naming, simulation, figures, package) ===")
# -- naming: labels that overstated the evidence or blurred the biology are gone from every document
_TERMS = ["high-confidence", "SOCE machinery", "SOCE-machinery", "core entry components", "SOCE regulator",
          "SOCE-regulator", "discordance", "did most of the work", "We first measure SOCE", "independently supported",
          "Integrated mechanism", "Which splicing changes", "observed SOCE phenotype", "central observation"]
for _t in _TERMS:
    close(f"'{_t}' is absent from the manuscript", 0, int(_t.lower() in SEARCH_TXT.lower()), tol=0)
    close(f"'{_t}' is absent from the supplementary document", 0, int(_t.lower() in _supp2.lower()), tol=0)
for _s in ["stringent-filter unannotated splicing candidate", "nineteen-gene SOCE-associated set",
           "core SOCE-pathway genes", "Ca²⁺-entry regulators",
           "The name describes the filter and does not imply validated specificity",
           "a set of candidates whose sensitivity has been checked against known targets, not a set with established specificity"]:
    check("present", _s, True)

# -- the claim: a single-plate discovery observation, not paired with the RNA data
_NEW_TITLE = ("Calcium-regulatory RNA candidates after TDP-43 knockdown in SH-SY5Y cells: reanalysis of public "
              "RNA-seq data and a single-plate Ca²⁺-readdition observation")
close("title states the RNA-candidate framing and the single plate", 1, int(_main_doc.paragraphs[0].text == _NEW_TITLE), tol=0)
for _s in ["this observation motivated the RNA analyses but was not paired with the RT-qPCR or RNA-seq samples",
           "the Fura-2 observation comes from a single culture plate and served to motivate the RNA analyses, not to establish a phenotype",
           "The RT-qPCR mRNA measurements and the one-plate Fura-2 response therefore differed in direction",
           "they are not paired and do not show opposing changes within the same cells",
           "Working hypotheses and discriminating experiments",
           "Possible contributors to the single-plate Fura-2 difference",
           "none of which was tested here",
           "no confirmed RNA-processing event in the core SOCE-pathway genes explained it"]:
    check("present", _s, True)
close("Results 3.1 heading names the single plate", 1,
      int("3.1 SOCE-associated mRNAs are higher, and the Ca²⁺-readdition amplitude is lower in one Fura-2 plate" in SEARCH_TXT), tol=0)

# -- Discussion order: RNA findings first, the single-plate hypotheses after them, tissue and experiments last
_disc = SEARCH_TXT.split("## 4. Discussion")[1].split("### Limitations")[0]
_order = [_disc.index(k) for k in ["This study asked which calcium-regulatory RNA changes",
                                   "RNA-processing findings depend on model and method",
                                   "CBARP and the other splicing candidates",
                                   "Possible contributors to the single-plate Fura-2 difference",
                                   "Comparison with ALS calcium models and cell state",
                                   "Relevance to disease tissue", "Working hypotheses and discriminating experiments"]]
close("Discussion paragraphs follow the intended order", 1, int(_order == sorted(_order)), tol=0)

# -- section numbering after the merge of the NMD and outlier results into Section 3.6
for _h in ["3.6 APA, NMD and outlier screens yield candidates but no confirmed event in the core SOCE-pathway genes",
           "3.7 Matched analyses do not establish enrichment", "3.8 TRPC1, SARAF and CBARP differ across ALS and neurological comparison cohorts"]:
    check("present", _h, True)
for _h in ["3.9 TRPC1", "Section 3.9", "3.6 Descriptive NMD", "A coverage-based APA screen yields"]:
    check("absent", _h, False)
_secs = {int(x) for x in re.findall(r"Sections? 3\.(\d+)", _body)}
close("every 'Section 3.x' reference points to an existing section (3.1-3.8)", 1, int(max(_secs) <= 8), tol=0)
_meth = {int(x) for x in re.findall(r"(?:Methods|Sections?) 2\.(\d+)", _body)}
close("every 'Methods 2.x' reference points to an existing section (2.1-2.16)", 1, int(max(_meth) <= 16), tol=0)

# -- detection-power simulation: values in the text come from the stored table, assumptions from the code
_pw = pd.read_csv(f"{P}/source_data/power_simulation_S1.csv")
_pv = lambda d, x: float(_pw[(_pw.reads_per_sample == d) & (_pw.true_delta_PSI.round(2) == x)].power.iloc[0])
for _lab, _d, _x, _v in [("10 reads, dPSI 0.10", 10, 0.10, 0.08), ("100 reads, dPSI 0.10", 100, 0.10, 0.28),
                         ("50 reads, dPSI 0.20", 50, 0.20, 0.58), ("100 reads, dPSI 0.20", 100, 0.20, 0.74),
                         ("50 reads, dPSI 0.30", 50, 0.30, 0.90), ("100 reads, dPSI 0.30", 100, 0.30, 0.97)]:
    close(f"simulated power, {_lab}", _v, _pv(_d, _x), tol=0.005)
close("80% power is not reached at 10 or 20 reads up to dPSI 0.30", 1,
      int(_pw[_pw.reads_per_sample.isin([10, 20])].power.max() < 0.8), tol=0)
close("80% power is reached between dPSI 0.20 and 0.30 at 50 and 100 reads", 1,
      int(all(_pv(d, 0.20) < 0.8 <= _pv(d, 0.30) for d in (50, 100))), tol=0)
_sim = re.search(r"def sim_power\(n_per_group, depth, dpsi, base_psi=([0-9.]+), nsim=(\d+), alpha=([0-9.]+), disp=([0-9.]+)\)",
                 io.open(f"{P}/code/05_meta_permutasyon_guc.py", encoding="utf-8").read())
close("simulation baseline PSI in the code (0.5)", 0.5, float(_sim.group(1)))
close("simulation runs in the code (2,000)", 2000, int(_sim.group(2)), tol=0)
close("simulation alpha in the code (0.05)", 0.05, float(_sim.group(3)))
close("simulation between-replicate SD in the code (0.05)", 0.05, float(_sim.group(4)))
_ev = pd.read_csv(f"{P}/supplementary/S3_rMATS_significant_events.csv.gz")
_rd = _ev[(_ev.dataset == "SH-SY5Y (GSE296712)") & (_ev.model == "JC")].mean_reads_per_sample
close("lower quartile of reads per sample, nominally significant SH-SY5Y events (9.7)", 9.7, _rd.quantile(0.25), tol=0.05)
close("median reads per sample, nominally significant SH-SY5Y events (21.3)", 21.3, _rd.median(), tol=0.05)
for _s in ["close to the lower quartile of the nominally significant events (9.7 reads per sample; median 21.3)",
           "The simulation is not an estimate of rMATS power",
           "standard deviation 0.05, limited to 0.01–0.99",
           "The simulation does not reproduce the rMATS model, its FDR correction or its |ΔPSI| threshold"]:
    check("present", _s, True)
close("Supplementary Figure S1 legend states the assumptions", 1,
      int("Example detection-power simulation" in _supp2 and "is not an estimate of rMATS power" in _supp2), tol=0)

# -- svaseq numbers live in Results 3.2 and Supplementary Results 5 (table recomputed from the CSV)
_sv = pd.read_csv(f"{P}/source_data/svaseq_sensitivity_SHSY5Y.csv").set_index("gene")
_SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def _qval(txt):
    m = re.match(r"([0-9.]+)\s*×\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", txt)
    return float(m.group(1)) * 10 ** int(m.group(2).translate(_SUP)) if m else float(txt)


_tb = _supp_doc.tables[0]
close("svaseq table lists nine genes", 9, len(_tb.rows) - 1, tol=0)
for _r in _tb.rows[1:]:
    _g = _r.cells[0].text
    close(f"svaseq table, {_g}: log2FC original", float(_sv.loc[_g, "log2FC_published"]), float(_r.cells[1].text.replace("−", "-")), tol=0.006)
    close(f"svaseq table, {_g}: log2FC with surrogate variables", float(_sv.loc[_g, "log2FC_with_SV"]), float(_r.cells[3].text.replace("−", "-")), tol=0.006)
    for _c, _col in ((2, "padj_published"), (4, "padj_with_SV")):
        _x, _y = float(_sv.loc[_g, _col]), _qval(_r.cells[_c].text)
        close(f"svaseq table, {_g}: {_col}", 1, int(abs(_x - _y) <= 0.06 * _x + 5e-4), tol=0)

# -- STIM1 transcripts: identifiers only in the supplement; every number from the source files
_el = pd.read_csv(f"{P}/source_data/STIM1_transcript_test_eligibility.csv")
close("all nine STIM1 transcript identifiers are listed in Supplementary Results 6", 9, sum(f in _supp2 for f in _el.feature_id), tol=0)
close("no transcript identifier remains in the main text", 0, len(re.findall(r"ENST0000\d+", SEARCH_TXT)), tol=0)
_isa = pd.read_csv(f"{P}/source_data/STIM1_isoformswitch_DEXSeq.csv").set_index("isoform_id")
close("STIM1 isoform switch q, ENST00000698913.1 (5.5e-38)", 5.5e-38, _isa.loc["ENST00000698913.1", "isoform_switch_q_value"], tol=0.06e-38)
close("STIM1 isoform switch q, ENST00000698912.1 (1.6e-10)", 1.6e-10, _isa.loc["ENST00000698912.1", "isoform_switch_q_value"], tol=0.06e-10)
close("the two changing isoforms are not PTC isoforms", 1, int((_isa.loc[["ENST00000698912.1", "ENST00000698913.1"], "PTC"] == False).all()), tol=0)
close("PTC isoform q 0.33 (ENST00000698919.1)", 0.33, _isa.loc["ENST00000698919.1", "isoform_switch_q_value"], tol=0.005)
close("PTC isoform q 0.77 (ENST00000698918.1)", 0.77, _isa.loc["ENST00000698918.1", "isoform_switch_q_value"], tol=0.005)
close("dIF of the two changing isoforms (+0.037, -0.032)", 1,
      int(abs(_isa.loc["ENST00000698912.1", "dIF"] - 0.037) < 0.001 and abs(_isa.loc["ENST00000698913.1", "dIF"] + 0.032) < 0.001), tol=0)
for _s in ["ENST00000698912.1", "ENST00000698913.1", "dIF +0.037 and −0.032", "not a single confirmation chain"]:
    close(f"Supplementary Results 6 contains '{_s}'", 1, int(_s in _supp2), tol=0)

# -- DESeq2 count sources: every statistic in the article is from the Salmon-based table of the primary comparison;
#    the featureCounts run of the same six libraries is a cross-check whose agreement is quoted in Methods 2.2
_fc = pd.read_csv(f"{P}/source_data/DESeq2_ctrl_vs_75_featureCounts.csv", index_col=0)
_sa = pd.read_csv(f"{P}/source_data/DESeq2_ctrl_vs_75_fullmap.csv", index_col=0)
_cm = _fc.index.intersection(_sa.index)
_fa, _fb = _fc.loc[_cm, "log2FoldChange"], _sa.loc[_cm, "log2FoldChange"]
_fok = _fa.notna() & _fb.notna()
close("genes counted in both DESeq2 runs (20,833)", 20833, int(_fok.sum()), tol=0)
close("log2FC correlation, featureCounts versus Salmon counts (r = 0.85)", 0.85, float(np.corrcoef(_fa[_fok], _fb[_fok])[0, 1]), tol=0.005)
for _g, _v in (("STIM1", 0.929), ("TRPC1", 0.958), ("ORAI1", 0.433), ("ATP2A3", 1.306), ("CBARP", -1.254),
               ("ORAI3", 2.056), ("ORAI2", -0.173), ("ATP2A2", -0.247), ("GAPDH", 0.612)):
    close(f"quoted DESeq2 log2FC of {_g} is the Salmon-based value", _v, float(_sa.loc[_g, "log2FoldChange"]), tol=0.0006)
close("the four RT-qPCR targets have the same direction in both count sets", 1,
      int(all(np.sign(_fc.loc[_g, "log2FoldChange"]) == np.sign(_sa.loc[_g, "log2FoldChange"])
              for _g in ("STIM1", "TRPC1", "ORAI1", "ATP2A3"))), tol=0)
for _s in ["every DESeq2 statistic reported here was computed on Salmon estimated counts summed to genes with the complete transcript-to-gene map",
           "gave similar log2 fold changes (Pearson r = 0.85 across the 20,833 genes counted in both)",
           "log2FC and p_adj are from DESeq2 on gene-level Salmon counts"]:
    check("present", _s, True)
check("absent", "Gene counts were generated with featureCounts", False)

# -- graphical abstract: one dot per region, no averaged delta, no disease ranking, CBARP visible
_gs = io.open(f"{P}/figures/graphical_abstract.svg", encoding="utf-8").read()
_gcode = io.open(f"{P}/code/fig_graphical_abstract.py", encoding="utf-8").read()
close("graphical abstract draws one dot per region", 1, int("one dot per region" in _gs), tol=0)
close("graphical abstract code does not average regional deltas", 0, int("cliffs_delta.mean()" in _gcode), tol=0)
close("graphical abstract shows the CBARP junction", 1, int("CBARP splicing" in _gs), tol=0)
close("graphical abstract makes no correlation-absence claim", 0, int("does not correlate" in _gs), tol=0)
close("graphical abstract states that the Fura-2 data are not paired with RT-qPCR", 1, int("not paired with RT-qPCR" in _gs), tol=0)

# -- journal graphical abstract: Elsevier format (5:2, at least 1328 x 531 px at 300 dpi), large type, message at the evidence level
from PIL import Image as _PILImage
_gj = f"{P}/figures/graphical_abstract_journal"
for _e in ("png", "tiff", "pdf", "svg"):
    close(f"journal graphical abstract exists as .{_e}", 1, int(os.path.exists(f"{_gj}.{_e}")), tol=0)
_ji = _PILImage.open(f"{_gj}.tiff")
close("journal graphical abstract is at least 1328 x 531 px", 1, int(_ji.size[0] >= 1328 and _ji.size[1] >= 531), tol=0)
close("journal graphical abstract has the 5:2 ratio", 2.5, _ji.size[0] / _ji.size[1], tol=0.001)
close("journal graphical abstract TIFF is 300 dpi RGB", 1, int(round(_ji.info["dpi"][0]) >= 300 and _ji.mode == "RGB"), tol=0)
_js = io.open(f"{_gj}.svg", encoding="utf-8").read()
for _s in ["RNA candidates for replicated follow-up, not a mechanism", "in 5 of 6 datasets", "not paired with mRNA", "one dot per region"]:
    close(f"journal graphical abstract says '{_s}'", 1, int(_s in _js), tol=0)
for _s in ["does not correlate", "Graphical abstract", "graphical abstract", "high-confidence"]:
    close(f"journal graphical abstract does not contain '{_s}'", 0, int(_s in _js), tol=0)
_jcode = io.open(f"{P}/code/fig_graphical_abstract_journal.py", encoding="utf-8").read()
close("journal graphical abstract does not average regional deltas", 0, int("cliffs_delta.mean()" in _jcode), tol=0)
_readme_txt = io.open(f"{P}/README.md", encoding="utf-8").read()
close("README names the journal graphical abstract", 1, int("figures/graphical_abstract_journal" in _readme_txt), tol=0)

# -- highlights
close("highlights: four items", 4, len(_highlights_doc.paragraphs) - 1, tol=0)
close("highlights: every item is at most 85 characters", 1, int(all(len(q.text) <= 85 for q in _highlights_doc.paragraphs[1:])), tol=0)

# -- figures: every picture has alt text; Supplementary Figure S9 exists and is cited
for _d, _lab in ((_main_doc, "manuscript"), (_supp_doc, "supplement")):
    _alts = [dp.get("descr") or "" for dp in _d._element.xpath(".//wp:docPr")]
    close(f"{_lab}: every picture has alt text ({len(_alts)} pictures)", 0, sum(len(a) < 20 for a in _alts), tol=0)
close("main manuscript embeds six figures", 6, len(_main_doc._element.xpath(".//wp:docPr")), tol=0)
close("supplement embeds nine figures", 9, len(_supp_doc._element.xpath(".//wp:docPr")), tol=0)
close("Supplementary Figure S9 is cited in the manuscript", 1, int("Supplementary Figure S9" in SEARCH_TXT), tol=0)
close("supplementary contents list the nine figures and six result sections", 1,
      int("S1–S9" in _supp_doc.paragraphs[8].text and "1–6" in _supp_doc.paragraphs[7].text), tol=0)
close("Supplementary Figure S5 alt text no longer calls the conditions units of inference", 0,
      int(any("units of inference" in (dp.get("descr") or "") for dp in _supp_doc._element.xpath(".//wp:docPr"))), tol=0)

# -- repository documents agree with the manuscript files
_leg = io.open(f"{P}/supplementary/Supplementary_Figure_Legends.md", encoding="utf-8").read()
close("stand-alone legends cover Supplementary Figures S1-S9", 9, len(re.findall(r"\*\*Supplementary Figure S\d\.\*\*", _leg)), tol=0)
for _t in ["units of inference", "high-confidence", "Replicate-level PSI for SOCE-related splicing candidates"]:
    close(f"stand-alone legends do not contain '{_t}'", 0, int(_t in _leg), tol=0)
_ld = _DocxDocument(f"{P}/supplementary/Supplementary_Figure_Legends.docx")
close("legend .docx and .md hold the same nine legends", 1,
      int(sum(bool(re.match(r"Supplementary Figure S\d\.", p.text)) for p in _ld.paragraphs) == 9), tol=0)
_rd = io.open(f"{P}/README.md", encoding="utf-8").read()
for _t in ["Figures 6-8", "Figures 6–8", "fig_supp_splicing.py", "eight\nsupplementary figures", "Supplementary Figures S1–S8",
           "high-confidence", "TDP-43 knockdown is associated with a lower"]:
    close(f"README does not contain the stale '{_t}'", 0, int(_t in _rd), tol=0)
for _t in ["Supplementary Figures S1–S9", "v1.0.8", "Which counts feed which result", "Evidence levels"]:
    close(f"README contains '{_t}'", 1, int(_t in _rd), tol=0)
for _f in ["CORRECTIONS_2026-09-29_v1.0.8.md", "source_data/power_simulation_S1.csv", "source_data/DESeq2_ctrl_vs_75_featureCounts.csv",
           "code/build_supplementary_legends.py"]:
    close(f"file named in the README exists: {_f}", 1, int(os.path.exists(f"{P}/{_f}")), tol=0)

# -- inventories
for _s in ["S16b. Donor-level re-analysis of multiple sclerosis", "The supplementary material also contains Figures S1–S9"]:
    check("present", _s, True)
close("supplementary document inventories S16b", 1, int("S16b. Donor-level re-analysis of multiple sclerosis" in _supp2), tol=0)
_main_len = len(re.sub(r"\s+", " ", _body.split("## 1. Introduction")[1].split("## Declarations")[0]).split())
out.append(f"    Introduction to Conclusion, including captions and Limitations: {_main_len} words")
close("Introduction to Conclusion is shorter than in v1.0.7 (13,396 words)", 1, int(_main_len < 13396), tol=0)

out.append(f"\n==== {ok} passed, {bad} failed ====")
os.makedirs(f"{P}/logs", exist_ok=True)
io.open(f"{P}/logs/consistency_check_reviewed.txt", "w", encoding="utf-8").write("\n".join(out))
print("\n".join(out[-3:]))
print("report:", f"{P}/logs/consistency_check_reviewed.txt")
if bad:
    print("\n".join([l for l in out if l.startswith("FAIL")]))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v3 — automated consistency check between MANUSCRIPT_v4_SUBMISSION.md and the package files.

Every check re-reads the number from the file that produced it and compares it with the
string in the manuscript. Writes 09_YAYIN_PAKETI/logs/consistency_check.txt
"""
import io, os, re
import numpy as np
import pandas as pd

M = "/Users/elmas/Desktop/MAKALE"
P = f"{M}/09_YAYIN_PAKETI"
TXT = io.open(f"{P}/manuscript/MANUSCRIPT_v4_SUBMISSION.md", encoding="utf-8").read()
out, ok, bad = [], 0, 0


def check(label, snippet, expected=True):
    global ok, bad
    present = snippet in TXT
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


out.append("=== Table 1: coverage pre-filter ===")
t1 = pd.read_csv(f"{P}/tables/Table1_coverage_prefilter.csv")
close("events removed range low", 18, round(t1.events_removed_pct.min()), tol=0.6)
close("events removed range high", 78, round(t1.events_removed_pct.max()), tol=0.6)
close("significant lost low", 33, round(t1.significant_lost_pct.min()), tol=0.6)
close("significant lost high", 76, round(t1.significant_lost_pct.max()), tol=0.6)
close("SH-SY5Y significant before filter", 7854,
      int(t1.loc[t1.dataset.str.startswith("SH-SY5Y"), "significant_before_filter"].iloc[0]), tol=0)

out.append("\n=== Table 2: robust SOCE events ===")
t2 = pd.read_csv(f"{P}/tables/Table2_robust_SOCE_splicing_events.csv").set_index("gene")
for g, d, lo, hi in [("STIMATE", 0.244, 0.095, 0.368), ("ORAI3", -0.269, -0.404, -0.107),
                     ("STIM2", -0.120, -0.165, -0.064), ("STIM1", 0.145, -0.002, 0.293)]:
    close(f"{g} delta_PSI", d, t2.loc[g, "delta_PSI"])
    close(f"{g} CI low", lo, t2.loc[g, "CI95_low"])
    close(f"{g} CI high", hi, t2.loc[g, "CI95_high"])
close("STIM1 exon start 1-based", 4088702, int(t2.loc["STIM1", "start_1based"]), tol=0)
close("STIM1 exon length", 37, int(t2.loc["STIM1", "exon_bp"]), tol=0)

out.append("\n=== Table 3 / S14: cryptic calls and null test ===")
t3 = pd.read_csv(f"{P}/tables/Table3_cryptic_events_eleven_comparisons.csv").set_index("comparison")
close("SH-SY5Y high-confidence events (regtools)", 165,
      int(t3.loc["SH-SY5Y 75 ng/mL", "high_confidence_events"]), tol=0)
close("SH-SY5Y genes", 113, int(t3.loc["SH-SY5Y 75 ng/mL", "genes"]), tol=0)
close("iPSC-MN TDP-43 high-confidence events", 18, int(t3.loc["iPSC-MN, TDP-43 KD", "high_confidence_events"]), tol=0)
close("iPSC-MN FUS high-confidence events", 26, int(t3.loc["iPSC-MN, FUS KD", "high_confidence_events"]), tol=0)
close("iPSC-MN TAF15 high-confidence events", 23, int(t3.loc["iPSC-MN, TAF15 KD", "high_confidence_events"]), tol=0)
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
close("burden minimum", 12, int(t3.high_confidence_events.min()), tol=0)
_mouse = t3.loc[["C2C12", "NSC34", "Mouse striatum"]]
close("Table 3 mouse rows: positive controls not assessed", 3,
      int((_mouse.positive_controls_high_confidence == "n/a (mouse)").sum()), tol=0)
close("burden maximum", 477, int(t3.high_confidence_events.max()), tol=0)
s14 = pd.read_csv(f"{P}/supplementary/S14_control_vs_control_null_test.csv")
_s14p = s14[s14.tier == "permissive"].set_index("dataset")
close("permissive null ratio iPSC colonies (Methods 2.5: 0.98 per real call)", 0.98,
      _s14p.loc["iPSC colonies", "null_to_real_ratio"])
close("permissive null ratio K562 total RNA (Methods 2.5: 2.01)", 2.01,
      _s14p.loc["K562 total RNA", "null_to_real_ratio"])
s14 = s14[s14.tier == "tier 2 (high-confidence)"].set_index("dataset")
close("null ratio iPSC colonies", 0.64, s14.loc["iPSC colonies", "null_to_real_ratio"])
close("null ratio K562 total RNA", 2.17, s14.loc["K562 total RNA", "null_to_real_ratio"])
close("null ratio mouse striatum", 0.83, s14.loc["Mouse striatum", "null_to_real_ratio"])

out.append("\n=== Table 4: transcript-family abundance (TPM) ===")
t4 = pd.read_csv(f"{P}/tables/Table4_transcript_family_abundance.csv")
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
             ("STIM (ER Ca2+ sensor)", 48), ("ORAI (CRAC channel)", 27), ("SOCE regulators", 30)]:
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
_t2 = pd.read_csv(f"{P}/tables/Table2_robust_SOCE_splicing_events.csv").set_index("gene")
close("STIM2 intron13 5' window contains the Table 2 STIM2 exon", 1,
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
close("WST-1 48 h knockdown", 61.5, x.loc["WST-1 viability 48 h | shTDP-43", "mean"], tol=0.05)
close("WST-1 24 h knockdown", 116.8, x.loc["WST-1 viability 24 h | shTDP-43", "mean"], tol=0.05)

out.append("\n=== cross-file agreement and the corrected correction family ===")
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
for _f in ["S5_high_confidence_cryptic_events", "S6_cryptic_positive_controls",
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
_t3 = pd.read_csv(f"{P}/tables/Table3_cryptic_events_eleven_comparisons.csv")
close("Table 3 rows = eleven comparisons", 11, len(_t3), tol=0)
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
_unnamed = sorted({g for g in _ex.gene if f"*{g}*" not in TXT.split("## References")[0]})
close("Ca2+ units excluding zero that Section 3.5 does not name", 0, len(_unnamed), tol=0)
out.append("    " + ("all named: " + ", ".join(sorted(set(_ex.gene)))) if not _unnamed
           else "    unnamed: " + ", ".join(_unnamed))
_s10b = pd.read_csv(f"{P}/supplementary/S10b_NMD_panel_level_tests.csv").set_index("panel")
close("S10b positive-control panel p", 0.62,
      float(_s10b.loc["Cryptic_positive_controls_16", "p_one_sided_MWU"]), tol=0.005)

out.append("\n=== manuscript strings that must be present ===")
for s in ["10,926 versus 176 reads", "three spinal cord levels", "six of the seven brain regions",
          "0.64 calls in iPSC colonies, 2.17 in K562 total RNA and 0.83 in mouse striatum",
          "interaction was −0.325 (p = 0.317, q = 0.508)", "chr11:4,088,702–4,088,738",
          "chr4:27,007,983–27,008,006", "three independent cultures", "10 µM cyclopiazonic acid",
          "Albarran L, Lopez JJ, Woodard GE, Salido GM, Rosado JA",
          "Of the 110 correlations tested against the junction-based marker",
          "116.8 ± 1.6% of control at 24 h", "Cutadapt v5.2",
          "59 units in 29 genes passed the same depth filter",
          "from 0.570 in controls to 0.819 in knockdown (Δ = +0.249",
          "`-p --countReadPairs` for paired-end libraries",
          "74 qualifying units in C2C12 and 131 in NSC34",
          "*SARAF* intron 5 index rises in the iPSC-derived motor neurons (+0.097) and the *Saraf* intron 5 index in NSC34 (+0.204)",
          "almost as many calls as the real comparison in iPSC colonies (0.98 per real call)",
          "this check is available only for the human comparisons",
          "isoform-level testing did not detect a *CBARP* isoform switch",
          "the three mouse comparisons have no conserved control",
          "*ORAI3* went from 8% to 30% of ORAI transcripts and *ORAI2* from 77% to 53%",
          "The STIM (+48%) and SOCE-regulator (+30%) pools also rose, whereas the SERCA and mitochondrial-uptake pools fell (−12% and −25%)",
          "the median gene expressed in both groups (mean TPM > 5) had 24% lower TPM",
          "The *STIM2* unit is not an independent observation",
          "Yoast RE, Emrich SM, Zhang X, et al.",
          "3\u2076 = 729 combinations for the three-versus-three comparisons and 2\u2074 = 16",
          "applied across the 220 informative correlations",
          "Viability was measured at 24 h and 48 h",
          "10 MS/5 control donors, 98 samples",
          "in SH-SY5Y the same unit is uninformative (0.000, interval \u22120.264 to +0.241)",
          "decreased at donor level across all sampled lesion types (\u03b4 = \u22120.840; q = 0.038)"]:
    check("present", s, True)
out.append("\n=== strings that must be gone ===")
for s in ["Sah P, et al.", "10,930", "four spinal cord regions", "p = 0.13)", "log2FC = −1.746",
          "single culture plate", "not independent biological replicates",
          "prioritized candidate", "full GENCODE v47 index", "Cutadapt v4.6",
          "the mouse datasets were not analysed at all",
          "all 3\u2076 = 729 replicate combinations", "their p values are therefore optimistic",
          "applied across all 231 correlations", "is not measurable in SH-SY5Y",
          "used in the first version of this analysis", "The original analysis treated eight contrasts",
          "the correction was applied only to the primary SH-SY5Y model",
          # 22 September 2026 audit
          "require independent biological replication", "Ca²⁺ imaging",
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
_body = TXT.split("## References")[0]
_first = {}
for _i in range(1, 10):
    _m = re.search(rf"Figure {_i}(?![0-9])", _body)
    close(f"Figure {_i} cited in the text", 1, int(_m is not None), tol=0)
    _first[_i] = _m.start() if _m else -1
close("figures first cited in numerical order", 1,
      int(all(_first[i] < _first[i + 1] for i in range(1, 9))), tol=0)
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

out.append("\n=== Tables 1-5 in the manuscript match tables/*.csv ===")
import subprocess as _sp
_r = _sp.run(["/usr/bin/python3", f"{P}/code/build_manuscript_docx.py", "--check"],
             capture_output=True, text=True)
close("Markdown tables regenerated from the CSV files are identical", 0, _r.returncode, tol=0)
for _n in range(1, 6):
    check(f"table {_n} block present", f"<!-- table:{_n} -->", True)

out.append(f"\n==== {ok} passed, {bad} failed ====")
os.makedirs(f"{P}/logs", exist_ok=True)
io.open(f"{P}/logs/consistency_check.txt", "w", encoding="utf-8").write("\n".join(out))
print("\n".join(out[-3:]))
print("report:", f"{P}/logs/consistency_check.txt")
if bad:
    print("\n".join([l for l in out if l.startswith("FAIL")]))

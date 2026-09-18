# Submission package — TDP-43 knockdown and store-operated Ca²⁺ entry

**Built:** 18 September 2026 · **Manuscript version:** v4 (submission draft)

This repository holds the material for the manuscript. The analyses of the underlying doctoral thesis live in a separate repository, [tdp43-thesis-reproducibility](https://github.com/elmasnuryilmaz/tdp43-thesis-reproducibility). Several thesis-era results were superseded during preparation of this manuscript and should not be quoted from there: the alternative-polyadenylation coverage was recomputed after excluding CIGAR reference skips (`samtools bedcov -j`), the NMD interaction was re-tested with the four intervention conditions as the unit of inference instead of eight dependent contrasts, the Salmon gene-level summaries were rebuilt with the complete GENCODE v47 transcript-to-gene map and reported as TPM, the multiple sclerosis tests were repeated with the donor as the unit of inference, and the *TRPC1* exon-skipping event is no longer presented as a finding. Where the two repositories disagree, this one is current.

Everything a journal or a reviewer needs is in this folder. Every number in the manuscript
was regenerated from the files here by the scripts in `code/`, and the consistency between
the two is checked automatically (`logs/qa_manuscript_vs_data.txt`, 94 checks, all passing).

```
09_YAYIN_PAKETI/
├── manuscript/      MANUSCRIPT_v4_SUBMISSION.md and .docx
├── figures/         Figure 1–9, PNG (300 dpi) and PDF, numbered as in the manuscript
├── tables/          Table 1–5, CSV, English headers
├── supplementary/   S1–S17
├── source_data/     recomputed intermediate data behind Tables 4–5 and Figures 6–8
├── code/            every script, including the v3 rebuild scripts
├── logs/            automated QA report
└── qa/              internal working files (thesis text dump) — NOT for submission
```

## Figures

| File | Manuscript | Content |
|---|---|---|
| `Figure1_TRPC1_robustness` | Figure 1 | replicate-level failure of the TRPC1 call |
| `Figure2_detection_power` | Figure 2 | power in a 3 + 3 design |
| `Figure3_robust_candidates` | Figure 3 | SOCE splicing candidates with bootstrap CIs |
| `Figure4_cryptic_discovery_and_NMD` | Figure 4 | cryptic discovery, Ca²⁺ panels, NMD interaction |
| `Figure5_specificity_and_APA` | Figure 5 | RBP specificity, corrected APA, STMN2 intron 2 coverage |
| `Figure6_functional_consequences` | Figure 6 | laboratory experiments (new in v4) |
| `Figure7_transcript_family_abundance` | Figure 7 | family abundance, TPM based (recomputed in v4) |
| `Figure8_TRPC1_disease_direction` | Figure 8 | TRPC1 across five diseases (all regions, donor-level MS) |
| `Figure9_NYGC_cryptic_STMN2` | Figure 9 | junction-level cryptic STMN2 in ALS tissue |

Figure numbers are **not** burned into the images; the file name carries the number.

## Supplementary files

| File | Content |
|---|---|
| `S1_laboratory_source_data.xlsx` | raw Ct, relative expression, Fura-2 and WST-1 values per replicate, primers, thermal profile, summary statistics |
| `S2_calcium_gene_panels.csv` | the four cumulative Ca²⁺ panels |
| `S3_rMATS_significant_events.csv.gz` | every rMATS event at FDR < 0.05 and \|ΔPSI\| ≥ 0.10, JC and JCEC, six datasets, with raw junction counts |
| `S4_matched_permutation_enrichment.csv` | enrichment against the covariate-matched null |
| `S5_high_confidence_cryptic_events.csv` | high-confidence unannotated splicing changes, all comparisons |
| `S6_cryptic_positive_controls.csv`, `S6b_positive_control_matrix.csv` | positive-control recovery (matrix headers rebuilt in v4) |
| `S7_SOCE_genes_annotation_free.csv` | SOCE genes in the junction-level analysis |
| `S8_cryptic_STMN2_ALS_vs_control.csv` | cryptic PSI by region, ALS versus control |
| `S9_cryptic_PSI_correlations_within_ALS.csv` | all 231 correlations within ALS samples |
| `S10_NMD_interaction_SOCE_panel.csv`, `S10b_NMD_panel_level_tests.csv` | four-condition NMD interaction and panel-level tests (recomputed in v4) |
| `S11_APA_candidate_gradients.csv` | depth-qualified coverage gradients after the `-j` correction |
| `S12_cryptic_counts_by_dataset.csv` | call counts per comparison with the RBP controls |
| `S13_STIM2_SOAR_exon_junction_level.csv` | SOAR exon at junction level |
| `S14_control_vs_control_null_test.csv` | split-control null test, three tiers |
| `S15_STIM2.1_exon_six_datasets.csv` | STIM2.1 meta-analysis |
| `S16_multiple_sclerosis_both_cohorts.csv`, `S16b_multiple_sclerosis_donor_level.csv` | MS analysis; the donor-level re-analysis is new in v4 |
| `S17_dataset_accessions.csv` | every accession with design, library type and sample-to-group assignment |

The invalid legacy table "cryptic-junction-positive and NMD-sensitive genes" was **removed**:
it listed genes with the superseded eight-contrast q-values and contradicted the corrected
analysis, in which no gene passes genome-wide FDR.

## How to reproduce

```bash
/usr/bin/python3 code/fig_v3_main.py        # Figures 1, 2, 3, 7
/usr/bin/python3 code/fig_v3_junction.py    # Figures 4, 5, 9
/usr/bin/python3 code/fig_v3_lab.py         # Figure 6
/usr/bin/python3 code/fig_v3_disease.py     # Figure 8
/usr/bin/python3 code/build_tables_v3.py    # Tables 1-5, S2, S4-S16
/usr/bin/python3 code/build_source_data.py  # S1
/usr/bin/python3 code/build_S3_rmats.py     # S3   (needs the rMATS output)
/usr/bin/python3 code/build_S17_datasets.py # S17
/usr/local/bin/Rscript  code/deseq_full.R   # DESeq2 on the complete transcript map
/usr/bin/python3 code/ms_donor_level.py     # donor-level MS analysis
/usr/bin/python3 code/qa_check_v4.py        # consistency check
```

Paths are absolute inside the scripts and point at `~/Desktop/MAKALE` and `~/Desktop/TEZ`;
they must be adapted before the package is published.

## Known gaps to close before submission

1. **Funding, author contributions and acknowledgements** are placeholders in the
   Declarations section of the manuscript.
2. **GitHub and Zenodo links** in Section 2.17 are placeholders.
3. ~~CPA concentration~~ — **resolved 18 September 2026: 10 µM is correct** (author
   confirmation). The Prism trace export was annotated 10⁻³ M in error; the relabelled
   version used in Figure 6D carries the correct value. The original export is kept as
   `source_data/representative_Fura2_traces_from_thesis.png` for provenance.
4. ~~Representative Fura-2 traces~~ — **done**: relabelled in English with the correct CPA
   concentration (`source_data/representative_Fura2_traces_relabelled.png`, produced by
   `code/relabel_traces.py`) and added as Figure 6D. Only the text annotations were
   replaced; the traces and axes are untouched.
5. ~~Alternative polyadenylation in the iPSC-derived motor neurons~~ — **done 18 September
   2026** once the external drive was mounted (`code/run_bedcov_v3_remaining.sh` →
   `code/recalc_apa_remaining_R.R`; results in `S11_APA_candidate_gradients.csv`). The
   *STMN2* positive control works in that model (+0.249), which is what the SH-SY5Y analysis
   could not show. The two mouse lines were recomputed as well (C2C12 and NSC34, 74 and 131
   qualifying units); all four comparisons are in `S11_APA_candidate_gradients.csv`.
6. **Abstract** is 404 words; trim to the target journal's limit.
7. The `qa/` folder holds internal working files and should not be uploaded.

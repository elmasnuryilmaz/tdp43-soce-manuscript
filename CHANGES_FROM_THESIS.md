# Analysis revisions relative to the doctoral thesis

The analyses behind this manuscript began as part of the first author's doctoral thesis (Ege
University, 2026). Five analyses were revised during preparation of the manuscript. This
document records what changed and why, so that numbers in the thesis and in this repository
can be told apart. Where they differ, the values in this repository are the current ones.

## 1. Alternative polyadenylation: reference skips excluded from coverage

`samtools bedcov` includes CIGAR reference skips (`N`) and deletions in its depth sum by
default, so an intron that a spliced read never traverses still accumulates coverage. A
synthetic read with CIGAR `10M1000N10M` returns a depth of 100 over a window lying entirely in
the `N` block without `-j`, and 0 with it. All coverage indices were recomputed with
`samtools bedcov -Q 30 -j`.

Effect: the *STMN2* intron 2 positive control falls from +0.489 to +0.145 in SH-SY5Y, and the
control then lies just below the pre-specified depth filter in that model. The corrected
analysis was extended to every comparison with available alignments (SH-SY5Y, iPSC-derived
motor neurons, C2C12, NSC34); the control behaves as intended in the iPSC-derived motor neurons
(+0.249). The mouse lines have no equivalent control, because the *STMN2* cryptic
polyadenylation site is absent from the mouse gene. Results:
`supplementary/S11_APA_candidate_gradients.csv`, which now lists the genomic windows of every unit.

## 2. NMD interaction: four conditions as the unit of inference

The first analysis treated eight difference-of-differences values as independent observations,
although the four NMD-inhibition conditions reuse the same two control and two TDP-43 knockdown
libraries. The test was repeated with the four conditions as the unit of inference (`df = 3`),
averaging replicates within each condition.

Effect: no gene passes genome-wide FDR correction (smallest q = 0.053). *CBARP* moves from
q = 0.0049 to q = 0.145 and is reported as a prioritised candidate rather than an established
NMD target; the predicted NMD sensitivity of *STIM1* is withdrawn (interaction −0.325,
p = 0.317). Results: `supplementary/S10_NMD_interaction_SOCE_panel.csv`,
`supplementary/S10b_NMD_panel_level_tests.csv`.

## 3. Transcript quantification: complete transcript-to-gene map, TPM

The Salmon index was built from the comprehensive GENCODE v47 transcript set (384,354
transcripts), but gene-level summaries had been aggregated with the transcript-to-gene map of
the basic annotation, which covers 157,588 of them and discards 26% of the transcripts per
million. Counts were re-aggregated with the complete map, DESeq2 was re-run, and family
abundance is now reported as length-corrected TPM rather than counts per million. Because a
few abundant transcripts take a larger share of the TPM total in the knockdown libraries (the
median expressed gene has 24% lower TPM there), changes between conditions are computed after a
median-of-ratios adjustment for library composition (`code/build_family_abundance.py`); with
that adjustment the ORAI pool rises by 27% and shifts towards ORAI3, rather than staying flat.

Effect on the genes discussed in the manuscript: *STIM1* log2FC 0.803 → 0.929, *TRPC1*
0.967 → 0.958, *ATP2A3* 1.202 → 1.306, *CBARP* −1.746 → −1.254, and *ORAI1* becomes
significant (0.364, not significant → 0.433, q = 2.4 × 10⁻⁴), which strengthens agreement with
the RT-qPCR measurement. Results: `tables/Table1_transcript_family_abundance.csv`,
`source_data/DESeq2_ctrl_vs_75_fullmap.csv`.

## 4. Multiple sclerosis cohort: donor as the unit of inference

Several samples in GSE138614 come from the same donor and had been treated as independent. Each
test was repeated after averaging within donors.

Effect: in normal-appearing white matter the reduction of *TRPC1* is δ = −0.482 at sample
level (p = 0.005, q = 0.10 after Benjamini–Hochberg correction) and δ = −0.771 at donor level
(uncorrected p = 0.030, seven multiple sclerosis donors versus five controls). The myelin- and
glia-adjusted comparison stays significant at sample level (δ = −0.418, p = 0.002) but not at
donor level (δ = −0.640, p = 0.055). The manuscript reports both levels and treats the donor
level as primary. Results: `supplementary/S16b_multiple_sclerosis_donor_level.csv`.

## 5. The *TRPC1* exon-skipping event is withdrawn as a finding

The event (chr3:142,792,824–142,792,967; ΔPSI = +0.108, FDR = 0.0459) rests on seven skipping
reads across six libraries, with the skipping form absent from four of them. Its bootstrap
confidence interval spans zero (−0.093 to +0.522), removing one control replicate reverses the
sign, it does not survive coverage pre-filtering, and LeafCutter does not call it. The
manuscript reports it as an illustration of why threshold-based splicing calls need read-level
verification, not as a result. Data: `figures/Figure2_TRPC1_robustness.*`,
`supplementary/S3_rMATS_significant_events.csv.gz`.

## Software versions

The versions stated in the manuscript were checked against the records of the run itself. One
correction followed: the trimming reports show **Cutadapt 5.2**, not the 4.6 recorded in the
thesis. HISAT2 2.2.2, SAMtools 1.23 (alignment) and 1.21 (junction and coverage steps),
featureCounts 2.1.1, Salmon 1.11.4, rMATS 4.3.0, LeafCutter 0.2.9, SUPPA2 2.3, regtools 1.0.0
and FRASER 1.14.1 were confirmed. Full list: `environment.yml`.

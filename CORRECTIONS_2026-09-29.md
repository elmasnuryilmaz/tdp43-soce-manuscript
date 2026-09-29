# Scientific corrections for release v1.0.6

The manuscript, supplementary document, source workbook and current figures were
revised together after checking the analysis outputs and clarifying the experimental
units with the author. This release supersedes the scientific descriptions in v1.0.5.

## Experimental units and laboratory interpretation

- The author confirmed that RT-qPCR n = 4 represents technical measurements of the
  same biological sample per group. Fura-2 n = 3 represents wells on one culture plate.
  The available WST-1 data comprise four wells per group from one experiment.
- Laboratory measurements are descriptive throughout the abstract, methods, results,
  discussion, captions, highlights, graphical abstract and source workbook. Inferential
  qPCR p values and significance stars were withdrawn. Public RNA-seq sample designs
  are unaffected. Technical SEM does not estimate biological variability.
- Original Ct, relative-expression, Fura-2 and WST-1 measurement cells are preserved.
  Technical-repeat labels and summaries were corrected. The remaining two WST-1
  experiments are unavailable; no measurements were imputed or invented.
- The calcium-readdition signal and its ratio to release do not isolate membrane
  influx or eliminate store-related effects. Mechanistic interpretations are framed as
  hypotheses. Acute SERCA reuptake is not offered as an explanation while CPA is present.

## RNA-processing analyses

- S13: fixed the upstream-exon boundary convention and mouse chromosome-prefix
  mismatch. The corrected table has 16 available junction records, eight upstream and
  eight downstream, across seven TDP-43 comparisons and the FUS/TAF15 controls.
  K562 total RNA has only a downstream record; FUS has only an upstream record.
  Both mouse cell models are now included. C2C12 has negative estimates on both
  flanks. Human downstream PSI uses junctions sharing an acceptor. The independently
  computed rMATS meta-analysis is unchanged.
- STIM1: nine DRIMSeq features comprise seven evaluable stageR transcripts (all
  confirmation-adjusted p = 1) and two missing feature p values. ENST00000698912.1
  and ENST00000698913.1 were not evaluated by stageR; these are the two significant
  features in the separate IsoformSwitchAnalyzeR/DEXSeq analysis. Neither is a predicted
  PTC isoform. The text no longer equates missing confirmation with a failed test or
  claims that an NMD substrate is absent. All nine eligibility records are published.
- NMD: the four contrasts reuse baseline libraries. Earlier t/sign tests did not include
  shared-baseline uncertainty; their p/q values and panel tests are withdrawn.
  S10, S10b and Figure S5 now describe estimates, ranges and panel medians. Full
  descriptive results for 19,145 genes are provided. CBARP is +1.52 log2 (four positive
  contrasts, range +1.19 to +2.23), which does not establish an NMD target.
- Cryptic-splicing reference genes are not a validated NMD-positive panel. Their
  negative group test cannot establish assay failure or make negative SOCE results
  conclusive. The earlier corresponding statements were removed.
- SARAF intron labels across species were ordinal matches, not sequence alignment;
  they identify gene-level candidates, not a replicated homologous APA event.
- FRASER and SUPPA2 statements were qualified to distinguish a limited or discordant
  result from proof of expected method behaviour or calibration. FUS/TAF15 specificity
  is not established by the small matched comparison.

## Expression, tissue comparisons and statistical interpretation

- The svaseq account now separates q < 0.05 from the manuscript's joint DE threshold
  q < 0.05 and absolute log2FC >= 1. Five named genes retain q < 0.05; only ORAI3 and
  CBARP meet both thresholds in both fits. The median 0.20 difference is absolute.
- NYGC, Alzheimer, Parkinson and multiple-sclerosis filters and normalizations now
  match their actual analysis code. Sample-level and donor-level units are distinguished.
- New S18d averages expression within donor and region before comparison/adjustment.
  ALS cerebellum has 158 samples from 147 donors. Its TRPC1 donor-level difference
  remains positive (delta +0.532, q = 2.27e-6); marker-adjusted differences are +0.248
  to +0.370 and remain significant in all four models. The multiple-testing families
  for S18d differ from Table 5 and are explicitly described.
- Attenuation to nonsignificance is not described as proof of no association, absence
  of a mechanism, or proof that cell loss explains an association. Nominal panel
  enrichment and MS sensitivity results retain their uncertainty.
- Unknown diagnoses in the Other Neurological Disorders group remain a limitation;
  STMN2 cryptic activity is not used to assign an unobserved diagnosis.

## Document and package integrity

S10b and S18/S18b/S18c/S18d are listed consistently. Figure captions, gene formatting,
normalization descriptions and several ambiguous sentences were corrected. The AI-use
statement reflects assistance with analysis code and figure preparation as well as
language. Zotero fields are retained. Old scripts remain historical where documented;
current analysis and figure entrypoints are listed in README.md.

Verification: the existing consistency suite plus the September 29 audit suite check
source-linked quantities, laboratory-cell preservation and the corrected interpretation.
The DOCX documents were rendered and inspected visually before packaging. Passing
these checks does not substitute for independent biological replication.

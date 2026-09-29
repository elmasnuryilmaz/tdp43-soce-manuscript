# Revisions for release v1.0.8

This release answers a second-round evaluation of release v1.0.7 (manuscript, supplement, highlights,
graphical abstract, figures, source files and README). The evaluation found no error in the measurements
themselves; its findings concern how far the text and figures claim more than the evidence supports, names that
suggested validation, and files that had drifted apart. No measured value, no RNA-seq result and no table entry
was changed; wording, labels, figures, one Methods paragraph and the repository documents were.

## 1. The claim is set to the evidence level

- The Fura-2 result is presented as a single-plate observation (three wells per group, one culture plate),
  descriptive, not paired with the RT-qPCR or RNA-seq samples, and not replicated. The reanalysis of public
  RNA-seq data is the main content of the article. Title, abstract, last paragraph of the Introduction,
  Results 3.1, Discussion, Conclusion, highlights and graphical abstract were rewritten to say so.
- The title is now *Calcium-regulatory RNA candidates after TDP-43 knockdown in SH-SY5Y cells: reanalysis of
  public RNA-seq data and a single-plate Ca²⁺-readdition observation*.
- "Discordance" between RT-qPCR and Fura-2 is described as a difference in direction between two unpaired
  experiments. The RT-qPCR and Fura-2 measurements are said to come from different cultures and time points.
- The Discussion follows the order: RNA findings, dependence on model and method, the *CBARP* candidate,
  possible contributors to the single-plate difference (none tested), tissue, and "Working hypotheses and
  discriminating experiments". The mechanism paragraph was condensed and the Limitations were shortened.
- Length, counted the same way for both releases (Introduction to Declarations, with captions and Limitations):
  13,460 words in v1.0.7, 13,128 now. Results are 269 words shorter, the Discussion 255 and the Limitations 76;
  Methods are 60 words longer (the simulation assumptions and the count sources), the Introduction 28 and the
  Conclusion 14.

## 2. Names that implied validation

- "High-confidence" unannotated splicing changes are now "stringent-filter unannotated splicing candidates":
  the control-versus-control null gave null-to-real call ratios of 0.64, 2.17 and 0.83, so the filter does not
  establish specificity. Renamed in the text, in Table 4, S5 (file renamed with `git mv`), S12 and S14, in the
  generating code and in the figure labels. Column names changed; the numbers did not.
- "SOCE machinery", "SOCE regulators" and "core entry components" are replaced by the fixed lists they stood for:
  the nineteen-gene SOCE-associated set, the twelve core SOCE-pathway genes, and the family "Ca²⁺-entry regulators".
- Wording that read as conclusion or as conversation was replaced ("did most of the work", "independently
  supported" for tools applied to the same libraries, question-form headings, "Integrated mechanism").

## 3. Detection-power simulation

Methods 2.3 and the legend of Supplementary Figure S1 now state the assumptions of the simulation exactly as the
code runs them: PSI 0.5 ± ΔPSI/2 with a between-replicate standard deviation of 0.05, binomial reads at 10, 20, 50
and 100 informative reads per sample, a two-sample t test at nominal p < 0.05 without correction, 2,000
simulations per cell. It is an example, not an estimate of rMATS power. The plotted values are stored in
`source_data/power_simulation_S1.csv`; the reads per sample of the nominally significant SH-SY5Y events (lower
quartile 9.7, median 21.3) are given for context.

## 4. Figures

- **Figure 2**: the redundant panels are gone. A shows the two representative recordings on common axes, aligned
  to the steepest point of the readdition rise (no plotted value smoothed or rescaled); B–D show ER release,
  readdition and their ratio for the three wells per group. The two original recordings are Supplementary Figure S9.
- **Figure 4**: the exon length is in every panel title; the 119-bp *STIM2* exon of the skipped-exon events and the
  24-nt SOAR exon (STIM2.1) of panel B are shown as different exons; a schematic explains PSI.
- **Figure 5**: rows carry the numbers of cases and controls; the legend says the effects are unadjusted and that
  the correction is within region.
- **Figure 6**: redrawn to match its legend. The ER-release change is drawn, the RT-qPCR and Fura-2 boxes are
  marked as different cultures, and *CBARP* has no arrow to the readdition response.
- **Graphical abstract**: one dot per region instead of an averaged Cliff's δ, no disease ranking, the *CBARP*
  junction panel added, the Fura-2 panel labelled one plate and not paired with RT-qPCR; the caption no longer
  says that anything "does not correlate".
- **Supplementary Figures**: S1 (assumptions on the figure), S4 (panels use different thresholds), S6 (ordinal
  intron numbers, coverage gradients, control genes shaded), S7 (grouped by disease, multiple-sclerosis re-analyses
  separated), S8 (the detection threshold is written on the panel; Supplementary Table S18b uses another).
  Supplementary Figure S9 is new.
- Every picture in the manuscript and the supplement has alt text.

## 5. Structure

- NMD, the outlier screen (FRASER) and the coverage-based APA screen share Results 3.6 (they were 3.6 and 3.7).
  The detail moved to Supplementary Results 1–3. Sections 3.7 and 3.8 follow (were 3.8 and 3.9).
- svaseq numbers moved from Methods 2.2 to Results 3.2 and Supplementary Results 5 (a table recomputed from
  `source_data/svaseq_sensitivity_SHSY5Y.csv`). The nine *STIM1* transcript identifiers and the IsoformSwitchAnalyzeR
  results moved to Supplementary Results 6 (`source_data/STIM1_isoformswitch_DEXSeq.csv`).
- Supplementary Table S16b is listed in the inventories; the supplement and its contents list cover S1–S9.
- Highlights: four items (73, 77, 72 and 51 characters); the third now states the recurrent *CBARP* splice-junction finding.

## 6. Which counts give which DESeq2 result

Methods 2.2 had described featureCounts as the count source, whereas every DESeq2 statistic in the article comes from
the Salmon counts of the primary comparison (summed to genes with the complete transcript-to-gene map). Methods 2.2
now says so and reports the featureCounts run of the same six libraries as a cross-check: the log2 fold changes of
the 20,833 genes counted in both correlate with Pearson r = 0.85, and the four RT-qPCR targets have the same
direction. The README has the input-to-output table. New files: `source_data/DESeq2_ctrl_vs_75_featureCounts.csv`,
`source_data/DESeq2_featureCounts_vs_Salmon_calcium_genes.csv`, `code/compare_featurecounts_salmon.py`.
The summation of the Salmon `quant.sf` files to genes is not scripted in the repository and is stated as such.

## 7. Repository documents

- `supplementary/Supplementary_Figure_Legends.md` and `.docx` were written for an earlier layout (an old four-gene
  Supplementary Figure S3; the NMD conditions called "units of inference"). They are regenerated from the
  supplement by `code/build_supplementary_legends.py`.
- `README.md`: title, figure and supplementary tables, the evidence-level table, the removal of the stale
  "Figures 6–8" and `fig_supp_splicing.py` references, the count-source table, the rebuild order, and a single stated
  role for the graphical abstract (optional journal graphical abstract and repository summary, not a numbered figure).
- `DATA_AVAILABILITY.md` names release v1.0.8 and the new source files.
- The README sheet of `S1_laboratory_source_data.xlsx` (cell B5) no longer refers to a panel that does not exist and
  says that the label "SOCE" of the Fura2 sheet denotes the Ca²⁺-readdition amplitude. No measurement cell changed.

## 8. Checks

`code/qa_check_v4.py` compares the manuscript with the data files (terminology, title, power table, S3 quartiles, the
svaseq table, the *STIM1* identifiers, the count-source values, alt text, section numbering, the word count) and now
compares the numbers of the manuscript tables with the tables regenerated from the CSV files. `code/qa_audit_v118.py`
still confirms cell by cell that the laboratory measurements of the workbook are unchanged.

## 9. What these revisions do not change

The laboratory data are what they were: one Fura-2 plate, one WST-1 experiment of four wells (two further
experiments are unavailable), four biological RT-qPCR replicates with a single reference gene. Independent Fura-2
replicates cannot be produced from the existing files, and the RNA analyses remain candidates. The file checks are
consistency checks; they are not a statement that the science is sufficient for acceptance.

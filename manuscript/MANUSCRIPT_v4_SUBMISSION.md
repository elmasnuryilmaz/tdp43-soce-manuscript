# TDP-43 knockdown is associated with reduced store-operated Ca²⁺ entry and altered calcium-regulatory RNA profiles in SH-SY5Y cells

**Elmasnur Yılmaz¹, Yasemin Eraç¹\***

¹ Department of Pharmacology, Faculty of Pharmacy, Ege University, İzmir, Türkiye

\* Correspondence: yasemin.erac@ege.edu.tr

## Abstract

**Background.** Nuclear clearance of TDP-43 occurs in about 97% of amyotrophic lateral sclerosis (ALS) and half of frontotemporal lobar degeneration cases. Ca²⁺ handling is disturbed in these diseases, but whether TDP-43 loss affects the store-operated Ca²⁺ entry (SOCE) machinery, and at which level, is unresolved.

**Methods.** We reanalysed six public RNA-seq comparisons of TDP-43 depletion (38 human and mouse libraries), testing candidate splicing events with coverage filters, bootstrap intervals, leave-one-out checks and independent tools. Annotation-free junction analysis covered eleven comparisons, including FUS and TAF15 knockdown controls; alternative polyadenylation (APA) and nonsense-mediated decay (NMD) were explored. In SH-SY5Y cells with shRNA knockdown we measured target mRNAs by RT-qPCR, cytosolic Ca²⁺ by Fura-2/AM and viability by WST-1. Candidate genes were examined in ALS post-mortem and other neurological cohorts, using junction-level cryptic *STMN2* inclusion as a marker of TDP-43 dysfunction.

**Results.** *CBARP* was the most reproducible splicing candidate across species; *STIM2*, *STIMATE* and *ORAI3* were supported in the primary model. After knockdown, *TRPC1*, *STIM1*, *ORAI1* and *ATP2A3* mRNAs increased 1.7- to 3.2-fold, whereas SOCE fell from 1.542 ± 0.282 to 0.245 ± 0.083 Δ(F340/F380) (p = 0.0115; n = 3 per group) and endoplasmic reticulum (ER) Ca²⁺ release did not differ (p = 0.299). Composition-adjusted transcript summaries showed ORAI3 rising about four-fold and SARAF by about half, while the dominant ORAI2 and ATP2A2 fell modestly. Annotation-free analysis recovered canonical cryptic targets (including *STMN2* and *UNC13A*), none in FUS or TAF15 knockdown, but found no high-confidence cryptic junction in the core SOCE panel of the primary model; APA and NMD analyses confirmed no SOCE-machinery event. In ALS tissue *TRPC1* and *SARAF* were increased and *CBARP* decreased in six of seven brain regions, whereas *TRPC1* moved the other way in all comparison cohorts. Cryptic *STMN2* inclusion tracked TDP-43 pathology regionally but did not correlate with *TRPC1* within ALS regions.

**Conclusions.** TDP-43 knockdown is associated with reduced SOCE and altered expression and splicing of calcium-regulatory genes in SH-SY5Y cells, but no single causal RNA switch was identified, and the ALS-tissue associations are not shown to arise from TDP-43 loss. *CBARP* and *SARAF* are priorities for follow-up.

**Keywords:** TDP-43; store-operated Ca²⁺ entry; calcium homeostasis; alternative splicing; SARAF; CBARP; ALS; neurodegeneration

## 1. Introduction

TAR DNA-binding protein 43 (TDP-43, encoded by *TARDBP*) is a ubiquitously expressed RNA-binding protein that governs splicing, transcript stability and transport. Nuclear clearance with cytoplasmic accumulation of TDP-43 is found in approximately 97% of amyotrophic lateral sclerosis (ALS) cases and about half of frontotemporal lobar degeneration cases (Neumann et al., 2006; Arai et al., 2006; Ling et al., 2013), and in limbic-predominant age-related TDP-43 encephalopathy (Nelson et al., 2019), making TDP-43 dysfunction one of the few molecular events shared across otherwise distinct neurodegenerative syndromes.

The principal consequence of nuclear TDP-43 loss is failure to repress cryptic splicing. The best-characterised targets, *STMN2* and *UNC13A* (Ling et al., 2015; Klim et al., 2019; Melamed et al., 2019; Brown et al., 2022; Ma et al., 2022), are now the basis of antisense oligonucleotide strategies (Baughn et al., 2023), establishing that TDP-43-dependent mis-splicing is both measurable and, in principle, correctable. This has motivated systematic searches for further TDP-43-dependent RNA processing events with functional consequences.

Disturbed intracellular Ca²⁺ homeostasis is a longstanding feature of ALS. Vulnerable motor neurons are notable for low cytosolic Ca²⁺-buffering capacity (Grosskreutz et al., 2010), and Ca²⁺-dependent excitotoxicity, mitochondrial dysfunction and endoplasmic reticulum (ER) stress have all been implicated. Store-operated Ca²⁺ entry (SOCE) — the STIM1/ORAI-mediated influx triggered by ER store depletion (Putney, 2010; Prakriya and Lewis, 2015), with contributions from canonical transient receptor potential (TRPC) channels (Ambudkar et al., 2007) — is a principal route of Ca²⁺ entry in non-excitable and many excitable cells, and in a human neuroblastoma cell line it depends on TRPC1 (Selvaraj et al., 2012). STIMATE promotes the conformational switch that activates STIM1 (Jing et al., 2015); ORAI2 and ORAI3 can form heteromers with ORAI1 and restrain SOCE (Vaeth et al., 2017; Yoast et al., 2020), and SARAF facilitates slow Ca²⁺-dependent inactivation of SOCE (Palty et al., 2012) and, in SH-SY5Y cells, restrains store-independent arachidonate-regulated Ca²⁺ entry (Albarran et al., 2016). Disease-associated TDP-43 has been reported to disrupt VAPB–PTPIP51 tethering and thereby ER–mitochondrial Ca²⁺ transfer (Stoica et al., 2014), but whether TDP-43 loss affects the SOCE pathway itself, and whether it does so through RNA processing, has not been examined systematically.

Here we address that question at three levels. We first ask whether TDP-43 depletion produces reproducible splicing and isoform-usage changes in genes of the Ca²⁺ homeostasis machinery, applying robustness criteria that are stricter than the field default. We then ask what functional phenotype accompanies TDP-43 knockdown by measuring SOCE directly. Finally we ask whether the transcript-level changes we identify are recapitulated in ALS brain and how their direction compares with that in selected neurological comparison cohorts.

A recurring theme in what follows is that the level of analysis determines the answer. Individual splicing events in these data are frequently underpowered; the transcript-level changes are best treated as candidate explanations for, rather than demonstrated causes of, the functional phenotype.

## 2. Materials and Methods

### 2.1 Public RNA-seq datasets

Six comparisons with paired controls and publicly available raw reads were selected for the expression and event-level splicing analyses: SH-SY5Y (GSE296712), induced pluripotent stem cell (iPSC) colonies (GSE230647; Hruska-Plochan et al., 2024), iPSC-derived motor neurons (GSE77702; Kapeli et al., 2016), mouse striatum (GSE27394/GSE27218; Polymenidou et al., 2011), and mouse C2C12 and NSC34 (GSE171714/SRP314028; Šušnjar et al., 2022) — 38 libraries in total, counting the pooled technical fragments of GSE27394 as one library per biological replicate. Two further comparisons, K562 total RNA and poly(A)+ mRNA (ENCODE ENCSR372DZW/ENCSR455TNF and ENCSR129RWD/ENCSR134JRE; Van Nostrand et al., 2020), were added for the annotation-free junction analysis only (Section 2.5) and were not processed with rMATS or DESeq2. Inclusion required a condition reducing TARDBP expression or TDP-43 protein, at least two biological replicates and open raw reads, and for the human datasets the ability to evaluate *STMN2* and *UNC13A* as positive controls. TDP-43-repressed cryptic exons are largely not conserved between species (Ling et al., 2015), and the *STMN2* and *UNC13A* events are absent from the mouse genes (Melamed et al., 2019; Ma et al., 2022), so the mouse comparisons have no equivalent literature positive control. GSE296712, a subseries of the dataset reported by Bryce-Smith et al. (2025), was designated the primary human model because the laboratory experiments use the same cell line.

Reads were downloaded with SRA-Toolkit v3.2.1 (Leinonen et al., 2011), quality-assessed with FastQC v0.12.1 (Andrews, 2010), trimmed with Trim Galore v0.6.11 (Cutadapt v5.2, Martin, 2011; Illumina TruSeq adapter auto-detected, Phred Q20, minimum length 20 bases) and aligned with HISAT2 v2.2.2 (Kim et al., 2019) in `--dta` mode with a known-splice-site file, to the GRCh38 primary assembly or, for mouse data, to the mm10 index distributed with HISAT2; GENCODE v47 (basic annotation) and vM25 (Frankish et al., 2019) were used for read counting and for junction annotation; the transcript-level quantification of Section 2.4 uses the comprehensive GENCODE v47 annotation. Alignments were sorted and indexed with SAMtools v1.23 (Danecek et al., 2021) and summarised with MultiQC v1.25.1 (Ewels et al., 2016); the junction and coverage analyses of Sections 2.5 and 2.7 were run with SAMtools v1.21. Library strandedness was determined per dataset and rMATS was re-run accordingly (fr-secondstrand for SH-SY5Y; fr-firststrand for iPSC colony, iPSC-MN and mouse single-end).

For GSE27394, two technical sequencing fragments (SRR107072, SRR107073) showed lower alignment rates; these are fragments of the same biological replicate and counts were pooled within biological replicates, yielding four control and four knockdown replicates.

### 2.2 Differential expression

Gene counts were generated with featureCounts v2.1.1 (Liao et al., 2014) at meta-feature level against the same annotation (`-s 0 -g gene_name`, with `-p --countReadPairs` for paired-end libraries; multi-mapping and multi-overlapping reads not counted) and modelled with DESeq2 v1.42.1 (Love et al., 2014) within each dataset, with Benjamini–Hochberg correction (Benjamini and Hochberg, 1995). Genes with p_adj < 0.05 and |log2FC| ≥ 1 were called differentially expressed. Sensitivity to unmodelled variation was assessed with svaseq (sva v3.50.0; Leek, 2014); the principal findings were unchanged.

### 2.3 Alternative splicing and robustness assessment

Splicing was quantified with rMATS Turbo v4.3.0 (Shen et al., 2014) with `--novelSS` for skipped exons, alternative 5′ and 3′ splice sites, mutually exclusive exons and retained introns, using both junction-count (JC) and junction+exon-body (JCEC) models. Events were called at a false discovery rate (FDR) < 0.05 with an absolute change in percent spliced in (|ΔPSI|) ≥ 0.10.

Because rMATS applies FDR correction across all events irrespective of read support, we additionally applied a coverage pre-filter **before** correction, retaining events with ≥ 10 informative reads per sample on average and ≥ 5 informative reads in every sample, and recomputing Benjamini–Hochberg q-values within the filtered set. For the *TRPC1* event of Section 3.1, and for the events that met the thresholds and the coverage criteria in twelve SOCE-machinery genes (*STIM1*, *STIM2*, *STIMATE*, *SARAF*, *CRACR2A*, *CRACR2B*, *ORAI1–3*, *TRPC1*, *ATP2A2* and *ATP2A3*), we further computed (i) per-replicate PSI from raw inclusion and skipping junction counts (IJC/SJC), (ii) replicate-level bootstrap 95% confidence intervals for ΔPSI (20,000 resamples for the *TRPC1* analysis of Section 3.1; 10,000 for the SOCE candidate events of Section 3.2), and (iii) leave-one-out ΔPSI, removing each replicate in turn. Group assignment (b1 = knockdown, b2 = control) was verified from the labelled group files of GSE27394. Exon coordinates taken from rMATS output are reported here as 1-based, inclusive intervals. For the 24-nucleotide *STIM2* exon whose inclusion produces STIM2.1 (Section 3.3), the rMATS ΔPSI of the six datasets were combined by fixed-effect inverse-variance meta-analysis; the variance of each estimate was the between-replicate variance of PSI divided by the number of replicates, summed over the two groups; heterogeneity was assessed with Cochran's Q and I² (Supplementary Table S15).

Detection power was estimated by simulation over the observed coverage distribution for a 3 + 3 design.

### 2.4 Isoform usage

Transcript-level quantification used Salmon v1.11.4 (Patro et al., 2017) against an index built from the comprehensive GENCODE v47 transcript set (384,354 transcripts). Transcript estimates were aggregated to genes with the transcript-to-gene map of the same comprehensive annotation. The map of the basic annotation covers only 157,588 of those transcripts and would discard 26% of the transcripts per million at the aggregation step, so it was not used for any gene-level summary. Differential isoform usage was assessed with IsoformSwitchAnalyzeR v2.6.0 (Vitting-Seerup and Sandelin, 2019), DRIMSeq v1.34.0 (Nowicka and Robinson, 2016) and stageR v1.28.0 (Van den Berge et al., 2017), including prediction of open reading frames, premature termination codons, sensitivity to nonsense-mediated decay (NMD) and protein-domain consequences. As an annotation-independent comparison, LeafCutter v0.2.9 (Li et al., 2018) was run on the same alignments, and SUPPA2 v2.3 (Trincado et al., 2018) was run on the three control and three 75 ng/mL libraries as a check of method concordance.

### 2.5 Annotation-free junction-level analysis

Splice junctions were extracted directly from the alignments in two independent ways, so that no conclusion rests on a single extraction rule. The first set was produced with regtools v1.0.0 (Cotto et al., 2023; `junctions extract -s XS -a 8 -m 50 -M 1000000`). The second was produced with an in-house extractor that parses N operations from the CIGAR string and additionally requires a mapping quality (MAPQ) ≥ 30, primary non-duplicate alignment, an 8 bp matched anchor on both sides, and an intron length between 50 and 500,000 bp. Junctions are reported as 1-based, intron-inclusive intervals. On the primary SH-SY5Y comparison the two raw sets share 217,550 junctions, counted before strand-independent merging and before the read-support filters below; the sets differ in size because regtools applies no mapping-quality threshold and admits longer introns. All eleven comparisons were run on the regtools set; the SH-SY5Y comparison was run on both, and every conclusion drawn from it holds in each.

Junctions were classified against GENCODE v47 (human) or vM25 (mouse) as annotated, novel combination (both splice sites annotated, the junction not), novel site (one site unannotated), or fully novel. Local splicing variations (LSVs) were defined by shared donor and by shared acceptor, and the share of an LSV's reads carried by each junction defined its percent-spliced-in. Differences between groups were tested with a beta-binomial likelihood-ratio test using a dataset-wide overdispersion parameter chosen by profile likelihood, with Benjamini–Hochberg correction across all tested junctions, and replicate-level bootstrap 95% confidence intervals (2,000 resamples). An LSV required at least 20 reads per sample on average and non-zero coverage in every sample; a junction required at least 5 reads in at least one sample.

The calling threshold was calibrated against a control-versus-control null: in datasets with four control replicates, the controls were split 2 + 2, the analysis repeated unchanged, and every call counted as a false positive. Under a permissive definition (unannotated junction, ΔPSI ≥ 0.05, q < 0.05, control PSI ≤ 0.05, bootstrap lower bound above zero) the null produced almost as many calls as the real comparison in iPSC colonies (0.98 per real call) and twice as many in K562 total RNA (2.01; Supplementary Table S14), so raw call counts are not interpretable. We therefore adopted an effect-based definition — a **high-confidence unannotated splicing change** requires a novel splice site, ΔPSI ≥ 0.20, q < 0.05, a bootstrap lower bound above 0.05, at least 20 reads in the knockdown group, and non-zero counts in every knockdown replicate. A requirement of zero counts in controls was deliberately not imposed, because it removes genuinely cryptic exons with basal leakiness, *STMN2* among them. Positive-control recovery, not call count, is used throughout as the measure of whether an analysis worked; because the sixteen literature controls are human cryptic events, this check is available only for the human comparisons.

Eleven comparisons were analysed: SH-SY5Y at 75 and at 25 ng/mL doxycycline, iPSC colonies, iPSC-derived motor neurons with TDP-43, FUS or TAF15 knockdown, K562 total RNA and poly(A)+ mRNA, C2C12, NSC34, and mouse striatum. For mouse striatum, technical sequencing fragments were pooled within biological replicates, giving four per group.

### 2.6 Aberrant splicing outlier detection

FRASER v1.14.1 (Mertes et al., 2021) was run on the nine-sample SH-SY5Y doxycycline series (0, 25 and 75 ng/mL), counting split and non-split reads from the alignments and modelling ψ5, ψ3 and splicing efficiency θ. FRASER treats aberrant splicing as a rare, sample-specific deviation from a cohort norm; that assumption is not met when depleted samples form two-thirds of the cohort, and results are therefore reported as per-sample outlier burden with the limitation stated explicitly.

### 2.7 Alternative polyadenylation

Alternative and cryptic polyadenylation (APA) were assessed from coverage rather than from junctions. Windows were built for 696 genes: those of the 732-gene expanded Ca²⁺ panel for which annotated introns and terminal exons could be retrieved, together with the cryptic positive controls. Introns were taken as the gaps between the merged exons of all basic-annotation transcripts of a gene and numbered in the direction of transcription; these unit numbers need not match canonical intron numbers (the *STMN2* unit 'intron 2', for example, corresponds to the canonical intron 1, which contains cryptic exon 2a), and the genomic windows of every reported unit are listed in Supplementary Table S11. Two 500 bp windows were placed at the 5′ and 3′ ends of every such intron of at least 1,200 bp (50 bp buffer from each splice site), and the terminal exon of each gene of at least 400 bp was split into proximal and distal halves, both in the direction of transcription. An exon missing from the basic annotation that lies inside a window adds to that window's coverage, so a unit containing an alternatively spliced exon also reports that exon's inclusion. Depth was summed with `samtools bedcov` (MAPQ ≥ 30) while excluding CIGAR reference skips and deletions with `-j`. Excluding reference skips matters for an intronic index: counted as depth, the skips of spliced reads fill intronic windows, and in SH-SY5Y they would raise the *STMN2* intron 2 index difference from +0.145 to +0.489. Equivalent windows were built for 687 orthologous mouse genes. The analysis was run in the four knockdown comparisons whose alignments were available — SH-SY5Y, iPSC-derived motor neurons, C2C12 and NSC34 — and results are reported only for units whose input windows passed the depth filter.

The **intronic polyadenylation index** is the 5′ window's share of the two intronic windows. Premature termination within an intron raises the 5′ window relative to the 3′ window; intron retention raises both equally and leaves the index unchanged, so the index responds to a coverage gradient rather than to retention. The **distal 3′ untranslated region (3′UTR) usage index** is the distal half's share of the terminal exon, following the logic of DaPars (Xia et al., 2014). Both indices are reported as point estimates with bootstrap intervals obtained by complete enumeration of the replicate combinations rather than by sampling: 3⁶ = 729 combinations for the three-versus-three comparisons and 2⁴ = 16 for the two-versus-two iPSC-derived motor neurons, where the interval is correspondingly coarse. No p or q values are produced. A unit was retained only if, in every sample, the mean per-base depths of its two windows summed to at least 3.

### 2.8 Nonsense-mediated decay inhibition

Gene-level counts were obtained from GSE307054 (Sinha et al., 2025; i3Neurons; TDP-43 knockdown crossed with knockdown of XRN1, UPF1 and SMG6 in four combinations, two replicates each) and normalised with the authors' size factors. Sequencing batch is confounded with TDP-43 status in that design, so the interaction was formed as a difference of within-batch differences:

*interaction(c) = [log2(TDP-43 knockdown + NMD inhibition c) − log2(TDP-43 knockdown)] − [log2(control + NMD inhibition c) − log2(control)]*

for each of the four NMD-inhibition conditions c, giving four condition-level interaction values per gene; the two replicates within each condition are averaged. A positive interaction indicates a transcript that is produced upon TDP-43 loss and degraded by NMD. Genes with at least 10 normalised counts were tested (n = 19,145) with a four-condition one-sample t-test and an exact sign test. Because the TDP-43 factor is confounded with sequencing batch, these results are treated as exploratory.

### 2.9 Ca²⁺ gene panels and enrichment testing

Four cumulative gene panels were assembled from curated SOCE/TRP components (n = 51), channels and transport systems (n = 117), curated Ca²⁺-handling genes (n = 258) and an expanded Ca²⁺-associated set from KEGG and Gene Ontology (n = 732; Supplementary Table S2). Symbols were harmonised to HGNC/MGI and orthology verified via Ensembl Compara. Where results are summarised for the SOCE machinery, this refers to a fixed set of nineteen genes — *STIM1*, *STIM2*, *ORAI1–3*, *TRPC1*, *SARAF*, *STIMATE*, *CBARP*, *CRACR2A*, *CRACR2B*, *SELENOK*, *ATP2A1–3*, *MCU*, *MCUB*, *MICU1* and *MICU2*; the core SOCE/TRP panel is Tier 1.

Enrichment of significant splicing events within panels was tested both by hypergeometric test against all testable genes and against a **covariate-matched empirical null**: for each panel gene, the 25 nearest background genes in the space of log event count (a proxy for gene length and testability) and log total read support (a proxy for expression) were identified, and 5,000 permutations were drawn from these matched pools. The matched permutation p-value is reported as the primary result, because detection power for splicing events scales with gene length and expression and Ca²⁺ channel genes lie at the extreme of both distributions. Neither p value is corrected across the 24 combinations of dataset and panel.

### 2.10 Transcript-family abundance

Transcripts per million (TPM) from Salmon were summed to gene level for the primary SH-SY5Y comparison (0 vs 75 ng/mL doxycycline, n = 3 + 3), so that members of a family are compared on a length-corrected scale. TPM sums to the same total in every library, so a large gain by a few abundant transcripts lowers the TPM of every other gene: in the knockdown libraries *CHGA* rose from 1,857 to 10,675 TPM and ribosomal-protein and mitochondrially encoded transcripts also took a larger share, so that the median gene expressed in both groups (mean TPM > 5) had 24% lower TPM than in controls. Before conditions were compared, each library's TPM was therefore divided by a median-of-ratios factor computed over the 13,907 genes with TPM > 1 in all six libraries, the normalisation principle of DESeq2; after this adjustment gene-level changes agreed with the DESeq2 fold changes (median difference −0.02 log2 units). For each functional family (STIM, ORAI, sarco/endoplasmic reticulum Ca²⁺-ATPase [SERCA], TRPC, SOCE regulators, mitochondrial Ca²⁺ uptake and plasma-membrane Ca²⁺-ATPase [PMCA]) we computed each member's share of the family total in control cells and the change in its adjusted TPM. Fold changes and adjusted p-values for the same genes come from DESeq2 applied to the gene-level Salmon counts of the same six libraries. These summaries remain transcript-level measures and are not estimates of absolute molecule number, protein abundance or complex stoichiometry.

### 2.11 Patient tissue and cross-disease comparison

ALS post-mortem RNA-seq was taken from the New York Genome Center (NYGC) ALS Consortium/Target ALS collection (GSE153960; Prudencio et al., 2020; 1,641 samples with metadata after filtering). After a low-expression filter (counts per million [CPM] > 0.40), comparisons used log2 CPM with per-sample median centring, Mann–Whitney U tests, Benjamini–Hochberg correction within each region and Cliff's δ as effect size (Cliff, 1993). Ten regions were analysed at the expression level: seven brain regions (cerebellum; frontal, occipital and temporal cortex; lateral and medial motor cortex; hippocampus) and three spinal cord levels (cervical, thoracic, lumbar). The junction-level analysis of Section 3.10 covers eleven regions, adding motor cortex of unspecified subdivision. Within the same cohort, the "Other Neurological Disorders" group was compared against the **same** non-neurological controls as ALS, so that control-side confounders (RNA quality, post-mortem interval, collection site, library batch) are shared between the two comparisons; case numbers allowed this comparison in three regions.

Independent cohorts were analysed identically: Alzheimer's disease (GSE125583, fusiform gyrus, 219 cases/70 controls; Srinivasan et al., 2020), Parkinson's disease (GSE68719, prefrontal cortex, Brodmann area 9 [BA9], 29/44; Dumitriu et al., 2016) and multiple sclerosis (MS; GSE138614, white matter, 10 MS/5 control donors, 98 samples; Elkjaer et al., 2019; and GSE123496, five brain regions, 5/5 donors; Voskuhl et al., 2019). GSE138614 is white matter, whereas the ALS, Alzheimer's and Parkinson's comparisons are predominantly grey matter; the tissue compartment therefore differs and is not controlled for. For multiple sclerosis, TRPC1 contributes to SOCE in oligodendrocyte precursor cells (Paez et al., 2011), so a fall in TRPC1 could reflect loss of oligodendrocyte-lineage cells; two additional controls were therefore applied: analysis restricted to normal-appearing white matter, in which myelin markers are preserved, and regression of TRPC1 on MBP, PLP1 and GFAP with testing of the residuals. Several GSE138614 samples come from the same donor, so the lesion-type and normal-appearing white matter comparisons were run both at sample level and after averaging within donors; both are reported, and the donor-level values, which count each donor once, are treated as primary.

Cell-composition markers (GFAP, AIF1, SNAP25, RBFOX3, MBP, PLP1) were carried through all comparisons.

Junction-level quantification of the truncated *STMN2* transcript was obtained without access to the original alignments: the same SRA study (SRP270799) is present in recount3 (Wilks et al., 2021) as a junction count matrix of 12,653,421 junctions across 2,256 libraries, from which the 215,321 junctions falling in 171 target genes were extracted. Cryptic PSI was computed at the shared exon-1 donor as cryptic junction reads divided by all reads leaving that donor, requiring at least 20 reads at the donor. The cryptic junction coordinates were not taken from the literature but from our own de novo discovery in SH-SY5Y (Section 3.3); their agreement with published positions is therefore an independent check rather than an assumption. Libraries were linked to the cohort metadata through the CGND identifier.

### 2.12 Cell culture and TDP-43 knockdown

SH-SY5Y cells were maintained in DMEM/F12 with L-glutamine and 15 mM HEPES (Gibco, cat. no. 11330032) supplemented with 10% fetal bovine serum (FBS) and 1% MEM non-essential amino acids, at 37 °C in 5% CO₂, and were tested for mycoplasma by PCR and by a luminescence assay. TARDBP was silenced with a pLKO.1-TRC lentiviral short hairpin RNA (shRNA; Dharmacon TRC Lentiviral shRNA, cat. no. RHS3979); a non-targeting shRNA in the same vector served as control. Lentiviral particles were produced in HEK293T cells co-transfected with 1 µg transfer plasmid, 750 ng psPAX2 (Addgene 12260) and 250 ng pMD2.G (Addgene 12259) using LipoFectMax (ABP Biosciences; 3 µL per µg DNA); supernatants were collected 48 and 72 h after transfection, cleared by centrifugation, passed through a 0.45 µm filter and stored at −80 °C. SH-SY5Y cells were seeded at 375,000 per well in 6-well plates and transduced overnight with 20% (v/v) viral supernatant in medium containing 8 µg/mL polybrene (Sigma-Aldrich, H9268); the medium was replaced 24 h after transduction. Selection with 2 µg/mL puromycin (Cayman Chemical, 13884; concentration chosen from a kill curve) began no earlier than 24 h after transduction and was complete on day 5, when puromycin-treated non-transduced cells had all died. Three groups were compared: non-transduced cells, a non-targeting (scrambled) shRNA control and shTDP-43. TARDBP knockdown was measured against both control groups; the target-gene reverse-transcription quantitative PCR (RT-qPCR), Fura-2 and WST-1 experiments compared shTDP-43 cells with the non-targeting shRNA control. Lentiviral work followed institutional biosafety rules.

### 2.13 RT-qPCR

Total RNA was isolated with the Total RNA Purification Kit (Norgen Biotek, cat. no. 17200) with on-column DNase I treatment, and samples with A260/A280 > 1.8 were used. cDNA was synthesised from 500 ng total RNA with oligo(dT) priming and the OneScript Plus cDNA Synthesis Kit (abm, cat. no. G236). Reactions of 10 µL, containing 1 µL cDNA and 0.5 µM of each primer, were run with iTaq Universal SYBR Green Supermix (Bio-Rad, cat. no. 1725121) on a LightCycler 480 II (Roche). Targets and GAPDH, the single reference gene, were run on the same plate with no-template and no-reverse-transcriptase controls, and specificity was checked from a single melt-curve peak. Relative expression was calculated by the 2^−ΔΔCt^ method (Livak and Schmittgen, 2001), with non-transduced cells as the calibrator for TARDBP and the non-targeting shRNA control for the target genes (n = 4 biological replicates per group; TARDBP knockdown and the target genes were measured on separate RNA sets, 14 April – 21 May and 3 June – 10 July 2026). Primer sequences, product sizes, annealing temperatures and the thermal cycling profile are given in Supplementary Table S1.

### 2.14 Cytosolic Ca²⁺ measurement

Non-targeting shRNA control and shTDP-43 cells were measured 72 h after transduction. Cells were seeded in 24-well plates at 40,000 cells per well; the three samples of each group came from three independent cultures. Cells were loaded with 5 µM Fura-2/AM (Invitrogen, F1221) and 0.02% Pluronic F-127 (Invitrogen, P3000MP) in HEPES-buffered saline (HBS) containing 1% bovine serum albumin (BSA) for 60 min at 25 °C in the dark and washed three times for 15 min in 1% BSA/HBS. HBS contained 135 mM NaCl, 5.9 mM KCl, 1.2 mM MgCl₂, 1.5 mM CaCl₂, 11.6 mM HEPES, 5 mM NaHCO₃ and 11.5 mM D-glucose (pH 7.3). Cytosolic free Ca²⁺ was followed spectrofluorometrically as the F340/F380 ratio (excitation 340 and 380 nm, emission 510 nm; Grynkiewicz et al., 1985) on a cuvette-based QM8/2005 spectrofluorometer (Photon Technology International) with peristaltic perfusion, following the protocol of Selli et al. (2009). Loaded cells were transferred to Ca²⁺-free HBS containing 1 mM EGTA; SERCA was then inhibited with 10 µM cyclopiazonic acid (Sigma-Aldrich, cat. no. C1530) and the transient rise in F340/F380 was recorded as ER Ca²⁺ release. CaCl₂ was then added to a final concentration of 1.5 mM and the subsequent rise was recorded as store-operated Ca²⁺ entry. Both were quantified as Δ(F340/F380) relative to the respective preceding baseline (n = 3 per group).

### 2.15 Viability

Viability was measured at 24 h and 48 h after seeding with the WST-1 assay (Premix WST-1, Takara Bio, cat. no. MK400). Cells were seeded at 10,000 per well in 100 µL in 96-well plates; at each time point 10 µL reagent was added for 3 h at 37 °C, and absorbance was read at 450 nm against a 620 nm reference on a Varioskan Flash reader (Thermo Scientific). The absorbance of cell-free wells containing medium and reagent was subtracted, and values were normalised to the mean of the non-targeting shRNA control at the same time point. The experiment was performed three times; the values and statistics reported are those of one experiment (n = 4 wells per group).

### 2.16 Statistics

Bioinformatic thresholds are given above. Laboratory data are mean ± SEM. TARDBP knockdown was assessed by one-way ANOVA on log2-transformed relative expression with Tukey's multiple comparison test; other laboratory comparisons used two-tailed Student's t-tests (GraphPad Prism 10). No correction for multiple comparisons was applied across the four RT-qPCR targets. Significance was set at p < 0.05.

### 2.17 Data and code availability

All RNA-seq datasets analysed are public and none were generated for this study; accessions are listed in Supplementary Table S17. The analysis package — scripts in execution order, including the coverage, NMD, transcript-family and donor-level scripts used for the results reported here, the four Ca²⁺ gene panels, per-event count tables, junction and LSV tables, and the code that draws every figure — is deposited at **https://github.com/elmasnuryilmaz/tdp43-soce-manuscript**, which is private during review and will be opened on acceptance; access for editors and reviewers is available on request. The repository will be archived with a Zenodo DOI on acceptance, and the complete 23 GB rMATS output will be deposited in the same archive. The package downloads the two externally hosted resources it needs (the recount3 junction matrix for SRP270799 and the GSE307054 count table) at run time. Every rMATS event that meets the thresholds of Section 2.3, with its raw junction counts, is provided as Supplementary Table S3. RT-qPCR, Fura-2 and WST-1 source data, including per-replicate values and the statistics behind every laboratory figure panel, are provided in Supplementary Table S1.

## 3. Results

### 3.1 TDP-43 depletion alters splicing across Ca²⁺ homeostasis genes, but most individual events are underpowered

Across the six comparisons, TDP-43 depletion produced widespread splicing changes; in the primary SH-SY5Y model 7,854 events met FDR < 0.05 and |ΔPSI| ≥ 0.10 before coverage filtering. Positive-control behaviour of *STMN2* and *UNC13A* was consistent with published direction of change. A local-event method run on the same libraries gave a more conservative picture: SUPPA2 tested 31,474 events in 5,345 genes and returned no event that met the shared FDR and |ΔPSI| ≥ 0.10 thresholds, and among the 9,369 skipped-exon events whose coordinates matched between the two tools ΔPSI agreed only weakly (Spearman ρ = 0.319; 62.4% directional concordance). SUPPA2 builds its null from the between-replicate ΔPSI distribution, which makes it conservative in a three-versus-three design, so the disagreement bounds how much weight any single event call can carry rather than showing that one tool is wrong.

However, a substantial proportion of nominally significant events rested on very low read support. Applying a coverage pre-filter before FDR correction removed 18–78% of tested events depending on dataset (Table 1) and reduced the number of significant calls by 33–76%, with the largest losses in the two single-end datasets (iPSC-derived motor neurons 70%; mouse striatum 76%).

Simulation over the observed coverage distribution showed why (Figure 1). For a 3 + 3 design at 10 informative reads per sample — a coverage typical of the lower quartile of events — power to detect a true ΔPSI of 0.10 is only 0.08, and even at 100 reads per sample it reaches only 0.28. Attaining 80% power at realistic depth requires |ΔPSI| ≳ 0.30. Nominally significant events at low coverage are therefore expected to be substantially inflated in effect size, and individual events at the |ΔPSI| ≥ 0.10 threshold should not be interpreted without explicit read-level support.

We applied this reasoning to a *TRPC1* skipped-exon event (chr3:142,792,824–142,792,967; ΔPSI = +0.108, FDR = 0.0459) that passes the conventional thresholds and would be prioritised on them alone. Replicate-level inspection showed that the entire signal derived from seven skipping reads across six libraries, with the skipping form absent in four of six samples; the bootstrap 95% confidence interval spanned zero (−0.093 to +0.522); removing one control replicate reversed the sign of ΔPSI (−0.074); and the event did not survive coverage pre-filtering. LeafCutter did not call *TRPC1* in the same comparison. We therefore do not report this event as a finding, and present it instead as an illustration of why threshold-based prioritisation of splicing events requires read-level verification (Figure 2).

### 3.2 Robust splicing changes concentrate in SOCE regulators

Applying the robustness criteria to the twelve SOCE-machinery genes examined event by event (Methods 2.3) identified a small set of events supported by adequate coverage and bootstrap intervals excluding zero (Table 2, Figure 3).

In the primary SH-SY5Y model these were *STIMATE* (ΔPSI = +0.244; FDR = 9.0 × 10⁻⁷; 95% CI +0.095 to +0.368), *ORAI3* (−0.269; 1.0 × 10⁻²; −0.404 to −0.107) and *STIM2* (−0.120; < 1 × 10⁻¹⁶; −0.165 to −0.064). A *STIM1* event (chr11:4,088,702–4,088,738; 37 bp, frame-disrupting) reached ΔPSI = +0.145 (FDR = 4.0 × 10⁻⁴) with a confidence interval marginally including zero, and was positive in all three human datasets. Isoform-level testing supported *STIM1* (DRIMSeq gene-level q = 5.49 × 10⁻⁷) with two isoforms predicted to carry premature termination codons.

*CBARP*, encoding a suppressor of voltage-gated Ca²⁺ channel activity and Ca²⁺-evoked exocytosis (Béguin et al., 2014), was the most reproducible candidate: 32 events met both the threshold and the coverage criteria, in five of the six datasets — iPSC colonies, mouse striatum, SH-SY5Y, C2C12 and NSC34 — at the orthologous human chr19 and mouse chr10 loci, with |ΔPSI| from 0.11 to 0.74 (median 0.37) and 61 to 4,559 junction reads per event (median 246; 26 of the 32 above 100 reads; Supplementary Table S3). The locus is therefore affected in every model except the iPSC-derived motor neurons, but the direction is not consistent across models: all three qualifying SH-SY5Y events show reduced inclusion of the rMATS-defined form, whereas the events in mouse striatum and C2C12, all but one of those in NSC34 and most of those in iPSC colonies are positive; the junction-level analysis of Section 3.3 quantifies a different arm of the same locus and is positive in SH-SY5Y. What reproduces is the involvement of the locus, not one directional switch. *CBARP* also showed the largest expression change of any SOCE regulator in SH-SY5Y (log2FC = −1.254; p_adj = 3.1 × 10⁻²⁵).

LeafCutter, which does not use the reference annotation to define events, independently called *CBARP* (p_adj = 0.0133) and *STIM2* (p_adj = 0.0288) in SH-SY5Y. Three approaches — coverage-filtered rMATS, bootstrap interval estimation and LeafCutter — therefore converge on the same candidates.

### 3.3 Annotation-free analysis recovers cryptic exons genome-wide but finds no high-confidence event in the SOCE machinery of the primary model

rMATS and LeafCutter both constrain what can be found: the first by the event classes and annotation it works from, the second by clustering rules that discard sparse junctions. To remove those constraints we extracted splice junctions directly from the alignments and tested local splicing variations with a beta-binomial model, across eleven comparisons; the primary SH-SY5Y comparison was additionally repeated on a second, independently produced junction set (Methods 2.5).

The approach recovered the established TDP-43 cryptic targets without being given their coordinates (Figure 4A). In SH-SY5Y at 75 ng/mL doxycycline, 188,491 junctions passed the read-support filters with the mapping-quality-filtered extractor, 23,421 of them absent from GENCODE v47, and 117 high-confidence unannotated splicing changes in 87 genes were called (Supplementary Table S5). Twelve of the sixteen literature positive controls (Supplementary Table S6) were among them: *STMN2* (chr8:79,611,215–79,616,821; PSI 0.985 in knockdown versus 0.028 in control, ΔPSI +0.957, q = 2.0 × 10⁻³⁰⁴, 10,926 versus 176 reads), *UNC13A*, *ACTL6B*, *PFKP*, *HDGFL2*, *AGRN*, *ARHGAP32*, *ATG4B*, *ELAVL3*, *SETD5*, *RSF1* and *GPSM2*. The acceptor of the *STMN2* cryptic junction falls one base before the published start of the cryptic exon, and both flanking junctions of the *UNC13A* cryptic exon were recovered separately, placing that exon at chr19:17,642,414–17,642,541. On the regtools set the same comparison gave 165 high-confidence events in 113 genes and the same twelve positive controls plus *KALRN*; twelve were recovered at 25 ng/mL (Table 3; Supplementary Table S12).

Specificity was tested by running the same analysis on FUS and TAF15 knockdown in the same iPSC-derived motor neurons, at the same depth and with the same design (Figure 5A). Under the permissive definition the three comparisons produced call counts of the same order (141 for TDP-43, 124 for FUS and 126 for TAF15); TDP-43 knockdown recovered two of the sixteen positive-control genes (*STMN2* and *KALRN*) and FUS and TAF15 knockdown recovered **none**. This dataset is too shallow for the high-confidence threshold to recover a positive control in any of the three comparisons (18, 26 and 23 events, none of them a control gene), so the within-dataset contrast rests on the permissive definition. In the two deeper TDP-43 comparisons the same permissive analysis recovered thirteen of sixteen positive controls in SH-SY5Y and fifteen in iPSC colonies. The background rate of the procedure is not specific to TDP-43; the genes it identifies are.

The null test also sets the limits of interpretation. It could be run in the three datasets with four control replicates: for every real call the split-control null produced 0.64 calls in iPSC colonies, 2.17 in K562 total RNA and 0.83 in mouse striatum (high-confidence definition; Supplementary Table S14). These ratios are not calibrated false-discovery rates — they are the ratio of two call counts from analyses with different sample sizes and different numbers of tests — but they show that call counts are not interpretable on their own. In the remaining eight comparisons — the two SH-SY5Y doses, the three iPSC-derived motor neuron comparisons, K562 poly(A)+ mRNA, C2C12 and NSC34 — the controls could not be split, and positive-control recovery is the only available check. Across the eleven comparisons the burden of high-confidence unannotated splicing changes ranged from 12 to 477 events, and of the three datasets with a null, the one with the largest burden also had the largest null.

Applied to the Ca²⁺ panels, the analysis returned an almost complete negative (Supplementary Table S7). **No gene of the core SOCE/TRP panel (Tier 1, n = 51) carried a high-confidence unannotated splicing change in any comparison except iPSC colonies.** The negative is informative where the analysis demonstrably worked — in SH-SY5Y at both doses and in K562 poly(A)+ mRNA, where positive controls were recovered at the high-confidence threshold, and less strongly in the TDP-43 knockdown of iPSC-derived motor neurons, where two were recovered under the permissive definition only. In K562 total RNA no positive control was recovered, and the three mouse comparisons have no conserved control, so the absence of calls there carries little information. The single exception was iPSC colonies, where *CBARP* and *TRPM3* were called — in the dataset whose null test produced 0.64 calls for every real call, so a lone call there is not evidence, although the *CBARP* junction is corroborated elsewhere (below). *STIM1*, *STIM2*, *ORAI1–3*, *TRPC1*, *SARAF*, *STIMATE* and the SERCA and mitochondrial uptake genes carried none in any comparison. Under the permissive definition, SH-SY5Y additionally showed events in *PLCD4* and in six genes of the expanded Ca²⁺ panel (*ADCY1*, *SYT7*, *PDE3B*, *AKT3*, *PIK3CB*, *PPP2R5C*), and novel combinations of annotated splice sites appeared in *ATP2A3* (both doses, ΔPSI +0.16 and +0.20) and *MCUB* (iPSC colonies, +0.27) (Figure 4B). Within the coverage, models and calling criteria used here, we did not detect high-confidence unannotated splicing changes in the store-operated entry machinery. This finding does not exclude lower-abundance or context-specific events, but it provides no support for a cryptic-splicing switch as the explanation for the observed SOCE phenotype.

What the annotation-free analysis did confirm were the annotated events. *CBARP* showed the largest local change of any SOCE gene in SH-SY5Y (ΔPSI +0.74, 95% CI +0.51 to +0.83, q = 2.7 × 10⁻¹²; +0.76, q = 5.6 × 10⁻¹⁷ in the regtools set), and one arm of the *CBARP* event uses a splice site absent from the annotation. This is the third method, after rMATS with `--novelSS` and LeafCutter, to converge on this locus; isoform-level testing did not detect a *CBARP* isoform switch (IsoformSwitchAnalyzeR gene-level q = 0.12; DRIMSeq q = 0.34). The same unannotated *CBARP* splice site was recovered independently in iPSC colonies (ΔPSI +0.33, q = 4.3 × 10⁻¹⁶, 201 versus 137 reads), where it met the cryptic criteria outright. *ORAI2*, the most abundant ORAI-family gene in the transcript estimates, showed a change of +0.16 (q = 1.9 × 10⁻³) in the mapping-quality-filtered set and +0.15 (q = 0.087) in the unfiltered set; we report it as supported but not established.

One specific mechanism could be tested directly. Inclusion of a 24-nucleotide exon in the STIM–ORAI activating region (SOAR) converts STIM2 into STIM2.1 (STIM2β), an isoform that inhibits SOCE (Miederer et al., 2015; Rana et al., 2015), so increased inclusion upon TDP-43 loss would be a simple route to reduced entry. Measured as an rMATS event, inclusion of the exon was higher in knockdown in five of six datasets, but the fixed-effect meta-analysis estimate was negligible (pooled ΔPSI +0.0013, 95% CI −0.022 to +0.024; p = 0.914; Supplementary Table S15). The three mouse datasets carry 86% of the weight, but the human datasets alone give the same answer (+0.031, 95% CI −0.030 to +0.091; p = 0.32), and the estimates are not heterogeneous (Cochran's Q = 4.84, 5 df, p = 0.44; I² = 0%). At junction level, inclusion of the SOAR exon (chr4:27,007,983–27,008,006) was positive in direction in all five TDP-43 comparisons in which it was measurable, but the effect was small and inconsistent: ΔPSI +0.149 in SH-SY5Y at 75 ng/mL (q = 0.033, 20 versus 3 reads), +0.116 at 25 ng/mL (q = 0.14), and +0.020 in iPSC colonies (q = 0.82), the best-powered comparison with 441 versus 281 reads (Supplementary Table S13). The junction-level estimate is larger than the rMATS estimate for the same exon in the same libraries (+0.149 versus +0.059; Supplementary Table S15) because the two are not the same quantity: rMATS combines both flanking junctions and normalises by effective length, whereas the junction-level PSI is the share of the single downstream junction among all junctions leaving that donor, which is the more sensitive and the noisier of the two for a lowly used exon. On either measure the direction is reproducible and the magnitude is too small to support an isoform switch, and the negative conclusion of the meta-analysis stands.

FRASER, run on the nine-sample doxycycline series, returned no genome-wide significant outlier and no difference in per-sample outlier burden between depleted and control libraries (4.8 versus 2.3 events at p < 10⁻⁵; p = 0.35). This is the expected behaviour of an outlier method applied to a cohort in which the aberrant state is the majority, and we report it as a limit of that approach at this cohort size rather than as evidence against aberrant splicing.

### 3.4 Exploratory NMD interactions prioritise CBARP but do not confirm an NMD target

Cryptic exons frequently introduce premature termination codons, which can make their transcripts susceptible to nonsense-mediated decay (NMD) and under-represented in steady-state RNA. An independent experiment in which TDP-43 knockdown is crossed with knockdown of XRN1, UPF1 and SMG6 (GSE307054) provides an exploratory interaction test (Figure 4C): a transcript pool produced upon TDP-43 loss and degraded by NMD would be expected to rise preferentially when both perturbations are present (Methods 2.8).

Because the four interventions reuse the same control and knockdown libraries, the eight difference-of-differences values that the design yields are not eight independent observations. We therefore took the four NMD-inhibition conditions as the unit of inference, averaging the two replicates within each condition and retaining the shared controls only within each contrast. Several positive-control and calcium-gene interactions were positive, but no gene passed genome-wide FDR correction in this conservative analysis (smallest q = 0.053). The result is therefore supportive at the level of a ranked exploratory signal, not a confirmatory NMD test.

*CBARP* had the largest SOCE-panel interaction (+1.52 log₂; positive in all four conditions; Supplementary Table S10), but its conservative condition-level q-value was 0.145 (two-sided t-test; exact one-sided sign-test q = 0.310). It is therefore a prioritised candidate rather than an established NMD target. *RYR2* and *MCU* also had positive interactions, without genome-wide significance.

Two negative results follow. First, the prediction from isoform-level testing that *STIM1* produces a premature-termination-codon-bearing, NMD-sensitive isoform was **not supported**: in the four-condition analysis the interaction was −0.325 (p = 0.317, q = 0.508), in the opposite direction. Second, in the same four-condition analysis none of the four Ca²⁺ panels showed collective NMD sensitivity (one-sided Mann–Whitney p = 0.44, 0.91, 0.96 and 0.91 for Tiers 1–4), and the sixteen positive-control genes were not NMD-sensitive as a group either (p = 0.62; Supplementary Table S10b), consistent with the absence of cryptic events in the SOCE machinery.

### 3.5 A coverage-based APA screen yields candidate gradients but no confirmed APA event in the SOCE machinery

Cryptic polyadenylation accounts for a substantial part of TDP-43-dependent RNA processing (Bryce-Smith et al., 2025) — the truncated *STMN2* transcript is generated this way — and had not been assessed for the Ca²⁺ gene set. We measured it from coverage, as the 5′ share of the two ends of each intron (an intronic polyadenylation index) and as the distal share of each terminal exon (a 3′UTR usage index), across the 696 qualifying genes of the expanded Ca²⁺ panel and the cryptic positive controls (Methods 2.7; Figure 5B; Supplementary Table S11).

In the primary SH-SY5Y model the positive control shifted only moderately: using all six libraries, the *STMN2* intronic index gave a group difference of +0.145 (bootstrap 95% CI +0.085 to +0.205). *STMN2* falls outside the depth-filtered summary table because one control library has a summed two-window depth of 2.97, just below the pre-specified threshold of 3; one knockdown library has a measured coverage of zero in the 3′ window, which by the definition of the index contributes a value of 1.0 rather than missing data, and excluding that library gives +0.121 (+0.068 to +0.173) (Figure 5C).

Among the depth-qualified units of the SOCE machinery in the full-panel analysis, two intronic units showed moderate decreases (*STIM1* intron 17, Δ = −0.173, bootstrap 95% CI −0.233 to −0.123; *STIM2* intron 13, Δ = −0.156, −0.239 to −0.037). The *STIM2* unit is not an independent observation: its 5′ window contains the alternatively spliced *STIM2* exon of Table 2 (chr4:27,021,494–27,021,612), whose inclusion falls in knockdown, and exon reads dominate that window (control index 0.96), so the decrease reports that splicing change. The *STIM1* unit is the intron immediately upstream of the alternatively included *STIM1* exon of Table 2 (its 3′ window ends 50 bp before that exon) and may be influenced by it as well. Both are candidate coverage gradients rather than localised poly(A) sites. Sparse intronic windows produced unstable point estimates and were excluded by the pre-specified depth filter.

The analysis therefore identifies at most one moderate candidate gradient in the SOCE machinery, in *STIM1*, and does not establish a direct APA event in any of its genes. No index change of 0.05 or more was found for *ORAI1–3*, *TRPC1*, *SARAF*, *STIMATE*, *CBARP* or the SERCA genes, and the two units of those genes whose intervals exclude zero are negligible in size (*SARAF* terminal exon −0.036, *ORAI2* terminal exon +0.012). Outside the SOCE machinery three units likewise exclude zero — the *ITPR3* terminal exon (−0.152), the *MICU3* terminal exon (+0.097) and the *ITPR1* terminal exon (−0.055) — so the ER-release and mitochondrial-uptake arms show the same kind of candidate gradient as *STIM1* rather than being spared. Because the positive control responded only moderately in this model, absence of a signal should not be interpreted as evidence that these genes are definitively spared.

We therefore repeated the analysis in the iPSC-derived motor neurons, where 59 units in 29 genes passed the same depth filter. Here the assay's positive control behaved as intended: the index of *STMN2* intron 2, which contains cryptic exon 2a, rose from 0.570 in controls to 0.819 in knockdown (Δ = +0.249, interval +0.208 to +0.290). Against that working positive control, the SOCE-machinery genes again showed only small shifts — the largest were *ATP2A2* intron 3 (−0.164) and, below 0.10, *SARAF* intron 5 (+0.097), *TRPC1* intron 1 (+0.083) and the *STIM2* terminal exon (−0.074) — and the SH-SY5Y *STIM1* intron 17 unit did not reach the depth filter in this shallower dataset (*STIM2* intron 13 gave −0.024, interval spanning zero). This comparison has only two replicates per group, so the enumerated bootstrap has sixteen draws and its intervals are coarse; the point estimates, not the intervals, carry the information.

The two mouse lines gave 74 qualifying units in C2C12 and 131 in NSC34. Neither has a positive control for this assay: the *STMN2* cryptic polyadenylation site is absent from the mouse gene (Melamed et al., 2019), so the *Stmn2* units cannot serve as one. No SOCE-machinery unit exceeded |Δ| = 0.30 in either line, and the two largest (*Atp2a2* intron 6, +0.285, and *Trpc1* intron 7, +0.260, both in C2C12) have intervals that include zero and are not reproduced in the other line. One unit is positive in both motor-neuron models: the *SARAF* intron 5 index rises in the iPSC-derived motor neurons (+0.097) and the *Saraf* intron 5 index in NSC34 (+0.204), while in SH-SY5Y the same unit is uninformative (0.000, interval −0.264 to +0.241). We record it as a candidate rather than a finding: the NSC34 interval is wide, the myoblast line shows nothing there, the human and mouse units are matched by number rather than by sequence alignment, and a coverage gradient is not a poly(A) site.

Taken together, these analyses support altered calcium-related RNA profiles but do not establish whether any SOCE-machinery gene is a direct target of cryptic splicing, APA or NMD-coupled degradation.

### 3.6 Enrichment of splicing changes in Ca²⁺ genes is model-dependent

Testing whether Ca²⁺ genes are enriched for splicing changes required care, because power to detect an event scales with gene length, exon number and expression, and Ca²⁺ channel genes are extreme on all three.

In SH-SY5Y, the expanded Ca²⁺ panel showed 133 of 383 testable genes with at least one significant event (34.7%) against a background of 30.7%, giving an uncorrected hypergeometric p = 0.051. Against the covariate-matched empirical null, however, the expected rate was 32.4% and the permutation p-value was 0.139 (Supplementary Table S4). The apparent enrichment in this model is therefore largely attributable to the length and expression properties of Ca²⁺ genes.

The result differed by model. In iPSC-derived motor neurons — the most disease-relevant system examined — enrichment survived matching for both the channel/transport panel (p = 0.046) and the curated Ca²⁺-handling panel (p = 0.032), although neither p value would survive correction across the 24 dataset–panel combinations tested. In iPSC colonies, by contrast, nominally significant enrichment (hypergeometric p = 0.018–0.022) disappeared entirely against the matched null (p = 0.54–0.64).

We therefore regard enrichment of splicing changes in Ca²⁺ homeostasis genes as a tentative, motor-neuron-associated observation rather than a general property of TDP-43 depletion.

### 3.7 TDP-43 knockdown reduces SOCE despite increased transcript levels of SOCE components

Lentiviral shRNA reduced TARDBP mRNA by 94.4% relative to the non-targeting shRNA control and 94.8% relative to non-transduced cells (one-way ANOVA with Tukey's test, p < 0.0001, n = 4).

Relative to the non-targeting shRNA control, mRNA levels of four SOCE-associated targets increased: *TRPC1* ≈ 1.8-fold (p < 0.001), *STIM1* ≈ 1.9-fold (p < 0.001), *ORAI1* ≈ 1.7-fold (p < 0.01) and *ATP2A3*/SERCA3 ≈ 3.2-fold (p < 0.0001; n = 4, two-tailed t-tests). All four directions matched the RNA-seq results in the same cell line (*STIM1* log2FC = 0.929, p_adj = 4.3 × 10⁻⁴⁷; *TRPC1* 0.958, 1.3 × 10⁻¹²; *ORAI1* 0.433, 2.4 × 10⁻⁴; *ATP2A3* 1.306, 1.4 × 10⁻²⁹).

Functionally, however, the SOCE response fell markedly (Figure 6). ER Ca²⁺ release did not differ significantly between groups (control 0.268 ± 0.042 versus shTDP-43 0.180 ± 0.061 Δ(F340/F380); p = 0.299), whereas store-operated Ca²⁺ entry decreased from 1.542 ± 0.282 to 0.245 ± 0.083 (p = 0.0115; n = 3 per group). Viability moved in opposite directions at the two time points measured: 116.8 ± 1.6% of control at 24 h (p = 2.4 × 10⁻⁴) and 61.5 ± 0.8% at 48 h (p = 2.1 × 10⁻⁷; n = 4 wells from one of three independent experiments).

The direction of the transcript changes and the direction of the functional response are therefore opposite. This is the central observation of the study.

### 3.8 Transcript-family abundance provides candidate explanations for the discordance

The apparent contradiction can be explored by comparing transcript-family abundance rather than fold-change alone (Figure 7; Table 4). Because a few abundant transcripts take a larger share of the TPM total in the knockdown libraries, unadjusted TPM places every gene about a quarter lower in knockdown than its true change; the values below are adjusted for library composition (Methods 2.10).

In control cells *ORAI2* and *ATP2A2* are the dominant members of their families (77% and 98% of the family total), whereas *ORAI1* and *ATP2A3* are minor members (15% and 1.5%); *STIM1* and *TRPC1* dominate their families (64% and 98%), and *SARAF* accounts for 82% of the SOCE-regulator pool. In knockdown the dominant *ORAI2* and *ATP2A2* transcripts fell modestly (DESeq2 log2FC −0.17 and −0.25; p_adj = 2.2 × 10⁻³ and 8.7 × 10⁻⁹), while *ORAI1*, *ATP2A3*, *STIM1*, *TRPC1* and *SARAF* rose and *ORAI3* rose about four-fold (log2FC 2.06; p_adj = 2.3 × 10⁻⁶⁴). At family level the ORAI pool therefore grew (+27%) and changed in composition: *ORAI3* went from 8% to 30% of ORAI transcripts and *ORAI2* from 77% to 53%. The STIM (+48%) and SOCE-regulator (+30%) pools also rose, whereas the SERCA and mitochondrial-uptake pools fell (−12% and −25%), the latter mainly through *MCU*, *MICU2* and *MCUB*. These are descriptive summaries of transcript abundance and cannot be converted directly into protein-complex stoichiometry.

Two of these changes point towards reduced entry. ORAI2 and ORAI3 both form heteromers with ORAI1 and restrain SOCE — deleting either increases it (Vaeth et al., 2017; Yoast et al., 2020) — so the rise in *ORAI3* is a candidate contributor to the phenotype, whereas the fall in *ORAI2* would, on the same evidence, be expected to increase entry. *SARAF*, which facilitates slow Ca²⁺-dependent inactivation of SOCE (Palty et al., 2012), rose by about half (log2FC 0.55). The mitochondrial panel also changes, but MCU, MCUB and MICU proteins have distinct and context-dependent roles (Lambert et al., 2019), so their summed transcript change does not establish a change in mitochondrial Ca²⁺ uptake. The STIM1:ORAI1 RNA ratio is not a measurement of the protein ratio at ER–plasma-membrane junctions (Hoover and Lewis, 2011) and is not used here as a mechanistic readout.

These transcript patterns are compatible with a compositional contribution to reduced SOCE, but protein abundance, membrane localisation, channel assembly and causal rescue were not measured.

The transcriptome state therefore provides testable candidate explanations for the SOCE phenotype, rather than a demonstrated mechanism.

### 3.9 TRPC1, SARAF and CBARP differ across ALS and neurological comparison cohorts

We next asked whether these transcript-level changes are present in patient tissue (Figure 8; Table 5).

In the NYGC ALS cohort, *TRPC1* was increased in six of the seven brain regions tested (Cliff's δ = +0.447 to +0.699, q < 0.05), the exception being occipital cortex, in which no gene reached significance in any comparison. *SARAF* was increased and *CBARP* decreased in the same six regions, and both also changed significantly in cervical and lumbar cord (*SARAF* δ = +0.48 and +0.50; *CBARP* δ = −0.46 and −0.45). *ORAI1* — a lower-relative-abundance member of its transcript family in the cellular model — showed no significant change in any ALS region. *ATP2A3*, the other low-abundance member that rises in the cellular model, moved in the opposite direction in tissue: it was significantly decreased in frontal and lateral motor cortex, hippocampus, and cervical and lumbar cord (δ = −0.36 to −0.57).

Notably, the ALS increase in *TRPC1* was confined to brain: none of the three spinal cord levels reached significance, and thoracic cord trended negative (δ = −0.321). Cervical and lumbar cord are well powered in this cohort, so the absence there is informative; thoracic cord and occipital cortex are not, since no gene reached significance in either.

The "Other Neurological Disorders" group within the same cohort was compared against the same controls. *TRPC1* decreased in all three testable regions (δ = −0.407 to −0.717), the opposite direction to ALS. Shared controls remove confounding specific to the control samples, but do not resolve differences in disease composition, brain region, case processing or cell composition between the case groups. GFAP and SNAP25 moved more strongly in the comparison group than in ALS while *TRPC1* fell; these markers make a simple neurodegeneration explanation less compelling but do not provide comprehensive cell-composition adjustment.

Three independent cohorts extended this. *TRPC1* was decreased in Alzheimer's disease (δ = −0.447; q = 1.6 × 10⁻⁷, with Braak-stage correlation ρ = −0.183, p = 1.8 × 10⁻³) and in Parkinson's disease (δ = −0.677; q = 1.8 × 10⁻⁴). In multiple sclerosis, *TRPC1* was decreased at donor level across all sampled lesion types (δ = −0.840; q = 0.038) and in lesions at sample level (δ = −0.594; q = 1.3 × 10⁻⁴); in the second multiple sclerosis cohort the direction was the same but the study was underpowered (five regions pooled, δ = −0.226; q = 0.49). Because TRPC1 contributes to SOCE in oligodendrocyte precursor cells (Paez et al., 2011), we asked whether the decrease simply reflected demyelination and loss of oligodendrocyte-lineage cells. In normal-appearing white matter, where the myelin markers MBP (δ = +0.051), PLP1 (−0.074), MOG (−0.257) and MAG (−0.299) were unchanged, *TRPC1* was lower at sample level (δ = −0.482, p = 0.005), although not after Benjamini–Hochberg correction (q = 0.10); because several samples come from the same donor, we repeated the test with the donor as the unit of inference, which gave a larger difference (seven multiple sclerosis donors versus five control donors, δ = −0.771, uncorrected p = 0.030). Regression of TRPC1 on MBP, PLP1 and GFAP attenuated the difference, which remained significant at sample level (δ = −0.562 → −0.418, p = 0.002) but not at donor level (δ = −0.640, p = 0.055). *TRPC1* was lower in every lesion type, including remyelinating and inactive lesions (donor-level δ = −1.000 for both). Multiple sclerosis showed marked astrogliosis (GFAP δ = +0.84 to +0.96), again with *TRPC1* falling (Supplementary Table S16).

Among the disease cohorts and regions surveyed here, ALS was the only setting in which *TRPC1* was significantly increased.

### 3.10 A junction-based marker of TDP-43 dysfunction does not correlate with TRPC1 within ALS regions

The comparisons above establish that *TRPC1*, *SARAF* and *CBARP* differ between ALS and control tissue. They do not establish that the differences are caused by TDP-43 loss. A per-sample marker of TDP-43 dysfunction is needed for that question; gene-level *STMN2*, the proxy used in earlier work, is inadequate on its own because in bulk tissue it also strongly tracks neuronal content.

The same cohort is present in recount3 as a junction count matrix, which allowed the truncated *STMN2* transcript to be quantified directly (Methods 2.11). This junction-based measure is more proximal to a known consequence of TDP-43 dysfunction than total *STMN2* expression, but remains sparse and can capture only one route of TDP-43-dependent RNA misprocessing. Its regional distribution is consistent with the expected anatomical pattern (Figure 9A): 70.5% of lumbar-cord and 63.5% of cervical-cord ALS samples carried the cryptic junction, compared with 33.7–43.5% across motor cortex, 28.3% in temporal cortex, 11.2% in frontal cortex, and **2.6% in cerebellum**. Against non-neurological controls from the same cohort, cryptic inclusion was increased in lumbar cord (Cliff's δ = +0.705, q = 9.6 × 10⁻¹⁴), cervical cord (+0.619, q = 1.9 × 10⁻¹⁰), medial motor cortex (+0.337, q = 0.012) and temporal cortex (+0.283, q = 1.6 × 10⁻³), and was unchanged in cerebellum (−0.040, q = 0.34; Supplementary Table S8).

The junction-based measure also helps assess the proxy problem. Within ALS samples, gene-level *STMN2* correlated strongly with the neuronal marker *SNAP25* (ρ = 0.66 to 0.93 across regions), whereas cryptic PSI did not (ρ = −0.28 to +0.05; Figure 9B). Thus, in this bulk-tissue cohort, total *STMN2* expression should not be interpreted as a specific measure of TDP-43 dysfunction; its prior positive correlations may substantially reflect neuronal content.

Using this junction-based marker, we found no significant correlation between cryptic *STMN2* inclusion and *TRPC1* in any region (ρ = −0.19 to +0.18; all q > 0.05). Of the 110 correlations tested against the junction-based marker, eight survived Benjamini–Hochberg correction, applied across the 220 informative correlations computed with the two proxies (the eleven self-correlations of gene-level *STMN2* were excluded from the family; Supplementary Table S9): gene-level *STMN2* fell with cryptic inclusion (cervical cord ρ = −0.244, q = 0.011), *SNAP25* fell and *GFAP* rose in medial motor cortex (ρ = −0.276 and +0.284, q = 0.028 and 0.023) — all three internally consistent — together with *ATP2A2* (medial motor cortex, −0.299, q = 0.015), *STIM1* (lumbar cord, −0.280, q = 0.0021), *SARAF* (lumbar cord, −0.194, q = 0.049) and *ORAI1* (+0.278 in cervical and +0.234 in lumbar cord, q = 0.0031 and 0.013). Of these five correlations in four SOCE genes, three run in the direction the cellular model predicts (the *ATP2A2* fall and both *ORAI1* rises) and two run against it (*STIM1* and *SARAF*).

We therefore report the patient-tissue findings as association with ALS and not as evidence that TDP-43 loss drives them. The regional pattern makes the point sharply: cryptic *STMN2* is highest in spinal cord, where *TRPC1* does not increase, and *TRPC1* rises in cortex, where cryptic inclusion is comparatively rare.

## 4. Discussion

This study asked whether TDP-43 knockdown is associated with altered store-operated Ca²⁺ entry and which RNA-level changes might explain that phenotype. The data support a marked reduction in the Fura-2 SOCE readout, but do not establish its molecular level or a causal RNA mechanism.

**A functional deficit with candidate transcript explanations.** TDP-43 knockdown reduced the Fura-2 SOCE response by approximately 84%, while the ER-release comparison was inconclusive at the current sample size. This supports a cellular phenotype, but three cultures per group is a small sample, the result has not yet been reproduced with a second silencing reagent, and the measurement does not localise the defect exclusively to plasma-membrane influx. After adjustment for library composition, the transcript summaries show a larger ORAI pool whose composition shifts from *ORAI2* towards *ORAI3*, modest decreases in the dominant *ORAI2* and *ATP2A2* transcripts, and an increase in *SARAF*. The *ORAI3* and *SARAF* increases are the changes most directly compatible with reduced entry (Palty et al., 2012; Yoast et al., 2020), whereas the *ORAI2* decrease would be expected to act in the opposite direction (Vaeth et al., 2017). These patterns are compatible with a compositional contribution, but TPM is not an absolute measurement and protein abundance, membrane localisation and complex assembly were not measured. ORAI2 and MCUB also have context-dependent regulatory roles, so none of these RNA changes can be assigned a functional direction in these cells without rescue or perturbation experiments.

**RNA-processing analyses are informative but not exclusionary.** Annotation-free junction analysis recovered canonical TDP-43 cryptic targets and none of those controls in FUS or TAF15 knockdown. In the coverage-based APA analysis the *STMN2* positive control shifted only moderately in the primary model but clearly in the iPSC-derived motor neurons (intronic index +0.249), where the SOCE-machinery genes nonetheless moved little; the conservative NMD analysis did not yield a genome-wide significant gene. These data therefore do not establish direct or indirect targeting of the SOCE machinery by any one RNA-processing route. *CBARP* remains a reproducible splicing candidate, but its NMD sensitivity is exploratory and requires isoform-level validation.

**Splicing changes are real but individually small.** Robust, reproducible splicing changes do occur in SOCE regulators. *CBARP* is the clearest: the same locus is affected in five of six datasets across two species, with a median |ΔPSI| of 0.37 and a median of 246 junction reads per event, independently confirmed by an annotation-free method, and accompanied by the largest expression change of any SOCE regulator — although the direction of the splicing change differs between the human cell line and the mouse lines, so it is the involvement of the locus rather than a single directional switch that reproduces. *STIM2*, *STIMATE* and *ORAI3* are supported in the primary model. What we could not demonstrate is that any of these individually accounts for the functional deficit, and the honest reading of the data is that no single event does.

We were unable to confirm a specific mechanistic hypothesis we had considered attractive. Inclusion of the 24-nucleotide exon that converts STIM2 into the inhibitory STIM2.1/STIM2β isoform (Miederer et al., 2015; Rana et al., 2015) was higher in knockdown in five of six datasets, but with a pooled ΔPSI of +0.0013 (95% CI −0.022 to +0.024; p = 0.914) the effect size is negligible. Within the sensitivity of these data, the results do not support a STIM2.1 isoform switch as the explanation for the SOCE reduction.

**Methodological implications.** Two of our results are cautionary and, we think, generalisable. First, threshold-based splicing calls at the |ΔPSI| ≥ 0.10 level are not interpretable without read-level support: at the coverage typical of these datasets, power at that effect size is below 0.30, so significant calls are enriched for inflated estimates. A *TRPC1* event that passed the conventional thresholds failed every robustness check we applied to it (Figure 2). Second, enrichment tests for splicing changes in curated gene panels require a covariate-matched null. Ca²⁺ channel genes are long, multi-exonic and highly expressed, and against a naive background this alone produces apparent enrichment; matching removed the signal in two of the three human models and left only a nominal signal, uncorrected for the number of tests, in iPSC-derived motor neurons.

**Relevance to disease tissue.** *TRPC1*, *SARAF* and *CBARP* differ between ALS and control samples in several regions, and *TRPC1* shows the opposite direction in the comparison cohorts; the Parkinson's disease decrease agrees with the lower TRPC1 levels reported in brain lysates from patients (Selvaraj et al., 2012). These are disease-associated transcript patterns. Bulk tissue, region and donor differences prevent a causal attribution to TDP-43 loss; the junction-level *STMN2* analysis did not provide a within-ALS correlation with *TRPC1*.

Two features of these data restrain interpretation. The regional pattern is discordant with ALS pathology: *TRPC1* rises in cortex, hippocampus and cerebellum but not in spinal cord, where motor neuron loss is most severe. Bulk spinal cord in advanced ALS is dominated by glia, and the cortical signal may reflect tissue in which vulnerable neurons are still present; but this remains an explanation to be tested rather than one supported here.

More importantly, *TRPC1* expression does not scale with the junction-based marker of TDP-43 dysfunction within ALS cases. With gene-level *STMN2* as the proxy, *TRPC1* correlates positively in ten of the eleven regions (seven with q < 0.05; Supplementary Table S9), the opposite of what the mechanistic model predicts. Gene-level *STMN2* correlates with *SNAP25* at up to ρ = 0.93 within regions and should not be treated as a specific measure of TDP-43 function in bulk tissue. Cryptic *STMN2* inclusion is a more direct but still incomplete marker; it lacks this strong correlation with *SNAP25*, is sparse outside spinal cord and reflects one RNA-processing route. Measured with this limitation in mind, *TRPC1* showed no significant relationship to cryptic *STMN2* inclusion in any region. Of the SOCE components tested, the *ATP2A2* fall in medial motor cortex and the *ORAI1* rise in cord move in the direction the cellular model predicts, while *STIM1* and *SARAF* move against it.

The patient-tissue findings are therefore an association with ALS and are not offered in support of causation. The regional dissociation is the clearest statement of the problem: cryptic *STMN2* is most abundant in spinal cord, where *TRPC1* does not change, and *TRPC1* rises in cortical regions where cryptic inclusion is comparatively rare. Whether the cortical *TRPC1* increase reflects a slower, indirect consequence of TDP-43 dysfunction, a cell-composition shift that our markers do not capture, or a process unrelated to TDP-43, cannot be settled with bulk data.

**Relation to existing work.** Previous work linking TDP-43 to Ca²⁺ has centred on ER–mitochondrial contacts, where disease-associated TDP-43 disrupts VAPB–PTPIP51 tethering and reduces IP₃-receptor-mediated mitochondrial Ca²⁺ uptake (Stoica et al., 2014). Our findings concern a different route — plasma-membrane store-operated influx — and the two are not mutually exclusive. The reduction we observe in mitochondrial uptake transcripts is compatible with convergent effects on the same Ca²⁺ network, but does not establish a change in mitochondrial flux. The context dependence of ORAI2 regulation (Vaeth et al., 2017) and the reported SARAF activity in SH-SY5Y (Albarran et al., 2016) also argue for direct perturbation experiments rather than directional inference from transcript abundance.

## 5. Limitations

SH-SY5Y is a neuroblastoma line and not a motor neuron. The functional measurements rest on three independent cultures per group for the Ca²⁺ measurements and four wells for viability; the number of biological replicates is small, and the phenotype requires replication with a second silencing reagent and a TDP-43 rescue. TDP-43 knockdown was confirmed at the mRNA level only; TDP-43 protein, and the protein levels, membrane localisation and complex assembly of candidate genes, were not measured. Relative expression rests on a single reference gene (GAPDH) and on the assumption of the 2^−ΔΔCt^ method that target and reference genes amplify with equal efficiency. The WST-1 assay does not separate proliferation from cell death, and the opposite directions at 24 h and 48 h show that it is not a direct viability measure. The experiment was performed three times, but the statistics reported compare the four wells of one experiment, so the very small p values reflect well-to-well rather than between-experiment variation; the reduced signal at 48 h also leaves open whether impaired cell health contributes to the lower SOCE response. Transcript-family summaries are length- and composition-adjusted short-read estimates, not absolute stoichiometry.

In the primary RNA-seq model, knockdown is induced with doxycycline and compared with untreated cells of the same line; GSE296712 contains no doxycycline-treated control without TDP-43 knockdown. Doxycycline inhibits mitochondrial translation and can alter nuclear gene expression and metabolism in cultured cells (Moullan et al., 2015), so part of the SH-SY5Y transcript response may reflect the inducer rather than TDP-43 loss, and the mitochondrial Ca²⁺-uptake family is the most exposed to this confound. The laboratory knockdown uses no doxycycline and is compared with a non-targeting shRNA control, and the four RT-qPCR targets changed there in the same direction as in the RNA-seq data; the iPSC-colony comparison is likewise made against a control knockdown.

Patient-tissue analyses are bulk and associative. Cryptic *STMN2* PSI is a more direct but incomplete marker of one consequence of TDP-43 dysfunction, and it is sparse: median inclusion is zero outside the spinal cord, so the correlation analysis has good power in cord and progressively less in cortical regions, and a small true effect in cortex would not be detected. Two further limits apply to the negative result. The measure reports the state of the tissue at death rather than the cumulative history of TDP-43 dysfunction that shaped its transcriptome, and a downstream change that persists after the lesion would not correlate with it. Truncated *STMN2* is also generated by cryptic polyadenylation, so its junction reports one RNA-processing route and may not represent all of TDP-43 function.

The ratio of null calls to real calls in the annotation-free analysis varies by dataset and, in two of the three datasets where a split-control null could be computed, approaches or exceeds one. That ratio is not a calibrated false-discovery rate; calibration would require label permutation and depth-matched simulation, which we did not perform. Call counts are therefore not comparable between datasets, and our conclusions rest on positive-control recovery and on the identity of the genes called rather than on how many were called. In the remaining eight comparisons the control replicates could not be split, so no null could be computed there. Junctions were extracted with strand information but aggregated by chromosome and coordinate, so at loci with an antisense partner a junction could in principle be assigned to either gene. Several SOCE genes have such a partner, *CBARP* among them: its junctions fall inside the divergent transcript *CBARP-DT*. For *CBARP* the same event is independently supported by the strand-aware rMATS and LeafCutter analyses, but transcription on the antisense strand also changes — *CBARP-DT* itself showed significant isoform switching in SH-SY5Y (IsoformSwitchAnalyzeR gene-level q = 0.008) — so for any single junction-level call at an antisense locus the gene assignment should be treated as provisional.

The alternative-polyadenylation analysis is coverage-based and does not localise poly(A) sites; 3′-end sequencing would be required for that purpose. In SH-SY5Y the *STMN2* positive control shifted only moderately (+0.145, bootstrap 95% CI +0.085 to +0.205) and falls just below the pre-specified depth filter, so the negative result in that model rests on an assay whose sensitivity could not be demonstrated there. Of the four comparisons in which the analysis could be run, the assay's own positive control is informative only in the iPSC-derived motor neurons, since it is attenuated in SH-SY5Y and the two mouse lines have no conserved control. The iPSC-derived motor neurons have only two replicates per group, so the interval around their working control is coarse. The analysis identified one moderate candidate gradient in SH-SY5Y, in *STIM1* intron 17, which lies next to an alternatively included *STIM1* exon and was not confirmed in another comparison, and no large APA shift in the SOCE machinery in any comparison. Sparse windows remain a sensitivity limit, so the negative result is not definitive evidence of sparing.

FRASER was applied but is not informative at this cohort size, because its model treats aberrant splicing as rare within the cohort and the depleted samples are the majority. MAJIQ (Vaquero-Garcia et al., 2016) itself was not run; the local-splicing-variation analysis we substituted covers the same annotation-free ground but does not reproduce MAJIQ's quantification model.

The nonsense-mediated decay experiment is a published gene-level dataset in a different cell model. Shared control libraries were reused across the four intervention contrasts; after treating the four conditions as the unit of inference, no gene passed genome-wide FDR correction. Because the variance between the four conditions does not capture the uncertainty of the shared control difference, and because the four interventions are distinct perturbations rather than biological replicates, even the conservative q values should not be read as fully calibrated. The analysis tests whether a gene-level transcript pool is NMD-sensitive, not whether the specific cryptic isoform identified here is the degraded species.

## 6. Conclusion

TDP-43 depletion in SH-SY5Y cells was associated with a large reduction in Fura-2-measured SOCE and with altered expression and splicing of calcium-regulatory genes. Transcript-family summaries and reproducible *CBARP* splicing identify candidate explanations, but they do not establish protein stoichiometry or a single causal RNA switch. Coverage-based APA and conservative NMD analyses do not provide confirmatory evidence that SOCE-machinery genes are direct targets of those pathways. In ALS tissue, *TRPC1*, *SARAF* and *CBARP* show disease-associated expression differences, whereas junction-level *STMN2* inclusion does not correlate with *TRPC1* within ALS samples. *CBARP* and *SARAF* are therefore priorities for targeted molecular and functional follow-up.

## References

Albarran L, Lopez JJ, Woodard GE, Salido GM, Rosado JA. Store-operated Ca²⁺ entry-associated regulatory factor (SARAF) plays an important role in the regulation of arachidonate-regulated Ca²⁺ (ARC) channels. *J Biol Chem*. 2016;291:6982–6988. doi:10.1074/jbc.M115.704940

Ambudkar IS, Ong HL, Liu X, Bandyopadhyay BC, Cheng KT. TRPC1: the link between functionally distinct store-operated calcium channels. *Cell Calcium*. 2007;42:213–223. doi:10.1016/j.ceca.2007.01.013

Andrews S. FastQC: a quality control tool for high throughput sequence data. Babraham Bioinformatics; 2010. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Arai T, Hasegawa M, Akiyama H, et al. TDP-43 is a component of ubiquitin-positive tau-negative inclusions in frontotemporal lobar degeneration and amyotrophic lateral sclerosis. *Biochem Biophys Res Commun*. 2006;351:602–611. doi:10.1016/j.bbrc.2006.10.093

Baughn MW, Melamed Z, López-Erauskin J, et al. Mechanism of STMN2 cryptic splice-polyadenylation and its correction for TDP-43 proteinopathies. *Science*. 2023;379:1140–1149. doi:10.1126/science.abq5622

Béguin P, Nagashima K, Mahalakshmi RN, et al. BARP suppresses voltage-gated calcium channel activity and Ca²⁺-evoked exocytosis. *J Cell Biol*. 2014;205:233–249. doi:10.1083/jcb.201304101

Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B*. 1995;57:289–300.

Brown AL, Wilkins OG, Keuss MJ, et al. TDP-43 loss and ALS-risk SNPs drive mis-splicing and depletion of UNC13A. *Nature*. 2022;603:131–137. doi:10.1038/s41586-022-04436-3

Bryce-Smith S, Brown AL, Chien MZYJ, et al. TDP-43 loss induces cryptic polyadenylation in ALS/FTD. *Nat Neurosci*. 2025;28:2190–2200. doi:10.1038/s41593-025-02050-w

Cliff N. Dominance statistics: ordinal analyses to answer ordinal questions. *Psychol Bull*. 1993;114:494–509. doi:10.1037/0033-2909.114.3.494

Cotto KC, Feng YY, Ramu A, et al. Integrated analysis of genomic and transcriptomic data for the discovery of splice-associated variants in cancer. *Nat Commun*. 2023;14:1589. doi:10.1038/s41467-023-37266-6

Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *Gigascience*. 2021;10:giab008. doi:10.1093/gigascience/giab008

Dumitriu A, Golji J, Labadorf AT, et al. Integrative analyses of proteomics and RNA transcriptomics implicate mitochondrial processes, protein folding pathways and GWAS loci in Parkinson disease. *BMC Med Genomics*. 2016;9:5. doi:10.1186/s12920-016-0164-y

Elkjaer ML, Frisch T, Reynolds R, et al. Molecular signature of different lesion types in the brain white matter of patients with progressive multiple sclerosis. *Acta Neuropathol Commun*. 2019;7:205. doi:10.1186/s40478-019-0855-7

Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*. 2016;32:3047–3048. doi:10.1093/bioinformatics/btw354

Frankish A, Diekhans M, Ferreira A-M, et al. GENCODE reference annotation for the human and mouse genomes. *Nucleic Acids Res*. 2019;47:D766–D773. doi:10.1093/nar/gky955

Grosskreutz J, Van Den Bosch L, Keller BU. Calcium dysregulation in amyotrophic lateral sclerosis. *Cell Calcium*. 2010;47:165–174. doi:10.1016/j.ceca.2009.12.002

Grynkiewicz G, Poenie M, Tsien RY. A new generation of Ca²⁺ indicators with greatly improved fluorescence properties. *J Biol Chem*. 1985;260:3440–3450.

Hoover PJ, Lewis RS. Stoichiometric requirements for trapping and gating of Ca²⁺ release-activated Ca²⁺ (CRAC) channels by stromal interaction molecule 1 (STIM1). *Proc Natl Acad Sci USA*. 2011;108:13299–13304. doi:10.1073/pnas.1101664108

Hruska-Plochan M, Wiersma VI, Betz KM, et al. A model of human neural networks reveals NPTX2 pathology in ALS and FTLD. *Nature*. 2024;626:1073–1083. doi:10.1038/s41586-024-07042-7

Jing J, He L, Sun A, et al. Proteomic mapping of ER–PM junctions identifies STIMATE as a regulator of Ca²⁺ influx. *Nat Cell Biol*. 2015;17:1339–1347. doi:10.1038/ncb3234

Kapeli K, Pratt GA, Vu AQ, et al. Distinct and shared functions of ALS-associated proteins TDP-43, FUS and TAF15 revealed by multisystem analyses. *Nat Commun*. 2016;7:12143. doi:10.1038/ncomms12143

Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nat Biotechnol*. 2019;37:907–915. doi:10.1038/s41587-019-0201-4

Klim JR, Williams LA, Limone F, et al. ALS-implicated protein TDP-43 sustains levels of STMN2, a mediator of motor neuron growth and repair. *Nat Neurosci*. 2019;22:167–179. doi:10.1038/s41593-018-0300-4

Lambert JP, Luongo TS, Tomar D, et al. MCUB regulates the molecular composition of the mitochondrial calcium uniporter channel to limit mitochondrial calcium overload during stress. *Circulation*. 2019;140:1720–1733. doi:10.1161/CIRCULATIONAHA.118.037968

Leek JT. svaseq: removing batch effects and other unwanted noise from sequencing data. *Nucleic Acids Res*. 2014;42:e161. doi:10.1093/nar/gku864

Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. The Sequence Read Archive. *Nucleic Acids Res*. 2011;39:D19–D21. doi:10.1093/nar/gkq1019

Li YI, Knowles DA, Humphrey J, et al. Annotation-free quantification of RNA splicing using LeafCutter. *Nat Genet*. 2018;50:151–158. doi:10.1038/s41588-017-0004-9

Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*. 2014;30:923–930. doi:10.1093/bioinformatics/btt656

Ling JP, Pletnikova O, Troncoso JC, Wong PC. TDP-43 repression of nonconserved cryptic exons is compromised in ALS-FTD. *Science*. 2015;349:650–655. doi:10.1126/science.aab0983

Ling S-C, Polymenidou M, Cleveland DW. Converging mechanisms in ALS and FTD: disrupted RNA and protein homeostasis. *Neuron*. 2013;79:416–438. doi:10.1016/j.neuron.2013.07.033

Livak KJ, Schmittgen TD. Analysis of relative gene expression data using real-time quantitative PCR and the 2^−ΔΔCt^ method. *Methods*. 2001;25:402–408. doi:10.1006/meth.2001.1262

Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol*. 2014;15:550. doi:10.1186/s13059-014-0550-8

Ma XR, Prudencio M, Koike Y, et al. TDP-43 represses cryptic exon inclusion in the FTD–ALS gene UNC13A. *Nature*. 2022;603:124–130. doi:10.1038/s41586-022-04424-7

Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet J*. 2011;17:10–12. doi:10.14806/ej.17.1.200

Melamed Z, López-Erauskin J, Baughn MW, et al. Premature polyadenylation-mediated loss of stathmin-2 is a hallmark of TDP-43-dependent neurodegeneration. *Nat Neurosci*. 2019;22:180–190. doi:10.1038/s41593-018-0293-z

Mertes C, Scheller IF, Yépez VA, et al. Detection of aberrant splicing events in RNA-seq data using FRASER. *Nat Commun*. 2021;12:529. doi:10.1038/s41467-020-20573-7

Miederer AM, Alansary D, Schwär G, et al. A STIM2 splice variant negatively regulates store-operated calcium entry. *Nat Commun*. 2015;6:6899. doi:10.1038/ncomms7899

Moullan N, Mouchiroud L, Wang X, et al. Tetracyclines disturb mitochondrial function across eukaryotic models: a call for caution in biomedical research. *Cell Rep*. 2015;10:1681–1691. doi:10.1016/j.celrep.2015.02.034

Nelson PT, Dickson DW, Trojanowski JQ, et al. Limbic-predominant age-related TDP-43 encephalopathy (LATE): consensus working group report. *Brain*. 2019;142:1503–1527. doi:10.1093/brain/awz099

Neumann M, Sampathu DM, Kwong LK, et al. Ubiquitinated TDP-43 in frontotemporal lobar degeneration and amyotrophic lateral sclerosis. *Science*. 2006;314:130–133. doi:10.1126/science.1134108

Nowicka M, Robinson MD. DRIMSeq: a Dirichlet-multinomial framework for multivariate count outcomes in genomics. *F1000Res*. 2016;5:1356. doi:10.12688/f1000research.8900.2

Paez PM, Fulton D, Spreuer V, Handley V, Campagnoni AT. Modulation of canonical transient receptor potential channel 1 in the proliferation of oligodendrocyte precursor cells by the golli products of the myelin basic protein gene. *J Neurosci*. 2011;31:3625–3637. doi:10.1523/JNEUROSCI.4424-10.2011

Palty R, Raveh A, Kaminsky I, Meller R, Reuveny E. SARAF inactivates the store operated calcium entry machinery to prevent excess calcium refilling. *Cell*. 2012;149:425–438. doi:10.1016/j.cell.2012.01.055

Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. *Nat Methods*. 2017;14:417–419. doi:10.1038/nmeth.4197

Polymenidou M, Lagier-Tourenne C, Hutt KR, et al. Long pre-mRNA depletion and RNA missplicing contribute to neuronal vulnerability from loss of TDP-43. *Nat Neurosci*. 2011;14:459–468. doi:10.1038/nn.2779

Prakriya M, Lewis RS. Store-operated calcium channels. *Physiol Rev*. 2015;95:1383–1436. doi:10.1152/physrev.00020.2014

Prudencio M, Humphrey J, Pickles S, et al. Truncated stathmin-2 is a marker of TDP-43 pathology in frontotemporal dementia. *J Clin Invest*. 2020;130:6080–6092. doi:10.1172/JCI139741

Putney JW. Pharmacology of store-operated calcium channels. *Mol Interv*. 2010;10:209–218. doi:10.1124/mi.10.4.4

Rana A, Yen M, Sadaghiani AM, et al. Alternative splicing converts STIM2 from an activator to an inhibitor of store-operated calcium channels. *J Cell Biol*. 2015;209:653–670. doi:10.1083/jcb.201412060

Selli Ç, Eraç Y, Kosova B, Tosun M. Post-transcriptional silencing of TRPC1 ion channel gene by RNA interference upregulates TRPC6 expression and store-operated Ca²⁺ entry in A7r5 vascular smooth muscle cells. *Vascul Pharmacol*. 2009;51:96–100. doi:10.1016/j.vph.2009.04.001

Selvaraj S, Sun Y, Watt JA, et al. Neurotoxin-induced ER stress in mouse dopaminergic neurons involves downregulation of TRPC1 and inhibition of AKT/mTOR signaling. *J Clin Invest*. 2012;122:1354–1367. doi:10.1172/JCI61332

Shen S, Park JW, Lu Z-x, et al. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. *Proc Natl Acad Sci USA*. 2014;111:E5593–E5601. doi:10.1073/pnas.1419161111

Sinha IR, Ye Y, Li Y, et al. Inhibition of nonsense-mediated decay in TDP-43 deficient neurons reveals novel cryptic exons. *bioRxiv* [preprint]. 2025. doi:10.1101/2025.06.28.661837

Srinivasan K, Friedman BA, Etxeberria A, et al. Alzheimer's patient microglia exhibit enhanced aging and unique transcriptional activation. *Cell Rep*. 2020;31:107843. doi:10.1016/j.celrep.2020.107843

Stoica R, De Vos KJ, Paillusson S, et al. ER–mitochondria associations are regulated by the VAPB–PTPIP51 interaction and are disrupted by ALS/FTD-associated TDP-43. *Nat Commun*. 2014;5:3996. doi:10.1038/ncomms4996

Šušnjar U, Škrabar N, Brown AL, et al. Cell environment shapes TDP-43 function with implications in neuronal and muscle disease. *Commun Biol*. 2022;5:314. doi:10.1038/s42003-022-03253-8

Trincado JL, Entizne JC, Hysenaj G, et al. SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions. *Genome Biol*. 2018;19:40. doi:10.1186/s13059-018-1417-1

Vaeth M, Yang J, Yamashita M, et al. ORAI2 modulates store-operated calcium entry and T cell-mediated immunity. *Nat Commun*. 2017;8:14714. doi:10.1038/ncomms14714

Van den Berge K, Soneson C, Robinson MD, Clement L. stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. *Genome Biol*. 2017;18:151. doi:10.1186/s13059-017-1277-0

Van Nostrand EL, Freese P, Pratt GA, et al. A large-scale binding and functional map of human RNA-binding proteins. *Nature*. 2020;583:711–719. doi:10.1038/s41586-020-2077-3

Vaquero-Garcia J, Barrera A, Gazzara MR, et al. A new view of transcriptome complexity and regulation through the lens of local splicing variations. *eLife*. 2016;5:e11752. doi:10.7554/eLife.11752

Vitting-Seerup K, Sandelin A. IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. *Bioinformatics*. 2019;35:4469–4471. doi:10.1093/bioinformatics/btz247

Voskuhl RR, Itoh N, Tassoni A, et al. Gene expression in oligodendrocytes during remyelination reveals cholesterol homeostasis as a therapeutic target in multiple sclerosis. *Proc Natl Acad Sci USA*. 2019;116:10130–10139. doi:10.1073/pnas.1821306116

Wilks C, Zheng SC, Chen FY, et al. recount3: summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biol*. 2021;22:323. doi:10.1186/s13059-021-02533-6

Xia Z, Donehower LA, Cooper TA, et al. Dynamic analyses of alternative polyadenylation from RNA-seq reveal a 3′-UTR landscape across seven tumour types. *Nat Commun*. 2014;5:5274. doi:10.1038/ncomms6274

Yoast RE, Emrich SM, Zhang X, et al. The native ORAI channel trio underlies the diversity of Ca²⁺ signaling events. *Nat Commun*. 2020;11:2444. doi:10.1038/s41467-020-16232-6

## Declarations

**Ethics approval.** No new human or animal material was collected for this study. All RNA-seq datasets are publicly available and were reanalysed under the terms of the repositories that host them. The laboratory experiments used an established cell line (SH-SY5Y) and no human participants or animals. Lentiviral work was carried out under the institutional biosafety rules of Ege University.

**Funding.** This work was supported by the Ege University Scientific Research Projects Coordination Unit (BAP), project 31901. E.Y. was supported by the Council of Higher Education (YÖK) 100/2000 PhD Scholarship Programme (priority area: Molecular Pharmacology and Drug Research) and by the TÜBİTAK BİDEB 2211-A National PhD Scholarship Programme. The funders had no role in study design, data collection and analysis, the decision to publish, or the preparation of the manuscript.

**Author contributions (CRediT).** Elmasnur Yılmaz: Investigation, Formal analysis, Software, Visualization, Writing – original draft. Yasemin Eraç: Conceptualization, Methodology, Resources, Funding acquisition, Supervision, Writing – review & editing. Both authors read and approved the final manuscript.

**Conflicts of interest.** The authors declare no competing interests.

**Acknowledgements.** This article is derived from the doctoral thesis of Elmasnur Yılmaz (Ege University, Graduate School of Natural and Applied Sciences, Department of Biotechnology, 2026). The authors thank the ALS Consortium of the New York Genome Center and Target ALS for making the post-mortem RNA-seq data available, the ENCODE Consortium for the K562 knockdown data, and the authors of GSE296712, GSE230647, GSE77702, GSE27394, GSE171714, GSE307054, GSE125583, GSE68719, GSE138614 and GSE123496 for depositing their raw data.

## Figures

**Figure 1. Detection power for splicing events in a 3 + 3 design,** simulated over the observed coverage distribution for true ΔPSI values of 0.05 to 0.30 at mean informative depths of 10 to 100 reads per sample.

**Figure 2. Threshold-based splicing calls require read-level verification.** (A) Replicate-level PSI for the *TRPC1* skipped-exon event, with IJC/SJC counts annotated. (B) Bootstrap distribution of ΔPSI (20,000 resamples); the 95% interval spans zero. (C) Leave-one-out ΔPSI; removal of one control replicate reverses the sign.

**Figure 3. SOCE splicing candidates in SH-SY5Y** with replicate-level bootstrap 95% confidence intervals (10,000 resamples). Green marks intervals that exclude zero; the *STIM1* event (grey) and the *TRPC1* event of Figure 2 (red; 20,000 resamples) are shown for comparison.

**Figure 4. Annotation-free junction analysis and exploratory NMD interaction.** (A) Cryptic junctions of the literature positive controls recovered de novo in SH-SY5Y (75 ng/mL doxycycline, mapping-quality-filtered junction set), with knockdown and control read counts. (B) Number of genes in each of the four cumulative Ca²⁺ panels that carry a high-confidence unannotated splicing change in the same comparison. (C) Condition-level TDP-43-specific NMD interaction for genes carrying cryptic junctions versus the rest of the transcriptome, with *CBARP* marked; the NMD analysis is exploratory because the TDP-43 factor is confounded with sequencing batch.

**Figure 5. Specificity of the cryptic calls and the coverage-based alternative-polyadenylation analysis.** (A) Positive-control genes recovered in each human comparison under the permissive definition, with FUS and TAF15 knockdown as controls; the mouse comparisons are omitted because the literature controls are human cryptic events that are not conserved in mouse. (B) Depth-qualified intronic and distal 3′-terminal usage estimates in SH-SY5Y; the *STIM2* intron 13 unit (†) contains the alternatively spliced *STIM2* exon of Table 2 and reports that splicing change. (C) Coverage at the two ends of *STMN2* intron 2 (the canonical intron 1, which contains cryptic exon 2a) in each replicate after reference skips are excluded; one knockdown replicate has a measured coverage of zero in the 3′ window.

**Figure 6. Functional consequences of TDP-43 knockdown in SH-SY5Y cells.** (A) TARDBP mRNA in non-transduced cells, non-targeting shRNA controls and shTDP-43 cells (2^−ΔΔCt^, n = 4; one-way ANOVA on log2 values with Tukey's test). (B) Relative mRNA of the four SOCE-associated targets in shTDP-43 cells, relative to the non-targeting shRNA control (n = 4; two-tailed t-tests). (C) WST-1 viability at 24 h and 48 h, as a percentage of the non-targeting shRNA control (n = 4 wells from one of three independent experiments). (D) Representative Fura-2 traces for a non-targeting shRNA control and a shTDP-43 sample: SERCA inhibition with 10 µM cyclopiazonic acid in Ca²⁺-free medium, then re-addition of 1.5 mM CaCl₂; the two traces are shown on different vertical scales and time windows, and the group comparison is in (E). (E) Group data for the two phases, as Δ(F340/F380) relative to the preceding baseline (n = 3 independent cultures per group). Bars are mean ± SEM with individual values; ** p < 0.01, *** p < 0.001, **** p < 0.0001.

**Figure 7. Transcript-family abundance of SOCE-related genes in SH-SY5Y.** Gene-level TPM, adjusted for library composition (Methods 2.10), per family member in control and TDP-43 knockdown cells for the ORAI, SERCA, STIM, SOCE-regulator and mitochondrial-uptake families, with the net change of each family pool; the TRPC and PMCA families are given in Table 4. Asterisks mark the RT-qPCR targets. These are transcript-level summaries, not protein stoichiometry.

**Figure 8. Direction of *TRPC1* expression change across five diseases.** Cliff's δ (case − control) for ALS (NYGC; all seven brain regions and three spinal cord levels), other neurological disorders (same controls, three testable regions), Alzheimer's disease, Parkinson's disease and multiple sclerosis, the last including the donor-level, lesion, normal-appearing white matter (sample and donor level) and myelin-adjusted analyses. Asterisks mark Benjamini–Hochberg q < 0.05 or, for the donor-level normal-appearing white matter and myelin-adjusted tests, uncorrected p < 0.05; pale bars are not significant.

**Figure 9. Junction-level cryptic *STMN2* in ALS post-mortem tissue.** (A) Proportion of samples carrying the cryptic junction by region, ALS versus non-neurological control. (B) Spearman correlations of the same target genes against two proxies for TDP-43 loss in the same samples: gene-level *STMN2* and cryptic *STMN2* PSI.

## Supplementary

**S1.** Laboratory source data: raw Ct values and per-replicate relative expression, Fura-2 amplitudes, WST-1 viability values, primer sequences and thermal profile, and the summary statistics behind every panel of Figure 6.

**S2.** Ca²⁺ gene panels (four cumulative sets).

**S3.** rMATS events meeting FDR < 0.05 and |ΔPSI| ≥ 0.10, JC and JCEC, six datasets, with the raw junction counts needed to reproduce the coverage pre-filter.

**S4.** Matched permutation enrichment results, all datasets and panels.

**S5.** High-confidence unannotated splicing changes in every comparison, including the independently produced mapping-quality-filtered SH-SY5Y junction set.

**S6.** Cryptic positive controls, the sixteen literature genes, and the dataset × gene recovery matrix for the human comparisons.

**S7.** SOCE genes in the annotation-free analysis.

**S8.** Cryptic *STMN2* PSI in ALS versus control by region (NYGC).

**S9.** Correlations of both TDP-43 proxies with target genes within ALS samples; the `in_correction_family` column marks the 220 informative tests over which the Benjamini–Hochberg correction was applied.

**S10.** Nonsense-mediated decay interaction for the SOCE panel and the panel-level tests, four-condition analysis.

**S11.** Depth-qualified intronic polyadenylation and 3′UTR usage estimates in all four comparisons, with the genomic windows of every unit; these are coverage gradients, not direct poly(A)-site calls.

**S12.** Machine-readable version of Table 3: cryptic event counts by comparison, with the FUS and TAF15 knockdown controls and the null-test ratios.

**S13.** *STIM2* SOAR exon measured at junction level, with both flanking junctions, in every comparison in which it was measurable.

**S14.** Control-versus-control null test of the cryptic calling procedure, under the permissive definition and three stricter thresholds.

**S15.** The same *STIM2.1*/SOAR exon measured as an rMATS event, with per-replicate PSI and bootstrap intervals, across six datasets — the estimate compared with S13 in Section 3.3.

**S16.** Multiple sclerosis analysis, both cohorts, with myelin adjustment and the donor-level re-analysis.

**S17.** Accession list for every dataset analysed, with design, library type and sample-to-group assignment.

## Tables

**Table 1.** Effect of coverage pre-filtering on event counts and significance, by dataset.

<!-- table:1 -->

| Dataset | Events tested | Events after filter | Removed (%) | Significant before filter | Significant after filter | Significant calls lost (%) |
|:--------------------------|---------:|---------:|--------:|---------:|---------:|----------:|
| SH-SY5Y (GSE296712) | 118,664 | 90,143 | 24.0 | 7,854 | 5,282 | 32.7 |
| iPSC colonies (GSE230647) | 392,236 | 302,415 | 22.9 | 15,512 | 10,055 | 35.2 |
| iPSC-derived motor neurons (GSE77702) | 49,392 | 28,463 | 42.4 | 1,465 | 443 | 69.8 |
| Mouse striatum (GSE27394) | 23,000 | 4,969 | 78.4 | 644 | 157 | 75.6 |
| C2C12 (GSE171714) | 225,016 | 184,869 | 17.8 | 5,468 | 3,140 | 42.6 |
| NSC34 (GSE171714) | 338,017 | 273,077 | 19.2 | 8,142 | 4,599 | 43.5 |

*Note.* Significant: FDR < 0.05 and |ΔPSI| ≥ 0.10. After the filter, Benjamini–Hochberg q values were recomputed within the retained events (Methods 2.3). FDR, false discovery rate; PSI, percent spliced in.

<!-- /table:1 -->

**Table 2.** Robust splicing events in the primary SH-SY5Y model among the twelve SOCE-machinery genes examined event by event (Methods 2.3), with bootstrap confidence intervals; all four are skipped-exon events. *CBARP* was not among the twelve; its events are in Supplementary Table S3.

<!-- table:2 -->

| Gene | Exon, GRCh38 (strand) | Length (bp) | Reading frame | ΔPSI | FDR | Bootstrap 95% CI | PSI, knockdown replicates | PSI, control replicates | Mean reads per sample | Minimum reads in a sample |
|:--------------|:-------------------------------------------|:------------|-----------------:|-----------:|-----------------:|-----------------:|:-------------------|-------------------:|------------:|-------------:|
| *STIMATE* | chr3:52,895,878–52,895,955 (−) | 78 | Preserved | +0.244 | 9.0 × 10⁻⁷ | +0.095 to +0.368 | 0.359; 0.427; 0.157 | 0.138; 0.024; 0.048 | 21.7 | 14 |
| *ORAI3* | chr16:30,953,185–30,953,460 (+) | 276 | Preserved | −0.269 | 0.010 | −0.404 to −0.107 | 0.465; 0.671; 0.346 | 0.778; 0.778; 0.733 | 39.7 | 8 |
| *STIM2* | chr4:27,021,494–27,021,612 (+) | 119 | Disrupted | −0.120 | < 1 × 10⁻¹⁶ | −0.165 to −0.064 | 0.125; 0.021; 0.053 | 0.192; 0.184; 0.184 | 24.7 | 9 |
| *STIM1* | chr11:4,088,702–4,088,738 (+) | 37 | Disrupted | +0.145 | 4.0 × 10⁻⁴ | −0.002 to +0.293 | 0.088; 0.345; 0.226 | 0.172; 0.000; 0.053 | 25.5 | 9 |

*Note.* SH-SY5Y, 0 versus 75 ng/mL doxycycline, three libraries per group. ΔPSI is knockdown minus control from rMATS junction counts; the confidence interval is a replicate-level bootstrap (10,000 resamples). Coordinates are 1-based and inclusive. CI, confidence interval; FDR, false discovery rate; PSI, percent spliced in.

<!-- /table:2 -->

**Table 3.** Cryptic events recovered by annotation-free junction analysis in eleven comparisons, with RNA-binding-protein specificity controls. Positive-control recovery is not assessed in the three mouse comparisons, because the sixteen literature controls are human cryptic events that are not conserved in mouse.

<!-- table:3 -->

| Comparison | Permissive calls (genes) | Positive controls, permissive | High-confidence calls (genes) | Positive controls, high-confidence | Positive-control genes, high-confidence | Tier 1 genes | SOCE-machinery genes | Split-control null: calls (ratio) |
|:------------------|----------:|---------:|----------:|---------:|:--------------------------------|:----------|:---------|:-----------|
| SH-SY5Y 75 ng/mL | 400 (288) | 13 | 165 (113) | 13 | *ACTL6B*, *AGRN*, *ARHGAP32*, *ATG4B*, *ELAVL3*, *GPSM2*, *HDGFL2*, *KALRN*, *PFKP*, *RSF1*, *SETD5*, *STMN2*, *UNC13A* | — | — | — |
| SH-SY5Y 25 ng/mL | 351 (264) | 13 | 122 (86) | 12 | *ACTL6B*, *AGRN*, *ARHGAP32*, *ATG4B*, *ELAVL3*, *GPSM2*, *HDGFL2*, *KALRN*, *PFKP*, *SETD5*, *STMN2*, *UNC13A* | — | — | — |
| iPSC colonies | 1,561 (1,013) | 15 | 477 (326) | 12 | *ACTL6B*, *ARHGAP32*, *ATG4B*, *CAMK2B*, *ELAVL3*, *GPSM2*, *HDGFL2*, *KALRN*, *PFKP*, *SETD5*, *STMN2*, *UNC13A* | *CBARP*, *TRPM3* | *CBARP* | 307 (0.64) |
| iPSC-MN, TDP-43 KD | 141 (108) | 2 | 18 (11) | 0 | — | — | — | — |
| iPSC-MN, FUS KD | 124 (80) | 0 | 26 (15) | 0 | — | — | — | — |
| iPSC-MN, TAF15 KD | 126 (88) | 0 | 23 (13) | 0 | — | — | — | — |
| K562 total RNA | 330 (229) | 0 | 23 (17) | 0 | — | — | — | 50 (2.17) |
| K562 poly(A)+ mRNA | 1,683 (1,030) | 6 | 145 (94) | 4 | *AGRN*, *ATG4B*, *HDGFL2*, *PFKP* | — | — | — |
| C2C12 | 566 (379) | n/a | 158 (110) | n/a | n/a | — | — | — |
| NSC34 | 643 (421) | n/a | 264 (166) | n/a | n/a | — | — | — |
| Mouse striatum | 43 (32) | n/a | 12 (9) | n/a | n/a | — | — | 10 (0.83) |

*Note.* Calls are unannotated splicing changes in the regtools junction set, with the number of genes carrying them in parentheses (Methods 2.5). Positive controls are counted among the sixteen literature cryptic genes (Supplementary Table S6); n/a, not assessed in mouse. The null column gives the high-confidence calls of the control-versus-control split and their ratio to the real calls; —, fewer than four control replicates. iPSC-MN, iPSC-derived motor neurons; KD, knockdown.

<!-- /table:3 -->

**Table 4.** Transcript-family abundance of SOCE-related genes in SH-SY5Y: unadjusted and composition-adjusted TPM (Methods 2.10), each member's share of its family in control cells, and DESeq2 fold changes.

<!-- table:4 -->

| Family | Gene | TPM, control | TPM, knockdown | Adjusted TPM, control | Adjusted TPM, knockdown | Share of family, control (%) | Adjusted change (%) | log2FC | p_adj |
|:--------------------------------------|:----------------|-----------------:|------------------:|-----------------:|------------------:|--------------------:|---------------:|-----------:|----------------:|
| STIM (ER Ca²⁺ sensor) | ***STIM1*** | 35.46 | 50.54 | 30.95 | 57.74 | 64.0 | +86.5 | +0.929 | 4.3 × 10⁻⁴⁷ |
|  | *STIM2* | 20.01 | 12.02 | 17.45 | 13.73 | 36.0 | −21.3 | −0.527 | 1.8 × 10⁻⁷ |
|  | Family total | 55.47 | 62.55 | 48.40 | 71.47 | 100.0 | +47.7 | — | — |
| ORAI (CRAC channel) | *ORAI1* | 10.60 | 11.70 | 9.25 | 13.37 | 15.2 | +44.5 | +0.433 | 2.4 × 10⁻⁴ |
|  | ***ORAI2*** | 53.78 | 35.86 | 46.92 | 40.98 | 77.2 | −12.7 | −0.173 | 2.2 × 10⁻³ |
|  | *ORAI3* | 5.24 | 20.06 | 4.57 | 22.92 | 7.5 | +401.4 | +2.056 | 2.3 × 10⁻⁶⁴ |
|  | Family total | 69.62 | 67.62 | 60.75 | 77.27 | 100.0 | +27.2 | — | — |
| SERCA (Ca²⁺ re-uptake into ER) | *ATP2A1* | 0.67 | 0.96 | 0.58 | 1.10 | 0.4 | +87.8 | +0.959 | 0.024 |
|  | ***ATP2A2*** | 162.38 | 103.10 | 141.71 | 117.80 | 98.1 | −16.9 | −0.247 | 8.7 × 10⁻⁹ |
|  | *ATP2A3* | 2.40 | 7.64 | 2.09 | 8.73 | 1.5 | +317.0 | +1.306 | 1.4 × 10⁻²⁹ |
|  | Family total | 165.45 | 111.70 | 144.39 | 127.63 | 100.0 | −11.6 | — | — |
| TRPC | ***TRPC1*** | 4.34 | 6.06 | 3.79 | 6.92 | 98.2 | +82.9 | +0.958 | 1.3 × 10⁻¹² |
|  | *TRPC3* | 0.06 | 0.02 | 0.05 | 0.02 | 1.3 | — | −0.938 | 0.47 |
|  | *TRPC4* | 0.00 | 0.00 | 0.00 | 0.00 | 0.0 | — | — | — |
|  | *TRPC5* | 0.02 | 0.00 | 0.02 | 0.00 | 0.4 | — | −1.633 | 0.45 |
|  | *TRPC6* | 0.00 | 0.00 | 0.00 | 0.00 | 0.1 | — | — | — |
|  | Family total | 4.42 | 6.08 | 3.86 | 6.95 | 100.0 | +80.3 | — | — |
| SOCE regulators | ***SARAF*** | 179.27 | 191.68 | 156.45 | 218.99 | 82.4 | +40.0 | +0.548 | 5.9 × 10⁻³⁷ |
|  | *STIMATE* | 14.31 | 12.01 | 12.49 | 13.72 | 6.6 | +9.9 | −0.050 | 0.72 |
|  | *CRACR2A* | 8.31 | 4.35 | 7.25 | 4.97 | 3.8 | −31.4 | −0.560 | 4.3 × 10⁻⁴ |
|  | *CRACR2B* | 1.18 | 1.55 | 1.03 | 1.77 | 0.5 | +72.0 | +0.925 | 0.017 |
|  | *CBARP* | 14.57 | 6.56 | 12.73 | 7.49 | 6.7 | −41.1 | −1.254 | 3.1 × 10⁻²⁵ |
|  | Family total | 217.64 | 216.15 | 189.95 | 246.96 | 100.0 | +30.0 | — | — |
| Mitochondrial Ca²⁺ uptake | *MCU* | 26.79 | 8.21 | 23.39 | 9.38 | 13.5 | −59.9 | −1.041 | 2.9 × 10⁻²³ |
|  | *MICU1* | 70.93 | 50.30 | 61.91 | 57.47 | 35.7 | −7.2 | −0.069 | 0.41 |
|  | *MICU2* | 29.50 | 12.21 | 25.72 | 13.95 | 14.8 | −45.8 | −0.566 | 6.9 × 10⁻⁷ |
|  | *MICU3* | 3.20 | 2.97 | 2.79 | 3.40 | 1.6 | +21.6 | +0.565 | 0.028 |
|  | *MCUR1* | 49.95 | 35.18 | 43.58 | 40.17 | 25.1 | −7.8 | +0.044 | 0.62 |
|  | *MCUB* | 18.28 | 4.69 | 15.96 | 5.36 | 9.2 | −66.4 | −1.655 | 1.5 × 10⁻³⁰ |
|  | Family total | 198.64 | 113.57 | 173.34 | 129.73 | 100.0 | −25.2 | — | — |
| PMCA (Ca²⁺ extrusion) | ***ATP2B1*** | 20.26 | 18.77 | 17.70 | 21.45 | 59.0 | +21.2 | +0.290 | 2.8 × 10⁻⁵ |
|  | *ATP2B2* | 2.79 | 6.54 | 2.43 | 7.48 | 8.1 | +207.3 | +1.346 | 8.6 × 10⁻³⁹ |
|  | *ATP2B3* | 0.31 | 0.97 | 0.27 | 1.11 | 0.9 | +306.6 | +2.584 | 2.2 × 10⁻²⁴ |
|  | *ATP2B4* | 10.99 | 9.21 | 9.58 | 10.52 | 32.0 | +9.8 | +0.308 | 5.0 × 10⁻⁶ |
|  | Family total | 34.34 | 35.49 | 29.99 | 40.56 | 100.0 | +35.2 | — | — |

*Note.* TPM (transcripts per million), mean of three libraries per group; adjusted TPM, after per-library median-of-ratios scaling for library composition (Methods 2.10). Share, percentage of the family total in control cells; bold, dominant member (> 50%). log2FC and p_adj from DESeq2 on the gene-level Salmon counts; —, not computed; CRAC, Ca²⁺ release-activated Ca²⁺.

<!-- /table:4 -->

**Table 5.** Cross-disease comparison of *TRPC1*, *SARAF* and *CBARP*.

<!-- table:5 -->

| Cohort | Region | n, case/control | *TRPC1* δ (q) | *SARAF* δ (q) | *CBARP* δ (q) |
|:----------------------|:------------------------|:----------|----------------:|----------------:|----------------:|
| ALS (NYGC, GSE153960) | Cerebellum | 158/38 | **+0.538 (1.5 × 10⁻⁶)** | **+0.599 (2.3 × 10⁻⁷)** | **−0.544 (1.2 × 10⁻⁶)** |
|  | Frontal cortex | 154/56 | **+0.486 (4.3 × 10⁻⁷)** | **+0.525 (6.1 × 10⁻⁸)** | **−0.511 (1.2 × 10⁻⁷)** |
|  | Motor cortex (lateral) | 82/18 | **+0.511 (8.8 × 10⁻³)** | **+0.415 (0.025)** | **−0.509 (8.9 × 10⁻³)** |
|  | Motor cortex (medial) | 81/19 | **+0.465 (0.013)** | **+0.402 (0.027)** | **−0.523 (8.0 × 10⁻³)** |
|  | Occipital cortex | 45/11 | +0.176 (0.67) | +0.277 (0.51) | −0.160 (0.70) |
|  | Temporal cortex | 25/24 | **+0.447 (0.038)** | **+0.503 (0.025)** | **−0.633 (0.019)** |
|  | Hippocampus | 29/11 | **+0.699 (2.1 × 10⁻³)** | **+0.956 (1.9 × 10⁻⁴)** | **−0.636 (5.0 × 10⁻³)** |
|  | Spinal cord (cervical) | 155/41 | +0.058 (0.63) | **+0.479 (1.1 × 10⁻⁵)** | **−0.461 (2.4 × 10⁻⁵)** |
|  | Spinal cord (lumbar) | 140/43 | +0.112 (0.33) | **+0.499 (3.9 × 10⁻⁶)** | **−0.455 (2.5 × 10⁻⁵)** |
|  | Spinal cord (thoracic) | 43/10 | −0.321 (0.49) | +0.102 (0.88) | −0.172 (0.77) |
| Other neurological disorders (NYGC, same controls) | Cerebellum | 49/38 | **−0.407 (1.0 × 10⁻²)** | −0.061 (0.76) | −0.088 (0.64) |
|  | Frontal cortex | 45/56 | **−0.717 (6.5 × 10⁻⁸)** | **−0.437 (5.3 × 10⁻⁴)** | +0.094 (0.49) |
|  | Temporal cortex | 35/24 | **−0.571 (2.5 × 10⁻³)** | **−0.500 (6.6 × 10⁻³)** | −0.283 (0.12) |
| Alzheimer's disease (GSE125583) | Fusiform gyrus | 219/70 | **−0.447 (1.6 × 10⁻⁷)** | — | — |
| Parkinson's disease (GSE68719) | BA9 | 29/44 | **−0.677 (1.8 × 10⁻⁴)** | — | — |
| Multiple sclerosis (GSE123496) | Corpus callosum | 5/5 | −0.520 (0.42) | −0.200 (0.73) | +0.760 (0.16) |
|  | Frontal cortex | 5/5 | −0.040 (1.0) | +0.040 (1.0) | −0.360 (1.0) |
|  | Hippocampus | 5/5 | −0.440 (1.0) | −0.280 (1.0) | +0.040 (1.0) |
|  | Internal capsule | 5/5 | +0.200 (0.78) | +0.120 (0.84) | −0.520 (0.54) |
|  | Parietal cortex | 5/5 | −0.280 (0.89) | −0.040 (1.0) | −0.200 (0.89) |
|  | Five regions pooled | 25/25 | −0.226 (0.49) | −0.075 (0.88) | −0.043 (0.88) |
| Multiple sclerosis (GSE138614) | Normal-appearing white matter | 21/25 | −0.482 (0.10) | −0.166 (0.74) | +0.051 (0.82) |
|  | Lesions | 52/25 | **−0.594 (1.3 × 10⁻⁴)** | −0.138 (0.52) | +0.327 (0.050) |
|  | All samples, averaged per donor | 10/5 | **−0.840 (0.038)** | −0.120 (0.77) | +0.120 (0.77) |

*Note.* δ, Cliff's delta, case minus control, with the Benjamini–Hochberg q value in parentheses; bold, q < 0.05; —, not tested. n counts samples, or donors in the donor-level row. The NYGC comparisons of each region share its non-neurological controls. BA9, Brodmann area 9; NYGC, New York Genome Center.

<!-- /table:5 -->

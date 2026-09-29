"""Apply the September 29 scientific audit and author-confirmed technical qPCR unit.

Preserve citation fields, section numbering and unedited content. Requires v117
as input; deliberately fails if a replacement is not uniquely located.
"""
from pathlib import Path
import copy
import re
from docx import Document
from docx.text.paragraph import Paragraph
from docx.oxml import OxmlElement

ROOT=Path(__file__).resolve().parents[1]
helper=(ROOT/'code/condense_results_v107.py').read_text()
exec(helper[helper.index('def replace_span('):helper.index('doc = Document(MANUSCRIPT)')])
main=ROOT/'manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx'
doc=Document(main); P=list(doc.paragraphs); touched=set()
fields_before=len(doc._element.xpath('.//w:instrText'))

def ch(i,old,new):
    assert old in P[i].text,(i,old[:90])
    replace_span(P[i],old,new); touched.add(i)

def full(i,new):
    assert not P[i]._p.xpath('.//w:instrText'),('citation field in full replacement',i)
    ch(i,P[i].text,new)

ch(5,'TRPC1, STIM1, ORAI1 and ATP2A3 mRNAs increased 1.7- to 3.2-fold,','TRPC1, STIM1, ORAI1 and ATP2A3 mRNA measurements were 1.7- to 3.2-fold higher (four technical RT-qPCR repeats per group; descriptive),')
ch(5,'The single-plate Fura-2 observation requires independent replication;','The laboratory observations require independent biological replication;')
old=P[19].text[P[19].text.index('Refitting with them'):]
ch(19,old,'After adjustment for two surrogate variables, 1,067 genes met both DE thresholds, compared with 1,694 in the original model; 850 were shared and all shared genes retained their direction. The median absolute difference in log2 fold change was 0.20 across the 14,012 genes with adjusted p values in both models. Including two surrogate variables reduced the residual degrees of freedom from four to two, and changed both effect estimates and their uncertainty. STIM1, TRPC1, ORAI3, SARAF and CBARP retained q < 0.05 in both models, but only ORAI3 and CBARP met both the q and fold-change thresholds in both. ORAI1, ATP2A3, ATP2A2 and STIM2 did not retain q < 0.05 after adjustment. The per-gene comparison is provided in source_data/svaseq_sensitivity_SHSY5Y.csv.')
ch(25,'Gene-level q values from this pipeline are the screening-stage values of the stage-wise procedure, and stageR supplies the confirmation stage that asks which individual transcripts carry the gene-level signal.',
   'DRIMSeq gene-level q values describe the screening stage, and stageR supplies transcript-level confirmation. Features with missing DRIMSeq p values are not eligible for confirmation. Separate IsoformSwitchAnalyzeR/DEXSeq isoform-usage results are identified as such and are not stageR confirmation values.')
ch(27,'so that no conclusion rests on a single extraction rule.','to assess sensitivity to the extraction rules.')
ch(27,'and every conclusion drawn from it holds in each.','and concordance of direction and significance is reported separately.')
ch(37,'so the interaction was formed as a difference of within-batch differences, with the four NMD-inhibition conditions as the unit of inference.',
   'so the interaction was summarised as a difference of within-batch differences for each of the four NMD-inhibition conditions. These contrasts share the same baseline libraries and are not independent biological replicates. We report point estimates descriptively, without inferential p or q values.')
ch(37,'The model, the gene filter and the tests are given','The contrast definition and gene filter are given')
ch(44,'After a low-expression filter (counts per million [CPM] > 0.40), comparisons used log2 CPM with per-sample median centring,',
   'Within each region and comparison, genes were retained if counts per million (CPM) exceeded 1 in at least max(10, floor[0.20 × number of samples]) samples. CPM was then recalculated on retained genes and transformed as log2(CPM + 1), followed by per-sample median centring. Comparisons used')
ch(45,'Independent cohorts were analysed identically:','The comparison cohorts were')
ch(45,'GSE138614 is white matter,','Alzheimer’s disease used the NYGC filtering and normalisation procedure. For Parkinson’s disease, the supplied gene-expression matrix was filtered at values > 1 in at least floor(0.20 × number of samples) samples, transformed as log2(value + 1) and median-centred. The multiple-sclerosis expression analyses used CPM > 0.40 in at least max(3, floor[0.20 × number of samples]) samples, followed by log2(CPM + 1) and median centring; the GSE138614 marker-adjustment analysis required CPM > 0.40 in at least 20 samples. GSE138614 is white matter,')
ch(45,'in which myelin markers are preserved,','with accompanying myelin-marker measurements,')
ch(46,'with Benjamini–Hochberg correction within each group and model.',
   'with Benjamini–Hochberg correction across regions within each group and model. A donor-level sensitivity analysis averaged normalised expression within donor and region before fitting the marker models or testing group differences. Its correction families were regions within each group, gene and adjustment model; these differ from Table 5’s within-region gene families (Supplementary Table S18d).')
ch(47,'to separate them from neuronal content,','to assess sensitivity to a neuronal-content marker,')
ch(51,'(n = 4; TARDBP knockdown and the target genes were measured on separate RNA sets, 14 April – 21 May 2026 and 3 June – 10 July 2026).',
   '(four technical measurements of the same biological sample per group). TARDBP and the target genes were measured on separate RNA sets; recorded assay dates were 14 April – 21 May 2026 and 3 June – 10 July 2026, respectively. These dates do not identify independent cultures or RNA isolations.')
ch(52,'Knockdown was verified in the RNA set collected in April and May','Lower TARDBP mRNA was observed in the RNA set assayed in April and May')
ch(54,'the RNA of Section 2.13 came from parallel cultures of the same transduction on day 5, after selection was complete.',
   'the RNA protocol of Section 2.13 used day-5 cultures after selection was complete. The available records do not establish sample-level pairing between the RNA sets and Fura-2 measurements.')
full(59,'Bioinformatic thresholds and correction families are specified above. All laboratory comparisons are descriptive. RT-qPCR comprises four technical measurements of the same biological sample per group; the TARDBP and target-gene assays used separate RNA sets. Fura-2 comprises three wells per group from one culture plate. WST-1 comprises four wells from one of three experiments; the other two experiments are not available in the package. We report means, technical-repeat or well-to-well SEM, and observed differences without inferential p values, biological confidence intervals or standardised effect sizes. These measurements do not estimate between-experiment variability. The Fura-2 readdition-to-release ratio is also summarised per well.')
full(62,'Mean TARDBP mRNA was 94.4% lower in the shTDP-43 sample than in the non-targeting shRNA sample and 94.8% lower than in the non-transduced sample (four technical RT-qPCR measurements per group; descriptive; Figure 1A).')
ch(63,'Relative to the non-targeting shRNA control, mRNA levels of four SOCE-associated targets increased: TRPC1 ≈ 1.8-fold (Holm-adjusted p = 0.00091), STIM1 ≈ 1.9-fold (adjusted p = 0.0011), ORAI1 ≈ 1.7-fold (adjusted p = 0.0011) and ATP2A3/SERCA3 ≈ 3.2-fold (adjusted p = 0.00038; n = 4, two-tailed t-tests; Figure 1B).',
   'Relative to the non-targeting shRNA sample, mean measurements of four SOCE-associated mRNAs were higher: TRPC1 ≈ 1.8-fold, STIM1 ≈ 1.9-fold, ORAI1 ≈ 1.7-fold and ATP2A3/SERCA3 ≈ 3.2-fold (four technical measurements per group; Figure 1B). These descriptive differences require biological replication.')
full(65,'Figure 1. Descriptive laboratory measurements in SH-SY5Y cells. (A) TARDBP mRNA in untransduced, non-targeting shRNA and shTDP-43 samples. (B) Relative TRPC1, STIM1, ORAI1 and ATP2A3 mRNA. Each RT-qPCR group comprises four technical measurements of the same biological sample; panels A and B use separate RNA sets. (C) WST-1 signal at 48 h, four wells from one experiment. Bars show means with technical-repeat or well-to-well SEM; dots show individual measurements. No inferential significance tests are assigned to these comparisons.')
ch(66,'The ratio pattern suggests that a smaller releasable store alone may not explain the within-plate difference, but it cannot establish that this response recurs across independent cultures.',
   'The ratio describes the within-plate pattern but does not isolate store depletion, influx or clearance, and cannot establish replication across independent cultures.')
ch(86,'the stage-wise confirmation did not single out any transcript: all seven STIM1 transcripts carried a confirmation-stage adjusted p of 1.0.',
   'none of seven evaluable transcripts passed stageR confirmation (all confirmation-stage adjusted p values = 1.0). Two additional transcripts in the DRIMSeq output, ENST00000698912.1 and ENST00000698913.1, had missing p values and were not evaluated by stageR (source_data/STIM1_transcript_test_eligibility.csv).')
ch(86,'Two STIM1 isoforms are predicted to carry premature termination codons, and neither is among the ones whose usage changes (isoform-level q = 0.33 and 0.77), whereas the two that do change are not predicted to carry them. The gene-level signal is therefore not attributable to a premature-termination-codon isoform.',
   'In the separate IsoformSwitchAnalyzeR/DEXSeq analysis, the two predicted premature-termination-codon (PTC) isoforms had q = 0.33 and 0.77. The two isoforms with significant usage changes were ENST00000698912.1 and ENST00000698913.1, neither predicted to carry a PTC; their stageR confirmation was unavailable. No confirmed PTC-associated isoform switch was identified.')
ch(94,'so the specificity of the procedure rests on the deeper comparisons, in which the same permissive analysis recovered thirteen of sixteen positive controls in SH-SY5Y and fifteen in iPSC colonies. The background call rate is therefore not specific to TDP-43, whereas the genes the analysis identifies are.',
   'so this matched comparison does not establish specificity. Deeper TDP-43 comparisons recovered thirteen of sixteen positive controls in SH-SY5Y and fifteen in iPSC colonies, supporting sensitivity to known targets in those models. The FUS/TAF15 comparison remains descriptive.')
ch(96,'observed SOCE phenotype','within-plate Fura-2 difference')
ch(98,'At junction level the exon (chr4:27,007,983–27,008,006) was positive in direction in all five TDP-43 comparisons in which it was measurable, but small and inconsistent:',
   'At junction level, both flanks of the exon (chr4:27,007,983–27,008,006 in humans) were considered. Sixteen flanking-junction records were measurable across seven TDP-43 comparisons and the FUS/TAF15 controls; direction was positive in six TDP-43 comparisons and negative in C2C12. The human downstream estimates included')
ch(98,'among all junctions leaving that donor, the more sensitive and the noisier measure for a lowly used exon. On either measure the magnitude is too small to support an isoform switch, and the negative conclusion of the meta-analysis stands.',
   'among junctions sharing the downstream acceptor. In C2C12, the upstream and downstream estimates were −0.026 and −0.019 (q = 0.54 and 0.76); in NSC34 they were +0.021 and +0.031 (q = 0.73 and 0.50). These junction summaries do not alter the independently computed rMATS meta-analysis, which did not detect a shared shift.')
full(101,'The TDP-43 × NMD-inhibition experiment (GSE307054) was examined descriptively because its four contrasts share baseline libraries and TDP-43 status is confounded with sequencing batch (Methods 2.8; Supplementary Results 2; Supplementary Figure S5). CBARP had a mean interaction of +1.52 log₂, positive in all four conditions (range +1.19 to +2.23; Supplementary Table S10). This is a candidate pattern, not evidence that CBARP is an NMD target. The literature cryptic-splicing reference genes were not selected as validated NMD-rescue controls, and gene-level counts can dilute an isoform-specific response. Their descriptive panel summary therefore does not establish assay sensitivity or make a negative SOCE result conclusive (Supplementary Table S10b).')
ch(108,'The two mouse lines gave 74 qualifying units in C2C12 and 131 in NSC34, neither with a positive control for this assay, because the',
   'The mouse datasets yielded 74 qualifying units in C2C12 and 131 in NSC34. Neither had a validated positive control for this coverage assay: the')
ch(108,', so the Stmn2 units cannot serve as one, so neither line can say whether an absent signal means an absent event.',
   '. Negative results therefore have limited interpretability.')
ch(108,'One unit is positive in both motor-neuron models, the SARAF intron 5 index (+0.097 in the iPSC-derived motor neurons and +0.204 in NSC34, against 0.000 in SH-SY5Y); we record it as a candidate rather than a finding, for the reasons set out in Supplementary Results 3, where the mouse values are given in full.',
   'Positive gradients occurred in units labelled intron 5 at the human SARAF locus (+0.097 in iPSC-derived motor neurons) and the mouse Saraf locus (+0.204 in NSC34); the human SH-SY5Y estimate was 0.000. The units were matched by ordinal label, not sequence alignment, so these observations do not demonstrate recurrence of the same RNA event across species. They remain gene-level candidates (Supplementary Results 3).')
ch(109,'the first because batch is confounded with TDP-43 status,','the first because its contrasts share controls and batch is confounded with TDP-43 status,')
full(110,'3.8 Matched analyses do not establish enrichment of splicing changes in Ca²⁺ genes')
full(114,'Two of the 24 dataset–panel combinations had nominal permutation p < 0.05 after matching. Neither survived adjustment for the number of comparisons. The panels are nested, and the number of nominal results does not establish enrichment or identify chance as their cause. In the primary model, matching for gene properties reduced the apparent excess without demonstrating a remaining enrichment (Supplementary Table S4).')
ch(118,'TRPC1 was unchanged','No statistically significant TRPC1 difference was detected') if 'TRPC1 was unchanged' in P[118].text else None
ch(118,'well-powered','larger') if 'well-powered' in P[118].text else None
ch(119,'(Supplementary Table S18).','(Supplementary Table S18). Only the ALS cerebellar expression comparison included repeated donor samples: 158 samples represented 147 donors. After averaging within donors, TRPC1 remained higher (δ = +0.532, q = 2.3 × 10⁻⁶) and remained significant in all four marker models (adjusted δ = +0.248 to +0.370; Supplementary Table S18d).')
ch(121,'Three independent diseases, in four datasets, extended this','TRPC1 was also examined in four independent datasets covering Alzheimer’s disease, Parkinson’s disease and multiple sclerosis')
ch(121,'which was underpowered','which included five donors per group')
ch(121,'Two controls argue against that explanation:','Two sensitivity analyses assessed that explanation:')
ch(121,'where the myelin markers are unchanged (δ = −0.482 at sample level, p = 0.005; −0.771 at donor level, p = 0.030), and it survives regression of TRPC1 on MBP, PLP1 and GFAP, although the adjusted difference is smaller and, at donor level, is not significant.',
   'where no statistically significant myelin-marker changes were detected (TRPC1 δ = −0.482 at sample level, p = 0.005, q = 0.10; −0.771 at donor level, uncorrected p = 0.030). Regression on MBP, PLP1 and GFAP attenuated the TRPC1 difference, which remained significant at sample level but not at donor level (p = 0.055). These analyses do not exclude a contribution from tissue composition.')
ch(128,'and with SNAP25 held constant no association remained','and after adjustment for SNAP25 the associations weakened and did not reach statistical significance')
ch(132,'This pattern suggests that store release alone may not explain the smaller readdition response, but it needs replication in independent cultures.',
   'This ratio is descriptive and does not isolate altered store depletion from influx or clearance; the response requires replication in independent cultures.')
ch(133,'Faster plasma-membrane extrusion or ER re-uptake could reduce the cytosolic peak even if channel influx were unchanged.',
   'Faster plasma-membrane extrusion could reduce the cytosolic peak even if channel influx were unchanged. SERCA was inhibited by CPA during the readdition protocol, so increased ER re-uptake cannot be assumed to explain that peak without evidence of residual SERCA activity or inhibitor washout.')
ch(136,'Because metabolic activity was','Because the WST-1 signal was') if 'Because metabolic activity was' in P[136].text else None
ch(137,'Their recovery, together with their absence from the corresponding FUS and TAF15 knockdowns, supports the biological specificity of the analysis.',
   'Their recovery supports detection of known targets in the deeper TDP-43 comparisons. In the matched motor-neuron experiment, STMN2 and KALRN were recovered only after TDP-43 knockdown, but two controls versus none after FUS or TAF15 knockdown did not establish a significant specificity difference.')
ch(137,'and the NMD interaction analysis yielded no genome-wide significant gene.','and the NMD interaction estimates were descriptive and did not confirm an NMD target.')
ch(137,'whereas the calcium phenotype is accompanied by','whereas the independent calcium-regulatory RNA analyses identify')
ch(139,'the functional decrease','the within-plate Fura-2 difference')
ch(139,'In the iPSC-derived motor neurons two panels survived matching, which is within what 24 nested combinations produce by chance.',
   'In the iPSC-derived motor neurons two panels retained nominal associations after matching, but neither met the multiple-comparison criterion.')
ch(141,'The present findings extend this framework to the plasma-membrane replenishment phase. The working model in Figure 6 is that TDP-43 loss reduces metabolic and ER-refilling capacity while simultaneously changing the composition and feedback control of the entry machinery;',
   'The working model in Figure 6 proposes hypotheses about the plasma-membrane replenishment phase. TDP-43 loss might alter cellular bioenergetics, ER refilling, channel composition and feedback control;')
ch(141,'Correcting ATP2A2 or cellular bioenergetics should improve both ER release and readdition, whereas reducing SARAF or restoring the ORAI balance should preferentially improve the readdition phase.',
   'ATP2A2 or bioenergetic perturbations could be tested with direct ER Ca²⁺ measurements and refilling after CPA washout. SARAF or ORAI perturbations could be tested for effects on the readdition response in independently replicated cultures.')
ch(141,'place the functional phenotype within a directly testable neurochemical framework.','test whether the candidate RNA changes contribute to a reproducible calcium response.')
ch(143,'Working model of the observed calcium phenotype','Working model linking descriptive calcium measurements to candidate hypotheses')
ch(143,'The Fura-2 comparison is descriptive.','RT-qPCR comprises four technical measurements per group; both the qPCR and Fura-2 comparisons are descriptive.')
ch(145,'TARDBP depletion was verified by RT-qPCR but not at the protein level,','Lower TARDBP mRNA was measured by RT-qPCR but depletion was not verified at the protein level,')
ch(145,'the well-level Welch p value and confidence interval calculated previously cannot support a population-level inference.',
   'it does not estimate between-experiment variability.')
ch(145,'Because metabolic activity was 38.5% lower,','Because the WST-1 signal was 38.5% lower,')
ch(145,'RT-qPCR was normalised to a single reference gene, GAPDH, whose Ct did not differ between groups.',
   'RT-qPCR used four technical measurements of the same biological sample per group and therefore does not provide independent biological confirmation of the RNA-seq changes. It was normalised to a single reference gene, GAPDH, whose mean Ct was similar between the measured groups.')
ch(145,'Normalising readdition to the release of the same well removes differences in store content but not their consequences, because a smaller release also means weaker STIM activation, so the ratio bounds rather than eliminates the contribution of the store.',
   'The readdition-to-release ratio is descriptive and cannot distinguish altered ER store depletion from changes in STIM activation, calcium influx or clearance.')
ch(148,'conservative NMD analyses','descriptive NMD summaries')
ch(154,'v1.0.5','v1.0.6') if P[154].text.count('v1.0.5')==1 else None
if 'v1.0.5' in P[154].text:
    for node in P[154]._p.xpath('.//w:t'): node.text=(node.text or '').replace('v1.0.5','v1.0.6')
    for rel in doc.part.rels.values():
        if rel.is_external and 'releases/tag/v1.0.5' in rel.target_ref: rel._target=rel.target_ref.replace('v1.0.5','v1.0.6')
ch(154,'including per-replicate values and the statistics behind every laboratory figure panel,','including technical-repeat values and the descriptive summaries behind every laboratory figure panel,')
full(156,'During preparation of this manuscript, the authors used Anthropic Claude and OpenAI Codex to assist with drafting, language revision, analysis-script development and revision, figure preparation, formatting and consistency checks. The authors reviewed and edited the output and take full responsibility for the content of the submitted article.')

suppath=ROOT/'supplementary/SUPPLEMENTARY_MATERIAL.docx'
sup=Document(suppath)
def sp(prefix,new):
    p=next((p for p in sup.paragraphs if p.text.startswith(prefix)),None)
    if p is None and prefix=='S10b.':
        anchor=next(p for p in sup.paragraphs if p.text.startswith('S10.'))
        node=OxmlElement('w:p'); anchor._p.addnext(node)
        p=Paragraph(node,anchor._parent); p.style=anchor.style; p.add_run(new); return
    assert p is not None,prefix
    assert not p._p.xpath('.//w:instrText')
    replace_span(p,p.text,new)
sp('The analyses below were carried out', 'Supplementary analyses provide the methods and descriptive results underlying the main-text summaries.')
sp('for each of the four NMD-inhibition conditions',
   'for each of four NMD-inhibition conditions c, after averaging the two replicates within each condition on the log2(normalised count + 1) scale. Genes were retained if mean normalised counts were at least 10 across the TDP-43-plus-NMD-inhibition samples and at least 5 in the TDP-43-only samples (19,145 genes). A positive interaction is compatible with preferential stabilisation after TDP-43 loss, but does not identify a particular transcript or establish NMD-mediated degradation.')
sp('Because the four interventions reuse',
   'The four interventions reuse the same control and TDP-43-knockdown baseline libraries. Averaging replicates does not remove the resulting covariance, and variation across the four contrasts does not include uncertainty in the shared baseline difference. We therefore report the mean, range and number of positive condition-level interactions descriptively, without t tests, sign tests, panel-enrichment tests or FDR claims. The full gene-level descriptive table is supplied as source_data/nmd_descriptive_all_genes.csv.gz.')
sp('CBARP',
   'CBARP had the largest mean interaction among the genes in Supplementary Table S10 (+1.52 log₂), with all four conditions positive (range +1.19 to +2.23). This pattern may motivate transcript-specific follow-up, but the shared controls, batch structure and gene-level measurement prevent it from establishing CBARP as an NMD target.')
sp('Two negative results follow.',
   'STIM1 had a mean interaction of −0.325 log₂, with one of four conditions positive (range −0.850 to +0.333). The separate SH-SY5Y isoform analysis identified no confirmed PTC-associated switch; it cannot explain or exclude a response in this i3Neuron experiment. The four calcium-panel medians and the 16 cryptic-splicing reference genes are summarised descriptively in Supplementary Table S10b. These reference genes were selected for cryptic splicing, not validated NMD rescue in this experiment. Their median interaction of −0.11 log₂ therefore does not demonstrate failure of assay sensitivity. Gene-level counts may also obscure an isoform-specific NMD response, so an absent gene-level interaction is not evidence that a SOCE transcript escapes NMD.')
sp('No SOCE-machinery unit exceeded',
   'No SOCE-machinery unit exceeded |Δ| = 0.30 in either mouse line. The largest estimates, Atp2a2 intron 6 (+0.285) and Trpc1 intron 7 (+0.260) in C2C12, had intervals including zero and were not reproduced in NSC34. Units labelled intron 5 had positive gradients at the human SARAF locus in iPSC-derived motor neurons (+0.097) and the mouse Saraf locus in NSC34 (+0.204); the human SH-SY5Y estimate was 0.000 (interval −0.264 to +0.241). Because human and mouse units were matched by ordinal number rather than sequence alignment, these are gene-level candidates, not confirmed recurrence of the same event. The NSC34 interval is wide and a coverage gradient does not localise a poly(A) site (Supplementary Table S11).')
for p in sup.paragraphs:
    if 'were unchanged' in p.text: replace_span(p,'were unchanged','showed no statistically significant differences')
    if 'the study was underpowered' in p.text: replace_span(p,'the study was underpowered','the study included five donors per group')
sp('S10.', 'S10. Descriptive TDP-43 × NMD-inhibition interactions for calcium-panel and reference genes: four condition estimates, mean, range and positive-condition count. No inferential p or q values are reported.')
sp('S10b.', 'S10b. Descriptive calcium-panel and cryptic-splicing-reference summaries of the NMD interactions. The cryptic-splicing references are not a validated NMD-positive control panel.')
sp('Supplementary Figure S5.', 'Supplementary Figure S5. Descriptive TDP-43 × NMD-inhibition interactions for genes with a cryptic junction and other genes, with CBARP marked. Each gene’s value is the mean of four contrasts sharing the same baseline libraries. No inferential significance is assigned; TDP-43 status is confounded with sequencing batch.')
entry='S18d. Donor-level NYGC sensitivity analysis: normalised expression averaged within donor and region before testing TRPC1, SARAF and CBARP and before fitting the TRPC1 marker-adjustment models. Correction families comprise regions within each group, gene and model.'
anchor=next(p for p in sup.paragraphs if p.text.startswith('S18c.'))
node=OxmlElement('w:p'); anchor._p.addnext(node); new=Paragraph(node,anchor._parent); new.style=anchor.style; new.add_run(entry)

# Supplementary inventories in the main text use the same descriptions.
for prefix in ('S10.','S10b.'):
    new=next(p.text for p in sup.paragraphs if p.text.startswith(prefix))
    p=next((p for p in doc.paragraphs if p.text.startswith(prefix)),None)
    if p is None:
        anchor=next(p for p in doc.paragraphs if p.text.startswith('S10.'))
        node=OxmlElement('w:p'); anchor._p.addnext(node)
        p=Paragraph(node,anchor._parent); p.style=anchor.style; p.add_run(new)
    else: replace_span(p,p.text,new)
anchor=next(p for p in doc.paragraphs if p.text.startswith('S17.'))
for prefix in ('S18.','S18b.','S18c.','S18d.'):
    txt=next(p.text for p in sup.paragraphs if p.text.startswith(prefix))
    node=OxmlElement('w:p'); anchor._p.addnext(node); new=Paragraph(node,anchor._parent)
    new.style=anchor.style; new.add_run(txt); anchor=new

# Remove accidental whole-clause italics and restore gene-only italics in revised
# body text. Nontext XML, citation fields and heading emphasis are preserved.
genes=set(re.findall(r'"([A-Za-z0-9]+)"',helper[helper.index('GENES ='):helper.index('# paragraph index')]))|{'RBFOX3','TBP','B2M','GAPDH','Atp2a2','Trpc1'}
pat=re.compile(r'\b('+'|'.join(sorted(genes,key=len,reverse=True))+r')\b')
def genes_only(p):
    for run in list(p.runs):
        if not run.text or run._r.xpath('./w:fldChar|./w:instrText|./w:drawing'): continue
        text=run.text
        # Only ordinary text runs can be split; references in field instructions remain intact.
        if len(run._r.xpath('./w:t'))!=1: continue
        parent=run._r.getparent(); at=parent.index(run._r)
        parts=re.split(pat,text)
        for part in parts:
            if not part: continue
            clone=copy.deepcopy(run._r)
            rr=__import__('docx.text.run',fromlist=['Run']).Run(clone,p)
            rr.text=part; rr.italic=part in genes
            parent.insert(at,clone); at+=1
        parent.remove(run._r)
for i in touched|{97}:
    if not P[i].style.name.startswith('Heading'): genes_only(P[i])
for p in sup.paragraphs:
    if not p.style.name.startswith('Heading') and not p.text.startswith(('TDP-43 knockdown','Elmasnur')): genes_only(p)
for d in (doc,sup):
    for p in d.paragraphs:
        if p.text.startswith('S11.'):
            for r in p.runs: r.bold=False
assert len(doc._element.xpath('.//w:instrText'))==fields_before,'citation fields changed'
doc.save(main); sup.save(suppath)

hp=ROOT/'highlights_Neurochemistry_International.docx'; h=Document(hp)
for p in h.paragraphs:
    if 'RT-qPCR' in p.text or 'mRNAs increased' in p.text:
        replace_span(p,p.text,'Technical RT-qPCR measurements showed higher calcium-regulatory mRNAs.')
h.save(hp)
print('Manuscript, supplementary text and highlights revised; citation fields preserved:',fields_before)

"""Apply the 25 September 2026 revision (third referee report) to the live DOCX.

Run once, against the v1.0.4 manuscript. Text edits operate on Word text nodes, so the
embedded Zotero citation fields are untouched. Figures:
  Figure 2  panel E added (code/fig2_common_scale_panel.py)
  Figure 4  new: splicing evidence in SOCE-related genes (code/fig_splicing_revision.py)
  Figure 5  the ALS tissue figure, renumbered from Figure 4
  Figure 6  redrawn so that arrows separate measurements from hypotheses
The former Figure 5 (CBARP events) is replaced; the CBARP locus is Supplementary Figure S3.
"""
from pathlib import Path

from docx import Document
from docx.shared import Emu, Inches

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx"
MAIN = ROOT / "figures/main"


def replace(paragraph, old, new):
    nodes = paragraph._p.xpath(".//w:t")
    parts = [node.text or "" for node in nodes]
    whole = "".join(parts)
    if whole.count(old) != 1:
        raise ValueError(f"Expected exactly one occurrence of {old!r}; got {whole.count(old)}")
    start = whole.index(old)
    end = start + len(old)
    offsets, position = [], 0
    for part in parts:
        offsets.append((position, position + len(part)))
        position += len(part)
    first = next(i for i, (a, b) in enumerate(offsets) if a <= start < b)
    last = next(i for i, (a, b) in enumerate(offsets) if a < end <= b)
    a, _ = offsets[first]
    z, _ = offsets[last]
    if first == last:
        nodes[first].text = parts[first][: start - a] + new + parts[first][end - a:]
    else:
        nodes[first].text = parts[first][: start - a] + new
        for i in range(first + 1, last):
            nodes[i].text = ""
        nodes[last].text = parts[last][end - z:]


def swap_image(paragraph, png, width_in):
    """Replace the picture in a paragraph, keeping its position, and rescale its height."""
    from PIL import Image
    blip = paragraph._p.xpath(".//a:blip")[0]
    rid = blip.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed")
    part = paragraph.part.related_parts[rid]
    part._blob = png.read_bytes()
    w, h = Image.open(png).size
    cx = int(Inches(width_in)); cy = int(cx * h / w)
    for ext in paragraph._p.xpath(".//wp:extent") + paragraph._p.xpath(".//a:ext"):
        ext.set("cx", str(cx)); ext.set("cy", str(cy))


def add_figure_after(doc, anchor, png, number, caption, width_in=6.15):
    doc.add_picture(str(png), width=Inches(width_in))
    picture = doc.paragraphs[-1]
    picture.paragraph_format.keep_with_next = True
    cap = doc.add_paragraph(style="Caption")
    cap.add_run(f"Figure {number}. ").bold = True
    cap.add_run(caption)
    anchor._p.addnext(cap._p)
    anchor._p.addnext(picture._p)


doc = Document(MANUSCRIPT)
p = list(doc.paragraphs)
assert p[92].text.startswith("Figure 5. CBARP splicing") and p[128].text.startswith("Figure 4. ALS tissue")

# ---------------------------------------------------------------- Abstract
replace(p[5], "while the Fura-2 Ca²⁺-readdition amplitude was approximately 84% lower after knockdown (n = 3). ER Ca²⁺ release amplitudes were 0.268 and 0.180 Δ(F340/F380), respectively.",
        "while the Fura-2 Ca²⁺-readdition amplitude was 84% lower after knockdown (n = 3 independent cultures; Welch’s t-test p = 0.035). ER Ca²⁺ release was 33% lower (p = 0.31), and readdition normalised to each culture’s own release remained 77% lower.")

# ---------------------------------------------------------------- Methods
replace(p[29], "(Supplementary Table S14). We therefore adopted an effect-based definition:",
        "(Supplementary Table S14). Because this definition produced about as many split-control calls as real calls, or more (Section 3.5), we adopted an effect-based definition:")
replace(p[46], "The public label does not identify the constituent disorders; these regions",
        "The public metadata do not identify the constituent disorders, which the NYGC provides on request; these regions")
replace(p[46], "Age, RNA integrity number and post-mortem interval were not included as covariates.",
        "Age, RNA integrity number and post-mortem interval were not included as covariates. Per-sample median centring removes sample-wide expression shifts and the shared controls equalise control-side confounders, but neither adjusts for case–control differences in these variables, so residual confounding by them cannot be excluded.")
replace(p[54], "(n = 4 each).",
        "(n = 4 each). Reporting guidelines for quantitative PCR recommend at least two validated reference genes; the single-reference design is therefore listed among the limitations.")
replace(p[56], "relative to the respective preceding baseline (n = 3).",
        "relative to the respective preceding baseline. Each group comprised three independent cultures (n = 3), and both phases were read from the same recording of each culture.")
replace(p[59], "At 48 h, n = 4; values are reported as mean ± SEM.",
        "At 48 h, n = 4 wells per group from one of the three experiments performed; values are reported as mean ± SEM.")
replace(p[61], "Fura-2 amplitudes (n = 3) and 48-h WST-1 well signals (n = 4) were summarised descriptively, without inferential tests. Significance was set at adjusted p < 0.05 for the four-target RT-qPCR family.",
        "Fura-2 amplitudes from three independent cultures per group were compared with Welch’s two-tailed t-test (SciPy v1.12.0), which does not assume equal variances, and are reported with the difference in means and its 95% confidence interval; these two pre-specified comparisons were not adjusted for multiplicity. The readdition-to-release ratio of each culture was summarised descriptively. The 48-h WST-1 wells (n = 4) come from one experiment and were summarised descriptively, because wells within an experiment are not independent replicates. Significance was set at p < 0.05, Holm-adjusted for the four-target RT-qPCR family.")

# ---------------------------------------------------------------- Results 3.1
replace(p[68], "The ER Ca²⁺-release means were 0.268 ± 0.042 and 0.180 ± 0.061 Δ(F340/F380) in the control and shTDP-43 groups, respectively. The readdition amplitudes were 1.542 ± 0.282 and 0.245 ± 0.083, respectively (n = 3). At 48 h, the WST-1 signal was 61.5 ± 0.8% of control, a 38.5% decrease (n = 4). These Fura-2 and WST-1 comparisons are descriptive.",
        "The readdition amplitude fell from 1.542 ± 0.282 to 0.245 ± 0.083 Δ(F340/F380) (mean ± SEM, n = 3 independent cultures per group; difference −1.30, 95% CI −2.40 to −0.19; Welch’s t-test p = 0.035). ER Ca²⁺ release fell less and not significantly, from 0.268 ± 0.042 to 0.180 ± 0.061 (difference −0.09, 95% CI −0.30 to +0.13; p = 0.31). Because both phases come from the same recording, readdition was also expressed relative to each culture’s own release. This ratio was 5.88 ± 1.10 in controls and 1.36 ± 0.01 after knockdown, 77% lower, and every shTDP-43 culture had a lower ratio than every control culture (Supplementary Table S1). A smaller releasable store therefore does not account for most of the decrease. Plotted on common axes, the two representative recordings had similar baselines before readdition (Figure 2E). At 48 h, the WST-1 signal was 61.5 ± 0.8% of control, a 38.5% decrease (n = 4 wells; descriptive).")
replace(p[71], "Panels C and D show individual measurements and mean ± SEM (n = 3; descriptive).",
        "Panels C and D show the three independent cultures per group and mean ± SEM (readdition, Welch’s t-test p = 0.035; ER release, p = 0.31). (E) The recordings in A and B replotted from their exported ratio values on common axes, with time aligned to the steepest point of the readdition rise.")
swap_image(p[70], MAIN / "Figure2_calcium_responses.png", 6.05)

# ---------------------------------------------------------------- Results 3.4
replace(p[87], "(Table 3, Supplementary Figure S3). A STIM1 event met the count and FDR criteria, but its interval included zero and is shown for context.",
        "(Table 3; Figure 4A).")
replace(p[89], ": 32 events met both the threshold and the coverage criteria in five of the six datasets (iPSC colonies, mouse striatum, SH-SY5Y, C2C12 and NSC34) at the orthologous human chr19 and mouse chr10 loci, with |ΔPSI| from 0.11 to 0.74 (median 0.37) and 61 to 4,559 junction reads per event",
        ". Thirty-two events met both the threshold and the coverage criteria in five of the six datasets (iPSC colonies, mouse striatum, SH-SY5Y, C2C12 and NSC34) at the orthologous human chr19 and mouse chr10 loci. Their |ΔPSI| ranged from 0.11 to 0.74 (median 0.37), with 61 to 4,559 junction reads per event")
replace(p[89], "The event-level evidence is shown in Figure 5. The locus is therefore affected in every model except the iPSC-derived motor neurons, but the direction is not consistent across models: all three qualifying SH-SY5Y events show reduced inclusion of the rMATS-defined form, whereas the events in mouse striatum and C2C12, all but one of those in NSC34 and most of those in iPSC colonies are positive; the junction-level analysis of Section 3.5 quantifies a different arm of the same locus and is positive in SH-SY5Y. What reproduces is the involvement of the locus, not one directional switch.",
        "The locus is therefore affected in every model except the iPSC-derived motor neurons. The sign of the rMATS events varies: all three qualifying SH-SY5Y events show reduced inclusion of their rMATS-defined form, whereas most events in the other four datasets are positive. The sign follows how each event is defined, and none of the three SH-SY5Y events contrasts the canonical exon 4–exon 5 junction with its alternative. Coverage and split reads from the alignments show the change directly (Supplementary Figure S3). In both human models the canonical junction nearly disappeared after knockdown, and exon 4 was instead joined to an annotated alternative 3′ splice site 197 nucleotides upstream of the exon 5 acceptor. This junction carried 17% of exon-4 donor reads in SH-SY5Y controls and 86% after knockdown, and 1% and 78% in iPSC colonies (Figure 4C). The mouse loci were not compared at junction level.")

# replace the former Figure 5 (CBARP events) with the new Figure 4
anchor = p[90]
for par in (p[91], p[92]):
    par._element.getparent().remove(par._element)
add_figure_after(doc, anchor, MAIN / "Figure4_SOCE_splicing.png", 4,
                 "Splicing evidence in SOCE-related genes. (A) Per-library PSI for the four SOCE-machinery skipped-exon events of Table 3 in the public SH-SY5Y comparison (three libraries per group); lines show group means, and the text gives the rMATS ΔPSI with its replicate-level bootstrap 95% interval. The STIM1 interval includes zero. (B) The 24-nucleotide SOAR exon that converts STIM2 into STIM2.1, measured as an rMATS event in six datasets. Squares are sized by inverse-variance weight, bars are 95% intervals and diamonds are fixed-effect pooled estimates for all six datasets and for the three human datasets (Supplementary Table S15). (C) CBARP: share of exon-4 donor reads joined to the alternative 3′ splice site at chr19:1,235,342 rather than to exon 5, per library, in SH-SY5Y and iPSC colonies. Reads were counted from the alignments with the rules of Methods 2.5; the full locus is shown in Supplementary Figure S3.")

# ---------------------------------------------------------------- Results 3.5 and 3.9
replace(p[99], "CBARP showed the largest local change of any SOCE gene in SH-SY5Y (ΔPSI +0.74",
        "CBARP showed the largest local change of any SOCE gene in SH-SY5Y, in the junction from exon 4 to the alternative 3′ splice site (ΔPSI +0.74")
replace(p[99], "one arm of the CBARP event uses a splice site absent from the annotation.",
        "one arm of the CBARP event uses a splice site absent from the annotation (Supplementary Figure S3).")
replace(p[100], "p = 0.914; Supplementary Table S15)", "p = 0.914; Figure 4B; Supplementary Table S15)")
replace(p[100], "but the effect was small and inconsistent: ΔPSI +0.149 in SH-SY5Y at 75 ng/mL (q = 0.033, 20 versus 3 reads), +0.116 at 25 ng/mL (q = 0.14), and +0.020 in iPSC colonies",
        "but the effect was small and inconsistent. ΔPSI was +0.149 in SH-SY5Y at 75 ng/mL (q = 0.033, 20 versus 3 reads), +0.116 at 25 ng/mL (q = 0.14) and +0.020 in iPSC colonies")
replace(p[100], "because the two are not the same quantity: rMATS combines both flanking junctions and normalises by effective length, whereas the junction-level PSI is the share of the single downstream junction among all junctions leaving that donor, which is the more sensitive and the noisier of the two for a lowly used exon.",
        "because the two are not the same quantity. rMATS combines both flanking junctions and normalises by effective length. The junction-level PSI is the share of the single downstream junction among all junctions leaving that donor, which makes it the more sensitive and the noisier measure for a lowly used exon.")
replace(p[121], "(Figure 4; Table 5)", "(Figure 5; Table 5)")
replace(p[125], "(q = 0.10); because several samples come from the same donor,",
        "(q = 0.10). Because several samples come from the same donor,")
replace(p[128], "Figure 4.", "Figure 5.")

# ---------------------------------------------------------------- Discussion and limitations
replace(p[133], "The ER-release mean decreased from 0.268 to 0.180 Δ(F340/F380), while the subsequent SOCE amplitude showed a much larger proportional decrease.",
        "ER release fell by 33% and not significantly, whereas readdition fell by 84% and remained 77% lower when each culture was normalised to its own release, so a smaller releasable store does not account for most of the SOCE decrease.")
replace(p[138], "UNC13A programs described", "UNC13A programmes described")
replace(p[138], "universal cryptic-splicing program.", "universal cryptic-splicing programme.")
replace(p[139], "(Figure 5)", "(Figure 4C; Supplementary Figure S3)")
replace(p[144], "The laboratory Fura-2 assay found a smaller readdition amplitude after TARDBP knockdown in SH-SY5Y cells. Public knockdown RNA-seq identified changes in calcium-regulatory transcripts and recurrent CBARP splicing. ATP2A2, ORAI-family composition and SARAF are potential contributors to store filling, entry and feedback; CBARP is placed on a parallel voltage-gated channel branch. Arrows indicate hypotheses for future perturbation, not demonstrated causal links.",
        "Blue boxes and solid arrows summarise what was measured in the laboratory cultures after TARDBP knockdown: higher TRPC1, STIM1, ORAI1 and ATP2A3 mRNA, a smaller Ca²⁺-readdition amplitude and a smaller, non-significant decrease in ER Ca²⁺ release. Orange boxes are candidate mechanisms suggested by the public RNA-seq data; dashed arrows mark hypotheses for future perturbation, not demonstrated causal links. CBARP is placed on a parallel voltage-gated channel branch that was not measured here.")
swap_image(p[143], MAIN / "Figure6_working_model.png", 6.15)
replace(p[146], "The Fura-2 (n = 3) and 48-h WST-1 (n = 4) results are descriptive.",
        "The Fura-2 comparison rests on three independent cultures per group, and the readdition signal was not characterised pharmacologically, for example with an ORAI channel inhibitor. The WST-1 signal is reported at 48 h only, from four wells of one experiment; it is not a cell count and cannot separate fewer cells from lower metabolic activity per cell. RT-qPCR was normalised to a single reference gene, GAPDH, whose Ct did not differ between groups. In the public inducible RNA-seq comparison GAPDH rose (log2FC +0.61) whereas TBP and B2M were stable; a similar rise in the laboratory cells would have led the target increases to be underestimated, and two-gene normalisation would be preferable.")

doc.save(MANUSCRIPT)
print(MANUSCRIPT)

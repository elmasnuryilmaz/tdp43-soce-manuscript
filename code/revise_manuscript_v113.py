"""Place the comparison-group finding where the manuscript already asks the question.

The paper asks, in two places, whether the tissue TRPC1 differences reflect TDP-43 loss or cell
composition: the composition caveat of the comparison-group paragraph and the cryptic STMN2 proxy
block. The new evidence goes into those two places, and the conclusion is stated once in the
Discussion. Nothing is added as a separate section.

  Methods 2.11   the sample count is 1,640 (the metadata listed one sample twice); groups are
                 single labels; the cell-composition regression; the cryptic indicator in the
                 comparison group
  Results 3.9    a composition paragraph after the ALS results; the comparison group's marker
                 changes replace a sentence whose reasoning the data contradict; two paragraphs
                 after the within-ALS proxy analysis
  Discussion     the tissue paragraph states the conclusion
  Limitations    the comparison group's diagnoses and the partial nature of marker adjustment

Every number comes from Supplementary Tables S18, S18b and S18c
(code/nygc_composition_and_cryptic_s18.py).
"""
import copy
import re
from pathlib import Path

from docx import Document
from docx.oxml.ns import qn
from docx.text.paragraph import Paragraph

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx"
src = (ROOT / "code/condense_results_v107.py").read_text()
exec(src[src.index("def replace_span("):src.index("doc = Document(MANUSCRIPT)")])

GENES = ["TRPC1", "SNAP25", "GFAP", "RBFOX3", "STMN2", "SARAF", "CBARP", "AIF1", "MBP", "PLP1"]
GENE_RE = re.compile(r"\b(" + "|".join(GENES) + r")\b")


def italicise_runs(paragraph):
    """Split non-italic runs so that each gene symbol is italic; field codes are left alone."""
    for r in list(paragraph.runs):
        ts = r._r.findall(qn("w:t"))
        if len(ts) != 1 or r._r.find(qn("w:tab")) is not None or r._r.find(qn("w:br")) is not None:
            continue
        t = ts[0]
        if not t.text or r.italic or not GENE_RE.search(t.text):
            continue
        parts, pos = [], 0
        for m in GENE_RE.finditer(t.text):
            if m.start() > pos:
                parts.append((t.text[pos:m.start()], False))
            parts.append((m.group(0), True)); pos = m.end()
        if pos < len(t.text):
            parts.append((t.text[pos:], False))
        for chunk, it in parts:
            nr = copy.deepcopy(r._r)
            for el in nr.findall(qn("w:t")):
                nr.remove(el)
            el = nr.makeelement(qn("w:t"), {}); el.text = chunk
            el.set("{http://www.w3.org/XML/1998/namespace}space", "preserve"); nr.append(el)
            if it:
                rpr = nr.find(qn("w:rPr"))
                if rpr is None:
                    rpr = nr.makeelement(qn("w:rPr"), {}); nr.insert(0, rpr)
                if rpr.find(qn("w:i")) is None:
                    rpr.append(rpr.makeelement(qn("w:i"), {}))
            r._r.addprevious(nr)
        r._r.getparent().remove(r._r)


def insert_after(paragraph, text):
    """A new paragraph with the template's paragraph properties and plain body text."""
    new = copy.deepcopy(paragraph._p)
    for child in list(new):
        if child.tag != qn("w:pPr"):
            new.remove(child)
    paragraph._p.addnext(new)
    par = Paragraph(new, paragraph._parent)
    par.add_run(text)
    return par


doc = Document(MANUSCRIPT)
p = list(doc.paragraphs)
assert p[44].text.startswith("ALS post-mortem RNA-seq was taken")
assert p[118].text.startswith("Notably, the ALS increase in TRPC1 was confined to brain")
assert p[119].text.startswith("The “Other Neurological Disorders” group within the same cohort")
assert p[126].text.startswith("Within ALS samples the proxy did not track")
assert p[137].text.startswith("Relevance to disease tissue.")
assert p[143].text.startswith("Bulk post-mortem expression can reflect cell composition")

# ---------------------------------------------------------------- Methods 2.11
# the count sits in the suffix of a Zotero citation: change the field code and its result together
n = 0
for tag in ("w:t", "w:instrText"):
    for el in p[44]._p.iter(qn(tag)):
        if el.text and "1,641 samples with metadata after filtering" in el.text:
            el.text = el.text.replace("1,641 samples with metadata after filtering",
                                      "1,640 samples with metadata after filtering")
            n += 1
assert n >= 2, n
replace_span(p[44], "case numbers allowed this comparison in three regions.",
             "case numbers allowed this comparison in three regions. Groups were used as single "
             "labels: 266 samples annotated with both ALS Spectrum MND and Other Neurological "
             "Disorders, and two annotated with Pre-fALS and Other Neurological Disorders, were "
             "assigned to neither group.")
replace_span(p[46], "were carried through all comparisons.",
             "were carried through all comparisons. To ask whether the TRPC1 differences reflected "
             "cell composition, TRPC1 was regressed within each region on a neuronal marker (SNAP25 "
             "or RBFOX3), alone or together with GFAP, by ordinary least squares with an intercept "
             "across all samples of the comparison, as for multiple sclerosis; the residuals were "
             "compared with the Mann–Whitney U test and Cliff’s δ, with Benjamini–Hochberg "
             "correction within each group and model.")
replace_span(p[47], "Libraries were linked to the cohort metadata through the CGND identifier.",
             "Libraries were linked to the cohort metadata through the CGND identifier. The same "
             "indicator was computed for the comparison group, on the samples of the expression "
             "comparison that reached the read threshold; where a sample had been sequenced as more "
             "than one library, the reads were pooled so that each sample counts once. Within that "
             "group, associations with cryptic PSI were tested by Spearman correlation and, to "
             "separate them from neuronal content, by partial Spearman correlation on SNAP25.")

# ---------------------------------------------------------------- Results 3.9
replace_span(p[119], "GFAP and SNAP25 moved more strongly in the comparison group than in ALS while "
                     "TRPC1 fell; these markers make a simple neurodegeneration explanation less "
                     "compelling but do not provide comprehensive cell-composition adjustment.",
             "In this group SNAP25 fell and GFAP rose in frontal and temporal cortex (SNAP25 "
             "δ = −0.57 and −0.55; GFAP δ = +0.59 and +0.64), consistent with neuronal loss and "
             "astrogliosis. Adjustment for SNAP25 removed most of the temporal-cortex decrease in "
             "TRPC1 (δ = −0.571 before and −0.169 after; −0.390 after adjustment for RBFOX3) but "
             "not the frontal or cerebellar decrease (−0.42 to −0.68 and −0.39 to −0.42 across the "
             "four models; Supplementary Table S18).")
composition = insert_after(p[118],
    "Because TRPC1 followed neuronal content in every group and region (Spearman ρ with SNAP25 = "
    "0.34 to 0.86; Supplementary Table S18c), we regressed it on cell-composition markers "
    "(Methods 2.11). The ALS increase was attenuated but not reversed: in the six regions in which "
    "it was significant, the adjusted δ remained positive under every marker combination (+0.26 to "
    "+0.66), and it remained significant in cerebellum, frontal cortex and medial motor cortex "
    "whichever markers were used, and in all six regions after adjustment for RBFOX3 "
    "(Supplementary Table S18).")
group_paragraph = insert_after(p[126],
    "The same indicator characterises the comparison group, whose constituent disorders are not "
    "public. In this group the cryptic junction was common in cortex: it exceeded 1% of the reads at "
    "the exon-1 donor in 21 of 42 frontal and 22 of 35 temporal cortex samples, against 2 of 154 and "
    "2 of 25 ALS samples and none of the controls of the same regions (Cliff’s δ versus controls = "
    "+0.73 and +0.71; q < 10⁻⁶), whereas in cerebellum no comparison-group or ALS sample exceeded 1% "
    "(Supplementary Table S18b). This is consistent with the TDP-43 proteinopathy of "
    "frontotemporal lobar degeneration or limbic-predominant age-related TDP-43 encephalopathy, "
    "although the diagnoses cannot be checked. TRPC1 was nevertheless lower in these regions. "
    "Within the group, TRPC1 fell as cryptic inclusion rose (ρ = −0.50 in frontal and −0.40 in "
    "temporal cortex), but SNAP25 fell with it (ρ = −0.49 and −0.47), and with SNAP25 held constant "
    "no association remained (partial ρ = −0.24, p = 0.12, and 0.00, p = 0.98; Supplementary "
    "Table S18c).")
synthesis = insert_after(group_paragraph,
    "Across the cohort, TRPC1 therefore did not behave as a readout of TDP-43 loss of function. It "
    "rose in ALS brain regions, in which the cryptic junction exceeded 1% of reads in at most 14% of "
    "samples, was unchanged in ALS spinal cord, where it did so in 41–67%, and fell in the cortex of "
    "the comparison group, where it did so in 50–63%; within groups it followed neuronal content "
    "rather than the junction.")

# ---------------------------------------------------------------- Discussion and Limitations
replace_span(p[137], "Because bulk-tissue expression is sensitive to region and cellular composition, "
                     "these observations identify disease-associated candidates rather than a direct "
                     "TDP-43-driven mechanism.",
             "The increase was attenuated but not removed by adjustment for neuronal and astrocytic "
             "markers, so it does not appear to be only a difference in cell composition; it did "
             "not, however, follow the cryptic STMN2 indicator of TDP-43 loss of function, and a "
             "comparison group that carried that indicator in cortex more often than the ALS "
             "samples showed the opposite change. The tissue increase is therefore a "
             "disease-associated candidate whose cause is not established, and the direction it "
             "shares with the cellular model does not by itself link the two.")
replace_span(p[143], "Bulk post-mortem expression can reflect cell composition, and age, RNA integrity "
                     "and post-mortem interval were not modelled as covariates.",
             "Bulk post-mortem expression can reflect cell composition, which regression on marker "
             "genes adjusts for only partially, and age, RNA integrity and post-mortem interval were "
             "not modelled as covariates. The disorders that make up the NYGC comparison group are "
             "not public, and the cryptic STMN2 junction indicates TDP-43 loss of function rather "
             "than a diagnosis.")

for par in (p[44], p[46], p[47], p[119], composition, group_paragraph, synthesis, p[137], p[143]):
    italicise_runs(par)

doc.save(MANUSCRIPT)
print("revised:", MANUSCRIPT.name)

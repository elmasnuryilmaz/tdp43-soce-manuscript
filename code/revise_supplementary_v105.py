"""Replace Supplementary Figure S3 (now the CBARP locus) in SUPPLEMENTARY_MATERIAL.docx.

The former S3 content (replicate-level PSI of STIMATE, ORAI3, STIM2 and STIM1) moved to
main Figure 4A. Run once, against the v1.0.4 supplementary document.
"""
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "supplementary/SUPPLEMENTARY_MATERIAL.docx"
src = (ROOT / "code/revise_manuscript_v105.py").read_text()
helpers = src[src.index("def replace("):src.index("def add_figure_after(")]
exec("from docx.shared import Inches\n" + helpers)

doc = Document(DOCX)
p = doc.paragraphs
assert p[37].text == "Supplementary Figure S3" and p[39].text.startswith("Supplementary Figure S3. Replicate-level PSI")
swap_image(p[38], ROOT / "figures/supplementary/Supplementary_Figure_S3_CBARP_locus.png", 6.1)
old = p[39].text[len("Supplementary Figure S3. "):]
replace(p[39], old,
        "CBARP exon 4–exon 5 region after TDP-43 depletion. (A) Sashimi plot of the public SH-SY5Y "
        "comparison (0 versus 75 ng/mL doxycycline, three libraries per group) and of iPSC colonies "
        "(four libraries per group), drawn from the alignments; the minus-strand gene is shown 5′→3′. "
        "Filled profiles are the mean per-base read depth of the libraries in each group, and each track "
        "has its own y-axis. Arcs show split reads summed over the libraries of a group, with width "
        "proportional to the junction’s share of reads at its shared splice site. Reads were counted with "
        "the junction rules of Methods 2.5 (MAPQ ≥ 30, primary non-duplicate alignments, 8 bp anchors); "
        "the SH-SY5Y counts are those of the local-splicing-variation test in Section 3.5. Junction a is "
        "the canonical exon 4–exon 5 junction (chr19:1,235,146–1,235,500); b joins exon 4 to the "
        "annotated alternative 3′ splice site at 1,235,342 (1,235,343–1,235,500); c leaves the 31-nt "
        "segment through the unannotated donor at 1,235,312 (1,235,146–1,235,311). The shaded interval "
        "lies between the alternative and canonical 3′ splice sites. Transcript models are from GENCODE "
        "v47; the dashed model is the exon structure implied by junction c. (B) ΔPSI of the 32 "
        "coverage-qualified rMATS CBARP events in five datasets (Supplementary Table S3). The sign refers "
        "to each event’s rMATS-defined inclusion form, which differs between events.")
doc.save(DOCX)
print(DOCX)

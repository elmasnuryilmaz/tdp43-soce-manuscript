"""Move the analyses that cannot carry their own weight into the supplementary material.

The structural scan of 27 September found four passages that occupy main-text space without
supporting the argument:

  the NMD interaction screen, whose own positive controls show no NMD sensitivity (median
  interaction −0.11 log2, p = 0.62), so its negative result for the SOCE genes is uninformative;
  the FRASER outlier analysis, which by its own design cannot work in a cohort where two thirds of
  the samples are the perturbed state;
  the polyadenylation screen in the two mouse lines, for which no positive control exists because
  the STMN2 cryptic polyadenylation site is absent from the mouse gene;
  the depth of the multiple sclerosis analysis, which is a side observation in a paper about ALS.

Each keeps a short statement in the main text, with the full text and every number moved to a new
"Supplementary results" section. Section numbers are left untouched, so no cross-reference breaks.
"""
from pathlib import Path

import docx
from docx import Document

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx"
SUPPLEMENTARY = ROOT / "supplementary/SUPPLEMENTARY_MATERIAL.docx"
src = (ROOT / "code/condense_results_v107.py").read_text()
exec(src[src.index("def replace_span("):src.index("doc = Document(MANUSCRIPT)")])

def replace_between(paragraph, start_marker, end_marker, new, keep_end=True):
    """Replace the text from start_marker to end_marker, reading the span off the paragraph itself.

    This avoids retyping long passages, where a single character would break the match.
    """
    text = paragraph.text
    i = text.index(start_marker)
    j = text.index(end_marker, i) + (len(end_marker) if keep_end else 0)
    replace_span(paragraph, text[i:j], new)


GENES = ["STMN2", "CBARP", "STIM1", "STIM2", "SARAF", "TRPC1", "ATP2A2", "RYR2", "MCU", "Stmn2",
         "Saraf", "Atp2a2", "Trpc1", "MBP", "PLP1", "MOG", "MAG", "GFAP"]

# ----------------------------------------------------------------- main text
doc = Document(MANUSCRIPT)
p = list(doc.paragraphs)
assert p[102].text.startswith("3.6 Exploratory NMD")

# 2.6 FRASER — the citation field stays in the first sentence
replace_span(p[32], " FRASER treats aberrant splicing as a rare, sample-specific deviation from a "
                    "cohort norm; because depleted samples form two-thirds of this cohort, results "
                    "are reported as per-sample outlier burden.",
             " Because that method treats aberrant splicing as a rare deviation from a cohort norm "
             "while two thirds of this cohort are depleted samples, the run is reported as a limit "
             "of the approach and is described in Supplementary Results 1.")

# 2.8 NMD — keep the first sentence with its field, replace the rest
replace_span(p[37], " Sequencing batch is confounded with TDP-43 status in that design, so the "
                    "interaction was formed as a difference of within-batch differences:",
             " Sequencing batch is confounded with TDP-43 status in that design, so the interaction "
             "was formed as a difference of within-batch differences, with the four NMD-inhibition "
             "conditions as the unit of inference. The model, the gene filter and the tests are "
             "given in Supplementary Results 2, and the results are treated as exploratory "
             "throughout.")

# 3.5 FRASER result
replace_span(p[101], "FRASER, run on the nine-sample doxycycline series, returned no genome-wide "
                     "significant outlier and no difference in per-sample outlier burden between "
                     "depleted and control libraries (4.8 versus 2.3 events at p < 10⁻⁵; p = 0.35). "
                     "This is the expected behaviour of an outlier method in a cohort where the "
                     "aberrant state is the majority, and we report it as a limit of that approach "
                     "at this cohort size rather than as evidence against aberrant splicing.",
             "An outlier-based method, FRASER, returned no genome-wide significant outlier and no "
             "difference in per-sample outlier burden between depleted and control libraries, which "
             "is how such a method behaves when the aberrant state is the majority of the cohort "
             "(Supplementary Results 1).")

# 3.6 NMD — four paragraphs become one
replace_span(p[103], "Cryptic exons frequently introduce premature termination codons, which can "
                     "make their transcripts susceptible to nonsense-mediated decay (NMD) and "
                     "under-represented in steady-state RNA. An independent experiment in which "
                     "TDP-43 knockdown is crossed with knockdown of XRN1, UPF1 and SMG6 "
                     "(GSE307054) provides an exploratory interaction test (Supplementary Figure "
                     "S5): a transcript pool produced upon TDP-43 loss and degraded by NMD would be "
                     "expected to rise preferentially when both perturbations are present "
                     "(Methods 2.8).",
             "Cryptic exons frequently introduce premature termination codons, which can make their "
             "transcripts susceptible to nonsense-mediated decay (NMD) and under-represented in "
             "steady-state RNA. An independent experiment in which TDP-43 knockdown is crossed with "
             "knockdown of XRN1, UPF1 and SMG6 (GSE307054) allows this to be tested as an "
             "interaction (Methods 2.8; Supplementary Figure S5; Supplementary Results 2). No gene "
             "passed genome-wide correction in the conservative four-condition analysis (smallest "
             "q = 0.053). CBARP had the largest interaction of the SOCE panel (+1.52 log₂, positive "
             "in all four conditions) but a condition-level q of 0.145, so the screen ranks it "
             "rather than establishing it as an NMD target. The analysis carries little weight in "
             "either direction: the sixteen literature positive controls were not NMD-sensitive as "
             "a group either (median interaction −0.11 log₂; p = 0.62), so this dataset does not "
             "demonstrate the sensitivity that would make a negative result for the SOCE genes "
             "informative.")

# 3.7 mouse lines
replace_between(p[113], ". No SOCE-machinery unit exceeded", "is not a poly(A) site.",
                ", so neither line can say whether an absent signal means an absent event. No "
                "SOCE-machinery unit exceeded |Δ| = 0.30 in either line. One unit is positive in "
                "both motor-neuron models, the SARAF intron 5 index (+0.097 in the iPSC-derived "
                "motor neurons and +0.204 in NSC34, against 0.000 in SH-SY5Y); we record it as a "
                "candidate rather than a finding, for the reasons set out in Supplementary "
                "Results 3, where the mouse values are given in full.")

# 3.9 multiple sclerosis
replace_span(p[125], " Because TRPC1 contributes to SOCE in oligodendrocyte precursor cells",
             " The decrease could reflect demyelination rather than a neuronal change, because "
             "TRPC1 contributes to SOCE in oligodendrocyte precursor cells")
replace_between(p[125], ", we asked whether the decrease simply reflected",
                "(Supplementary Table S16).",
                ". Two controls argue against that explanation: the decrease is present in "
                "normal-appearing white matter, where the myelin markers are unchanged "
                "(δ = −0.482 at sample level, p = 0.005; −0.771 at donor level, p = 0.030), and it "
                "survives regression of TRPC1 on MBP, PLP1 and GFAP, although the adjusted "
                "difference is smaller and, at donor level, no longer reaches significance. The "
                "lesion types, the myelin-adjusted values and the second cohort are given in "
                "Supplementary Results 4 and Supplementary Table S16.")

doc.save(MANUSCRIPT)

# the three paragraphs that are now carried by the supplementary
doc = Document(MANUSCRIPT)
p = list(doc.paragraphs)
for i in sorted([39, 38, 106, 105, 104], reverse=True):
    p[i]._element.getparent().remove(p[i]._element)
doc.save(MANUSCRIPT)
print("main text condensed")

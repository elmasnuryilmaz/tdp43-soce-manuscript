"""Fifth referee round (27 September 2026) — the points that are not already in the text.

Four of the report's points were already answered and are left as they are: the WST-1 wells are
summarised descriptively rather than tested (Methods 2.16), the store is emptied with
cyclopiazonic acid and not thapsigargin, the title and conclusion are associative, and the NMD
analysis is labelled exploratory in Methods 2.8, Section 3.6 and the Limitations.

Applied here:
  effect sizes for the two Fura-2 comparisons and for the per-culture ratio, so the result does
  not rest on the p value alone (Hedges' g, which corrects the small-sample bias of Cohen's d);
  the Figure 1C caption names the WST-1 wells as technical replicates of one experiment;
  the NMD and APA screens are called hypothesis-generating where they are summarised.
"""
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx"
src = (ROOT / "code/condense_results_v107.py").read_text()
exec(src[src.index("def replace_span("):src.index("doc = Document(MANUSCRIPT)")])

doc = Document(MANUSCRIPT)
p = list(doc.paragraphs)

# effect sizes in the Results
replace_span(p[68], "difference −1.30, 95% CI −2.40 to −0.19; Welch’s t-test p = 0.035",
             "difference −1.30, 95% CI −2.40 to −0.19; Welch’s t-test p = 0.035; Hedges’ g = −2.9, "
             "95% CI −5.3 to −0.5")
replace_span(p[68], "difference −0.09, 95% CI −0.30 to +0.13; p = 0.31",
             "difference −0.09, 95% CI −0.30 to +0.13; p = 0.31; Hedges’ g = −0.8, 95% CI −2.2 to +0.6")
replace_span(p[68], "This ratio was 5.88 ± 1.10 in controls and 1.36 ± 0.01 after knockdown, 77% lower,",
             "This ratio was 5.88 ± 1.10 in controls and 1.36 ± 0.01 after knockdown, 77% lower "
             "(Hedges’ g = −2.7, 95% CI −4.9 to −0.4),")

# how the effect sizes were obtained, and what they do and do not add at this sample size
replace_span(p[61], "these two pre-specified comparisons were not adjusted for multiplicity.",
             "these two pre-specified comparisons were not adjusted for multiplicity. Effect sizes "
             "are Hedges’ g, the small-sample-corrected standardised mean difference, with the "
             "interval of Hedges and Olkin; at three samples per group the point estimate is large "
             "but its interval is correspondingly wide, so it describes the size of the difference "
             "rather than establishing it.")

# the WST-1 wells are technical replicates
replace_span(p[67], "(C) WST-1 metabolic signal at 48 h (n = 4; descriptive summary only).",
             "(C) WST-1 metabolic signal at 48 h (n = 4 wells of a single experiment, that is "
             "technical replicates; summarised descriptively and not tested).")

# the two screens are hypothesis-generating
replace_span(p[114], "Taken together, these analyses support altered calcium-related RNA profiles "
                     "but do not establish whether any SOCE-machinery gene is a direct target of "
                     "cryptic splicing, APA or NMD-coupled degradation.",
             "Taken together, these analyses support altered calcium-related RNA profiles but do "
             "not establish whether any SOCE-machinery gene is a direct target of cryptic "
             "splicing, APA or NMD-coupled degradation. The NMD and APA screens are "
             "hypothesis-generating: the first because batch is confounded with TDP-43 status, the "
             "second because it measures coverage gradients rather than poly(A) sites.")

doc.save(MANUSCRIPT)
print("revised:", MANUSCRIPT.name)

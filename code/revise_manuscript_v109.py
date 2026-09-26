"""Apply the fourth referee round (26 September 2026) — the points that need no new experiment.

  M1 (text part)  with three cultures per group no rank test can reach two-sided p < 0.10, so the
                  comparison rests on a parametric assumption; stated in Limitations
  M3 (text part)  the Fura-2 timing is restored to Methods 2.14 and the seeding density is stated,
                  and the cell-number confounder is named in Limitations
  M4, M2          the RT-PCR validation of the CBARP junction and the SOAR exon in the authors' own
                  knockdown cells is named as the first follow-up experiment
  M5              cells per cuvette and what was not measured are stated
  m2              the two meanings of NA in Table 4 are separated and the SH-SY5Y CBARP entry is
                  explained by the read-support and novel-site requirements
  m4              amplification efficiencies were not determined; knockdown was not re-measured in
                  the target-gene RNA set
  m5              GSE307054 is marked as a preprint dataset
  m6              the abstract names the spinal cord result
  m8              the limit of the readdition-to-release ratio is stated

Left for the authors: the new experiments (M1, M2, M4, M6), whether the Fura-2 samples were
adherent or in suspension (M5), and the two further WST-1 experiments (M6).
"""
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx"
src = (ROOT / "code/condense_results_v107.py").read_text()
exec(src[src.index("def replace_span("):src.index("doc = Document(MANUSCRIPT)")])

doc = Document(MANUSCRIPT)
p = list(doc.paragraphs)
assert p[56].text.startswith("Non-targeting shRNA control and shTDP-43 cells were seeded")

# m6 — the abstract names the spinal cord result
replace_span(p[5], "In ALS brain, TRPC1 and SARAF increased and CBARP decreased in six brain regions",
             "In ALS tissue, TRPC1 and SARAF increased and CBARP decreased in six brain regions, and "
             "SARAF and CBARP changed in the same direction in cervical and lumbar cord")

# m5 — the NMD dataset is a preprint
replace_span(p[37], "Gene-level counts were obtained from GSE307054 (Sinha et al., 2025; i3Neurons;",
             "Gene-level counts were obtained from GSE307054, the dataset of a preprint "
             "(Sinha et al., 2025; i3Neurons;")

# m4 — amplification efficiencies and the separate RNA sets
replace_span(p[54], " Reporting guidelines for quantitative PCR recommend at least two validated "
                    "reference genes; the single-reference design is therefore listed among the "
                    "limitations.",
             " Amplification efficiencies were not determined for the individual assays, and the "
             "2^(−ΔΔCt) calculation therefore assumes near-equal efficiencies. Knockdown was "
             "verified in the RNA set collected in April and May and was not re-measured in the "
             "June and July set in which the four targets were quantified. Reporting guidelines "
             "for quantitative PCR recommend at least two validated reference genes, so the "
             "single-reference design and these two points are listed among the limitations.")

# M3 and M5 — timing, seeding and what the cuvette measurement cannot separate
replace_span(p[56], "Non-targeting shRNA control and shTDP-43 cells were seeded in 24-well plates "
                    "at 40,000 cells per well.",
             "Non-targeting shRNA control and shTDP-43 cells were prepared for measurement 72 h "
             "after transduction, while puromycin selection was still in progress; the RNA of "
             "Section 2.13 came from parallel cultures of the same transduction on day 5, after "
             "selection was complete. Both groups were seeded in 24-well plates at the same "
             "density, 40,000 cells per well, and the samples of a group were prepared and "
             "measured on the same day.")
replace_span(p[57], "Free extracellular Ca²⁺ after readdition was not measured independently.",
             "Free extracellular Ca²⁺ after readdition was not measured independently, and neither "
             "the cell number nor the dye loading of each cuvette was recorded, so a difference in "
             "either between the groups cannot be excluded.")

# M1 — what a three-versus-three design can and cannot support
replace_span(p[61], "The readdition-to-release ratio of each culture was summarised descriptively.",
             "With three samples per group the smallest two-sided p value attainable by a rank or "
             "permutation test is 0.10, so these comparisons rest on the parametric assumption of "
             "Welch’s test rather than on a distribution-free one. The readdition-to-release ratio "
             "of each culture was summarised descriptively.")

# m2 — the two meanings of NA, and why SH-SY5Y CBARP is not counted
replace_span(p[259], " iPSC-MN, iPSC-derived motor neurons; KD, knockdown.",
             " In the Tier 1, SOCE-machinery and positive-control columns NA means that no gene of "
             "that set carried a high-confidence call, whereas in the null column it means that the "
             "dataset has fewer than four control replicates. The CBARP junction of Section 3.5 is "
             "not counted for SH-SY5Y: the high-confidence definition requires a novel splice site "
             "carrying at least 20 knockdown reads, and the novel-site junction at this locus "
             "carried six, while the junction that changes most in SH-SY5Y uses two annotated "
             "sites and is already used in controls (17% of exon-4 donor reads). iPSC-MN, "
             "iPSC-derived motor neurons; KD, knockdown.")

# M1, M3, m8 — Limitations
replace_span(p[147], "The Fura-2 comparison rests on three independent cultures per group, and the "
                     "readdition signal was not characterised pharmacologically, for example with "
                     "an ORAI channel inhibitor.",
             "The Fura-2 comparison rests on three cultures per group. At that size no rank or "
             "permutation test can reach a two-sided p below 0.10, so the difference is supported "
             "by a parametric test alone and would be more convincing with five to six cultures per "
             "group. The readdition signal was not characterised pharmacologically, for example "
             "with a store-operated channel blocker such as BTP2, Synta66 or Gd³⁺, so it is defined "
             "by the depletion–readdition protocol rather than by pharmacology.")
replace_span(p[147], "it is not a cell count and cannot separate fewer cells from lower metabolic "
                     "activity per cell.",
             "it is not a cell count and cannot separate fewer cells from lower metabolic activity "
             "per cell. Because metabolic activity was 38.5% lower, a smaller or less healthy cell "
             "population could contribute to the smaller Fura-2 signal; the ratiometric readout is "
             "in principle independent of cell number, but neither cell number nor maximum response "
             "was measured in the cuvettes.")
replace_span(p[147], "Cuvette Fura-2 ratios report net cytosolic accumulation and do not separate "
                     "influx from extrusion or reuptake.",
             "Cuvette Fura-2 ratios report net cytosolic accumulation and do not separate influx "
             "from extrusion or reuptake. Normalising readdition to the release of the same culture "
             "removes differences in store content but not their consequences, because a smaller "
             "release also means weaker STIM activation, so the ratio bounds rather than eliminates "
             "the contribution of the store.")

doc.save(MANUSCRIPT)
print("revised:", MANUSCRIPT.name)

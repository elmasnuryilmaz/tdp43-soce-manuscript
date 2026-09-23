# Data availability

No new sequencing data were generated for this study. Every dataset analysed is public.
`supplementary/S17_dataset_accessions.csv` lists each one with its design, library type and
run-level group assignment for knockdown experiments and group definitions for patient cohorts; the table below is the short version.

## RNA-seq datasets

| Analysis | Dataset | Accession |
|---|---|---|
| expression + splicing + junction | SH-SY5Y, doxycycline-inducible TDP-43 knockdown (primary model) | GSE296712 / PRJNA1256902 |
| expression + splicing + junction | iPSC colonies, TDP-43 shRNA | GSE230647 / PRJNA962064 |
| expression + splicing + junction | iPSC-derived motor neurons (TDP-43, FUS, TAF15 knockdown) | GSE77702 / PRJNA311234 |
| expression + splicing + junction | Mouse striatum, TDP-43 ASO | GSE27394 / GSE27218 |
| expression + splicing + junction | C2C12 and NSC34 | GSE171714 / SRP314028 |
| junction only | K562 total RNA | ENCODE ENCSR372DZW / ENCSR455TNF |
| junction only | K562 poly(A)+ mRNA | ENCODE ENCSR129RWD / ENCSR134JRE |
| NMD interaction | i3Neurons, TDP-43 × NMD inhibition | GSE307054 / PRJNA1235234 |
| patient tissue | NYGC ALS Consortium / Target ALS | GSE153960 / SRP270799 |
| comparison cohort | Alzheimer's disease, fusiform gyrus | GSE125583 |
| comparison cohort | Parkinson's disease, BA9 | GSE68719 |
| comparison cohort | Multiple sclerosis, white matter | GSE138614 |
| comparison cohort | Multiple sclerosis, five brain regions | GSE123496 |

## Resources downloaded at run time

- recount3 junction matrix for SRP270799 (2.55 GB):
  `https://duffel.rail.bio/recount3/human/data_sources/sra/junctions/99/SRP270799/sra.junctions.SRP270799.ALL.MM.gz`
  with the matching `ALL.RR.gz`, `ALL.ID.gz` and `metadata/99/SRP270799/sra.sra.SRP270799.MD.gz`.
  The file name contains `ALL`; without it the URL returns 404. This study is **not** in
  Snaptron's `srav3h` compilation and cannot be queried there.
- GSE307054 gene-level count table and the authors' size factors.
- GENCODE v47 (human) and vM25 (mouse) annotation; GRCh38 primary assembly; the
  transcriptome-aware HISAT2 index `grch38_tran`.

## What is not in this repository

- FASTQ and BAM files (public; re-download with the accessions above).
- Reference genome, annotation and aligner indexes.
- The complete rMATS output (23 GB). `supplementary/S3_rMATS_significant_events.csv.gz`
  contains every event meeting the manuscript thresholds together with its event coordinates,
  form lengths and raw junction counts. Recomputing the pre-filtered FDR across all tested
  events requires the complete rMATS output.
- Laboratory raw instrument files. The values behind every panel of Figure 6 — Ct values,
  per-replicate relative expression, Fura-2 amplitudes, normalised WST-1 signal values, primers and
  the thermal profile — are in `supplementary/S1_laboratory_source_data.xlsx`.

## Paths

The scripts in `code/` carry absolute paths to the two working roots used during the study
(`~/Desktop/MAKALE` for outputs and `~/Desktop/TEZ` plus an external drive for inputs). They
must be edited before the code runs elsewhere; the intended order is in `README.md`.

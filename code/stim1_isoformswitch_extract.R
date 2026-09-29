# STIM1 rows of the IsoformSwitchAnalyzeR/DEXSeq analysis quoted in Supplementary Results 6.
# Needs the R installation that has IsoformSwitchAnalyzeR (R 4.4; see README, Requirements).
# Run from 09_YAYIN_PAKETI:  ISA_DIR=<folder with 03_switch_list_dexseq_tested.rds> Rscript code/stim1_isoformswitch_extract.R
suppressMessages(library(IsoformSwitchAnalyzeR))
isa_dir <- Sys.getenv("ISA_DIR", "/Users/elmas/Desktop/TEZ/output/reanalysis_corrected_full_2026-07-22/isoform_0_vs_75")
x <- readRDS(file.path(isa_dir, "03_switch_list_dexseq_tested.rds"))
f <- x$isoformFeatures
f <- f[f$gene_name == "STIM1" & !is.na(f$isoform_switch_q_value), ]
out <- unique(f[, c("isoform_id", "dIF", "isoform_switch_q_value", "PTC")])
out <- out[order(out$isoform_switch_q_value), ]
write.csv(out, "source_data/STIM1_isoformswitch_DEXSeq.csv", row.names = FALSE)
print(out, digits = 4)

# Sensitivity of the primary SH-SY5Y comparison to unmodelled variation.
# The published result uses ~ grup; here the same counts are refitted with the surrogate
# variables that svaseq estimates, and the two fits are compared.
# Run: Rscript code/svaseq_sensitivity.R   (from 09_YAYIN_PAKETI)
suppressMessages({library(DESeq2); library(sva)})
set.seed(1)

cnt <- read.csv("source_data/gene_counts_full.csv", row.names = 1, check.names = FALSE)
ctrl <- c("SRR33374996", "SRR33375001", "SRR33374995")
kd   <- c("SRR33374999", "SRR33374997", "SRR33375000")
cnt <- round(as.matrix(cnt[, c(ctrl, kd)]))
cnt <- cnt[rowSums(cnt) >= 10, ]                      # the filter of deseq_full.R
cd <- data.frame(row.names = colnames(cnt),
                 grup = factor(c(rep("ctrl", 3), rep("kd", 3)), levels = c("ctrl", "kd")))

dds <- DESeqDataSetFromMatrix(cnt, cd, ~ grup)
dds <- DESeq(dds, quiet = TRUE)
res0 <- results(dds, contrast = c("grup", "kd", "ctrl"))

dds2 <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt, cd, ~ grup))
dat <- counts(dds2, normalized = TRUE)
dat <- dat[rowMeans(dat) > 1, ]
mod  <- model.matrix(~ grup, colData(dds2))
mod0 <- model.matrix(~ 1, colData(dds2))
sv <- svaseq(dat, mod, mod0)
cat("\nestimated surrogate variables:", sv$n.sv, "\n")

if (sv$n.sv > 0) {
  for (i in seq_len(sv$n.sv)) colData(dds2)[[paste0("SV", i)]] <- sv$sv[, i]
  design(dds2) <- as.formula(paste("~", paste(c(paste0("SV", seq_len(sv$n.sv)), "grup"), collapse = " + ")))
  dds2 <- DESeq(dds2, quiet = TRUE)
  res1 <- results(dds2, contrast = c("grup", "kd", "ctrl"))

  keep <- intersect(rownames(res0)[!is.na(res0$padj)], rownames(res1)[!is.na(res1$padj)])
  a <- res0[keep, ]; b <- res1[keep, ]
  cat("genes with a p_adj in both fits:", length(keep), "\n")
  cat(sprintf("log2FC Pearson r  = %.4f\n", cor(a$log2FoldChange, b$log2FoldChange)))
  cat(sprintf("log2FC Spearman r = %.4f\n", cor(a$log2FoldChange, b$log2FoldChange, method = "spearman")))
  cat(sprintf("median |difference| in log2FC = %.4f\n", median(abs(a$log2FoldChange - b$log2FoldChange))))
  de0 <- keep[a$padj < 0.05 & abs(a$log2FoldChange) >= 1]
  de1 <- keep[b$padj < 0.05 & abs(b$log2FoldChange) >= 1]
  cat(sprintf("differentially expressed, published model: %d; with surrogate variables: %d; shared: %d (%.1f%% of the published set)\n",
              length(de0), length(de1), length(intersect(de0, de1)),
              100 * length(intersect(de0, de1)) / length(de0)))
  cat(sprintf("direction agrees for %.1f%% of the shared genes\n",
              100 * mean(sign(a[intersect(de0, de1), "log2FoldChange"]) == sign(b[intersect(de0, de1), "log2FoldChange"]))))

  cat("\ngene      published          with surrogate variables\n")
  for (g in c("TARDBP", "STIM1", "STIM2", "ORAI1", "ORAI2", "ORAI3", "TRPC1", "SARAF",
              "STIMATE", "CBARP", "ATP2A2", "ATP2A3", "STMN2")) {
    if (g %in% keep)
      cat(sprintf("%-9s log2FC %7.3f q %-9.2g log2FC %7.3f q %-9.2g\n", g,
                  a[g, "log2FoldChange"], a[g, "padj"], b[g, "log2FoldChange"], b[g, "padj"]))
  }
  out <- data.frame(gene = keep,
                    log2FC_published = a$log2FoldChange, padj_published = a$padj,
                    log2FC_with_SV = b$log2FoldChange, padj_with_SV = b$padj)
  write.csv(out, "source_data/svaseq_sensitivity_SHSY5Y.csv", row.names = FALSE)
  cat("\nwritten: source_data/svaseq_sensitivity_SHSY5Y.csv\n")
}

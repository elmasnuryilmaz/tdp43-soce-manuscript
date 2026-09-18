suppressMessages(library(DESeq2))
cnt <- read.csv("gene_counts_full.csv", row.names=1, check.names=FALSE)
ctrl <- c("SRR33374996","SRR33375001","SRR33374995"); kd <- c("SRR33374999","SRR33374997","SRR33375000")
cnt <- round(as.matrix(cnt[, c(ctrl, kd)]))
cnt <- cnt[rowSums(cnt) >= 10, ]
cd <- data.frame(row.names=colnames(cnt), grup=factor(c(rep("ctrl",3), rep("kd",3)), levels=c("ctrl","kd")))
dds <- DESeqDataSetFromMatrix(cnt, cd, ~ grup)
dds <- DESeq(dds, quiet=TRUE)
res <- as.data.frame(results(dds, contrast=c("grup","kd","ctrl")))
res$gene <- rownames(res)
write.csv(res[, c("gene","baseMean","log2FoldChange","lfcSE","pvalue","padj")], "DESeq2_ctrl_vs_75_fullmap.csv", row.names=FALSE)
cat("gen sayisi:", nrow(res), "\n")
for (g in c("STIM1","STIM2","ORAI1","ORAI2","ORAI3","TRPC1","SARAF","CBARP","ATP2A2","ATP2A3","MCUB","STIMATE","TARDBP","STMN2")) {
  if (g %in% rownames(res)) cat(sprintf("%-8s log2FC=%7.3f padj=%.3g\n", g, res[g,"log2FoldChange"], res[g,"padj"]))
}

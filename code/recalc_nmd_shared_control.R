#!/usr/bin/env Rscript
# Conservative re-analysis of GSE307054 NMD interaction.
# Unit of inference = the four NMD-inhibition conditions. Replicate-level
# values are averaged within condition; shared CON/TDPKD baselines are not
# counted as independent observations across conditions.
options(stringsAsFactors = FALSE)
base <- "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/nmd_GSE307054"
out <- "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
cnt <- read.csv(gzfile(file.path(base, "GSE307054_counts.csv.gz")), row.names=1,
                check.names=FALSE)
meta <- read.csv(gzfile(file.path(base, "GSE307054_meta.csv.gz")), row.names=1,
                 check.names=FALSE)
common <- intersect(colnames(cnt), rownames(meta))
cnt <- cnt[, common, drop=FALSE]
meta <- meta[common, , drop=FALSE]
norm <- sweep(as.matrix(cnt), 2, meta$sizeFactor, "/")
lg <- log2(norm + 1)

conds <- c("X", "XS", "XU", "US")
int <- matrix(NA_real_, nrow=nrow(lg), ncol=length(conds),
              dimnames=list(rownames(lg), conds))
for (cc in conds) {
  tx <- rowMeans(lg[, paste0("T", cc, "_", 1:2), drop=FALSE])
  cx <- rowMeans(lg[, paste0("C", cc, "_", 1:2), drop=FALSE])
  t0 <- rowMeans(lg[, paste0("TDPKD_", 1:2), drop=FALSE])
  c0 <- rowMeans(lg[, paste0("CON_", 1:2), drop=FALSE])
  int[, cc] <- (tx - t0) - (cx - c0)
}

keep <- rowMeans(norm[, unlist(lapply(conds, function(x) paste0("T", x, "_", 1:2))), drop=FALSE]) >= 10 &
        rowMeans(norm[, paste0("TDPKD_", 1:2), drop=FALSE]) >= 5
int <- int[keep, , drop=FALSE]

mu <- rowMeans(int)
sdv <- apply(int, 1, sd)
tstat <- mu / (sdv / 2)
p_two <- 2 * pt(-abs(tstat), df=3)
p_sign <- sapply(seq_len(nrow(int)), function(i) {
  k <- sum(int[i, ] > 0)
  binom.test(k, 4, p=0.5, alternative="greater")$p.value
})
bh <- function(p) {
  o <- order(p); q <- p[o] * length(p) / seq_along(p)
  q <- rev(cummin(rev(q))); out <- numeric(length(p)); out[o] <- pmin(q, 1); out
}
res <- data.frame(ens=rownames(int), interaction_log2=mu, sd=sdv,
                  n_positive=rowSums(int > 0), p_t4=p_two, q_t4=bh(p_two),
                  p_sign4=p_sign, q_sign4=bh(p_sign),
                  X=int[,"X"], XS=int[,"XS"], XU=int[,"XU"], US=int[,"US"],
                  check.names=FALSE)
res <- res[order(-res$interaction_log2), ]
write.table(res, file.path(out, "NMD_etkilesim_paylasimli_kontrol_t4.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

mapfile <- "/Volumes/10TBElmas/thesis_addendum_2026/ref/gene_id2name_v47.tsv"
if (file.exists(mapfile)) {
  mp <- read.delim(mapfile, header=FALSE, stringsAsFactors=FALSE)
  colnames(mp)[1:2] <- c("ens", "symbol")
  mp$ens <- sub("\\..*$", "", mp$ens)
  res$symbol <- mp$symbol[match(sub("\\..*$", "", res$ens), mp$ens)]
} else res$symbol <- res$ens
write.table(res, file.path(out, "NMD_etkilesim_paylasimli_kontrol_t4_sembol.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

pos <- c("STMN2","UNC13A","HDGFL2","ACTL6B","AGRN","ARHGAP32","PFKP","ATG4B","KALRN","CAMK2B","RSF1","SETD5")
soce <- c("STIM1","STIM2","ORAI1","ORAI2","ORAI3","TRPC1","TRPC3","TRPC4","TRPC5","TRPC6","SARAF","STIMATE","CBARP","ATP2A1","ATP2A2","ATP2A3","MCU","MCUB","MICU1","MICU2","MICU3","CRACR2A","CRACR2B","ITPR1","ITPR2","ITPR3","RYR1","RYR2","RYR3")
sub <- res[res$symbol %in% c(pos, soce), c("symbol","interaction_log2","n_positive","p_t4","q_t4","p_sign4","q_sign4","X","XS","XU","US")]
write.table(sub, file.path(out, "NMD_etkilesim_paylasimli_kontrol_SOCE_pozitif.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Genes tested:", nrow(res), "\n")
cat("Positive controls:\n"); print(sub[sub$symbol %in% pos, ], row.names=FALSE)
cat("SOCE panel:\n"); print(sub[sub$symbol %in% soce, ], row.names=FALSE)

#!/usr/bin/env Rscript
# v3 - corrected APA indices (samtools bedcov -j) for the three datasets that were never
# recomputed after the reference-skip correction: iPSC-MN, C2C12 and NSC34.
# Mirrors recalc_apa_full_core_R.R: coverage indices, complete enumeration of the
# replicate bootstrap, depth filter of >= 3 summed depth in every sample, no p/q values.
options(stringsAsFactors = FALSE)
D <- "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI"
MAN <- read.delim(file.path(D, "kod", "ornekler.tsv"))

GEN <- c("STMN2","UNC13A","HDGFL2","ACTL6B","AGRN","ARHGAP32","PFKP","ATG4B","KALRN",
         "ELAVL3","SETD5","RSF1","GPSM2","POLDIP3","STIM1","STIM2","ORAI1","ORAI2",
         "ORAI3","TRPC1","SARAF","STIMATE","CBARP","ATP2A1","ATP2A2","ATP2A3","MCU",
         "MCUB","MICU1","MICU2","MICU3","CRACR2A","CRACR2B","SELENOK","SELENON",
         "ITPR1","ITPR2","ITPR3","RYR1","RYR2","RYR3")

run_ds <- function(ds) {
  f <- file.path(D, "sonuclar", paste0("apa_bedcov_", ds, "_corrected_full.tsv"))
  if (!file.exists(f)) { cat(ds, ": coverage file missing\n"); return(NULL) }
  d <- read.delim(f, check.names = FALSE)
  samples <- colnames(d)[5:ncol(d)]
  grp <- MAN$group[match(samples, MAN$sample)]
  kd <- which(grp == "KD"); ct <- which(grp == "CTRL")
  cat(sprintf("%s: %d windows, %d KD + %d CTRL\n", ds, nrow(d), length(kd), length(ct)))
  z <- do.call(rbind, strsplit(d$isim, "\\|"))
  colnames(z) <- c("gene", "unit", "window", "strand")
  d <- cbind(d, z)
  d <- d[toupper(d$gene) %in% GEN, ]
  out <- list()
  set.seed(42)
  for (key in unique(paste(d$gene, d$unit, sep = "|"))) {
    g <- d[paste(d$gene, d$unit, sep = "|") == key, ]
    if (!all(c("I5", "I3") %in% g$window) && !all(c("Udist", "Uprox") %in% g$window)) next
    a <- if ("I5" %in% g$window) "I5" else "Udist"
    b <- if ("I3" %in% g$window) "I3" else "Uprox"
    ra <- as.numeric(g[g$window == a, samples]); rb <- as.numeric(g[g$window == b, samples])
    la <- as.numeric(g[g$window == a, "end"] - g[g$window == a, "start"])
    lb <- as.numeric(g[g$window == b, "end"] - g[g$window == b, "start"])
    aa <- ra / la; bb <- rb / lb; idx <- aa / (aa + bb)
    if (any(!is.finite(idx)) || any((aa + bb) < 3)) next
    delta <- mean(idx[kd]) - mean(idx[ct])
    cmb <- expand.grid(c(rep(list(kd), length(kd)), rep(list(ct), length(ct))))
    boot <- apply(cmb, 1, function(q)
      mean(idx[q[seq_along(kd)]]) - mean(idx[q[length(kd) + seq_along(ct)]]))
    out[[length(out) + 1]] <- data.frame(
      dataset = ds, gene = g$gene[1], unit = g$unit[1],
      measure = if (a == "I5") "IPA_index" else "UTR_distal_index",
      index_KD = mean(idx[kd]), index_CTRL = mean(idx[ct]), delta = delta,
      boot_low = unname(quantile(boot, .025)), boot_high = unname(quantile(boot, .975)),
      n_boot = nrow(cmb),
      sample_KD = paste(round(idx[kd], 4), collapse = ";"),
      sample_CTRL = paste(round(idx[ct], 4), collapse = ";"))
  }
  if (!length(out)) { cat(ds, ": no unit passed the depth filter\n"); return(NULL) }
  do.call(rbind, out)
}

args <- commandArgs(trailingOnly = TRUE)
dsets <- if (length(args)) args else c("iPSC_MN", "C2C12", "NSC34")
res <- do.call(rbind, lapply(dsets, run_ds))
o <- file.path(D, "sonuclar", if (length(args)) paste0("APA_corrected_", paste(dsets, collapse="_"), ".tsv")
                              else "APA_corrected_remaining_datasets.tsv")
write.table(res, o, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nwritten:", o, "-", nrow(res), "units\n\n")
res$abs <- abs(res$delta)
print(res[order(res$dataset, -res$abs),
          c("dataset","gene","unit","measure","index_KD","index_CTRL","delta",
            "boot_low","boot_high")], row.names = FALSE, digits = 3)
cat("\n-- units whose bootstrap interval excludes zero --\n")
sig <- res[sign(res$boot_low) == sign(res$boot_high), ]
print(sig[order(sig$dataset, -sig$abs),
          c("dataset","gene","unit","measure","delta","boot_low","boot_high")],
      row.names = FALSE, digits = 3)

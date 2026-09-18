#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
f <- "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar/apa_bedcov_SH_SY5Y_target_j.tsv"
d <- read.delim(f, check.names=FALSE)
samples <- colnames(d)[5:ncol(d)]
z <- do.call(rbind, strsplit(d$isim, "\\|", fixed=FALSE))
colnames(z) <- c("gene","unit","window","strand")
d <- cbind(d,z)
out <- list(); set.seed(42)
for (key in unique(paste(d$gene,d$unit,sep="|"))) {
  ix <- which(paste(d$gene,d$unit,sep="|")==key); g <- d[ix,]
  if (!all(c("I5","I3") %in% g$window) && !all(c("Udist","Uprox") %in% g$window)) next
  a <- if ("I5" %in% g$window) "I5" else "Udist"
  b <- if ("I3" %in% g$window) "I3" else "Uprox"
  ra <- as.numeric(g[g$window==a, samples]); rb <- as.numeric(g[g$window==b, samples])
  la <- as.numeric(g[g$window==a,"end"]-g[g$window==a,"start"])
  lb <- as.numeric(g[g$window==b,"end"]-g[g$window==b,"start"])
  aa <- ra/la; bb <- rb/lb; idx <- aa/(aa+bb)
  if (any(!is.finite(idx)) || any((aa+bb) < 3)) next
  delta <- mean(idx[1:3],na.rm=TRUE)-mean(idx[4:6],na.rm=TRUE)
  boot <- replicate(1000, mean(sample(idx[1:3],3,replace=TRUE),na.rm=TRUE)-mean(sample(idx[4:6],3,replace=TRUE),na.rm=TRUE))
  out[[length(out)+1]] <- data.frame(gene=g$gene[1], unit=g$unit[1],
    measure=if(a=="I5")"IPA_index" else "UTR_distal_index",
    index_KD=mean(idx[1:3]), index_CTRL=mean(idx[4:6]), delta=delta,
    boot_low=quantile(boot,.025,na.rm=TRUE), boot_high=quantile(boot,.975,na.rm=TRUE),
    sample_KD=paste(round(idx[1:3],4),collapse=";"),
    sample_CTRL=paste(round(idx[4:6],4),collapse=";"))
}
r <- do.call(rbind,out)
write.table(r, "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar/APA_target_j_SH_SY5Y_summary.tsv", sep="\t", quote=FALSE, row.names=FALSE)
print(r[order(r$gene,r$unit),], row.names=FALSE)

#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
f <- "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar/apa_bedcov_SH_SY5Y_corrected_full.tsv"
d <- read.delim(f, header=FALSE, check.names=FALSE)
colnames(d) <- c("chrom","start","end","isim","SRR33374999","SRR33374997","SRR33375000","SRR33374996","SRR33375001","SRR33374995")
samples <- colnames(d)[5:10]
z <- do.call(rbind, strsplit(d$isim, "\\|"))
colnames(z) <- c("gene","unit","window","strand")
d <- cbind(d,z)
genes <- c("STMN2","UNC13A","HDGFL2","ACTL6B","AGRN","ARHGAP32","PFKP","ATG4B","KALRN","ELAVL3","SETD5","RSF1","GPSM2","POLDIP3","STIM1","STIM2","ORAI1","ORAI2","ORAI3","TRPC1","SARAF","STIMATE","CBARP","ATP2A1","ATP2A2","ATP2A3","MCU","MCUB","MICU1","MICU2","MICU3","CRACR2A","CRACR2B","SELENOK","SELENON","ITPR1","ITPR2","ITPR3","RYR1","RYR2","RYR3")
d <- d[d$gene %in% genes,]
out <- list(); set.seed(42)
for (key in unique(paste(d$gene,d$unit,sep="|"))) {
  g <- d[paste(d$gene,d$unit,sep="|")==key,]
  if (!all(c("I5","I3") %in% g$window) && !all(c("Udist","Uprox") %in% g$window)) next
  a <- if ("I5" %in% g$window) "I5" else "Udist"; b <- if ("I3" %in% g$window) "I3" else "Uprox"
  ra <- as.numeric(g[g$window==a,samples]); rb <- as.numeric(g[g$window==b,samples])
  la <- as.numeric(g[g$window==a,"end"]-g[g$window==a,"start"]); lb <- as.numeric(g[g$window==b,"end"]-g[g$window==b,"start"])
  aa <- ra/la; bb <- rb/lb; idx <- aa/(aa+bb)
  if (any(!is.finite(idx)) || any((aa+bb)<3)) next
  delta <- mean(idx[1:3])-mean(idx[4:6])
  # Fast exact enumeration of the 3^6 bootstrap combinations (729 draws).
  cmb <- expand.grid(i1=1:3,i2=1:3,i3=1:3,j1=4:6,j2=4:6,j3=4:6)
  boot <- apply(cmb,1,function(q) mean(idx[q[1:3]])-mean(idx[q[4:6]]))
  out[[length(out)+1]] <- data.frame(gene=g$gene[1],unit=g$unit[1],measure=if(a=="I5")"IPA_index" else "UTR_distal_index",index_KD=mean(idx[1:3]),index_CTRL=mean(idx[4:6]),delta=delta,boot_low=quantile(boot,.025),boot_high=quantile(boot,.975),sample_KD=paste(round(idx[1:3],4),collapse=";"),sample_CTRL=paste(round(idx[4:6],4),collapse=";"))
}
r <- do.call(rbind,out)
write.table(r,"/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar/APA_corrected_full_core_summary.tsv",sep="\t",quote=FALSE,row.names=FALSE)
print(r[order(r$gene,r$unit),],row.names=FALSE)

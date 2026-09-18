#!/bin/bash
# BAM -> junction sayımları (regtools junctions extract eşdeğeri, saf samtools+awk)
# Çıktı: chrom \t intron_start(1-based) \t intron_end(1-based) \t strand \t count
# Filtreler: MAPQ>=30, primary+mapped, min anchor 8 bp, intron 50..500000 bp
set -euo pipefail
BAM="$1"; OUT="$2"; REGION="${3:-}"
samtools view -q 30 -F 3844 "$BAM" $REGION 2>/dev/null | awk -v OFS='\t' '
{
  pos=$4; cig=$6
  if (cig !~ /N/) next
  strand="."
  for (i=12;i<=NF;i++) if ($i ~ /^XS:A:/) { strand=substr($i,6,1); break }
  n=split(cig, ops, /[0-9]+/); m=split(cig, lens, /[MIDNSHP=X]/)
  ref=pos; anchor=0; delete J
  nj=0
  for (i=1;i<=m;i++) {
    L=lens[i]+0; O=ops[i+1]
    if (O=="M"||O=="="||O=="X") { ref+=L; anchor=L }
    else if (O=="D") { ref+=L }
    else if (O=="N") {
      if (L>=50 && L<=500000 && anchor>=8) { nj++; JS[nj]=ref; JE[nj]=ref+L-1 }
      else { nj=-1000 }
      ref+=L; anchor=0
    }
  }
  if (nj<=0) next
  # son blok anchor kontrolü: CIGAR sonundaki M uzunluğu
  lastM=0
  for (i=1;i<=m;i++) { O=ops[i+1]; if (O=="M"||O=="="||O=="X") lastM=lens[i]+0 }
  if (lastM<8) next
  for (k=1;k<=nj;k++) print $3, JS[k], JE[k], strand
}' | sort -k1,1 -k2,2n -k3,3n -k4,4 | uniq -c | awk -v OFS='\t' '{print $2,$3,$4,$5,$1}' > "$OUT"

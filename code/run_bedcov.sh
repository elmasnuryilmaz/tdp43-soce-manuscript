#!/bin/bash
set -uo pipefail
D=/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI
awk -F'\t' 'NR>1 && $4=="human"{print $1"\t"$2"\t"$3"\t"$5}' $D/kod/ornekler.tsv > $D/kod/human_bams.tsv
for ds in SH_SY5Y iPSC_koloni iPSC_MN K562_mRNA K562_totalRNA; do
  bams=(); smps=()
  while IFS=$'\t' read -r d s g b; do
    [ "$d" = "$ds" ] || continue
    bams+=("$b"); smps+=("$s")
  done < $D/kod/human_bams.tsv
  { printf "chrom\tstart\tend\tisim"; printf "\t%s" "${smps[@]}"; printf "\n"; } > $D/sonuclar/apa_bedcov_${ds}.tsv
  samtools bedcov -Q 30 -j "$D/kod/apa_pencereleri_human.bed" "${bams[@]}" >> $D/sonuclar/apa_bedcov_${ds}.tsv
  echo "$ds: $(( $(wc -l < $D/sonuclar/apa_bedcov_${ds}.tsv) - 1 )) pencere x ${#bams[@]} örnek"
done

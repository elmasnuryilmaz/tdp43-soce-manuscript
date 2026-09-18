#!/bin/bash
set -uo pipefail
D=/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI
for ds in NSC34 C2C12; do
  bams=(); smps=()
  while IFS=$'\t' read -r d s g sp b; do
    [ "$d" = "$ds" ] || continue
    bams+=("$b"); smps+=("$s")
  done < <(tail -n +2 $D/kod/ornekler.tsv)
  { printf "chrom\tstart\tend\tisim"; printf "\t%s" "${smps[@]}"; printf "\n"; } > $D/sonuclar/apa_bedcov_${ds}.tsv
  samtools bedcov -Q 30 -j "$D/kod/apa_pencereleri_mouse_chr.bed" "${bams[@]}" 2>/dev/null >> $D/sonuclar/apa_bedcov_${ds}.tsv
  echo "$ds: $(( $(wc -l < $D/sonuclar/apa_bedcov_${ds}.tsv) - 1 )) pencere x ${#bams[@]} örnek"
done
echo "BEDCOV4 TAMAM"

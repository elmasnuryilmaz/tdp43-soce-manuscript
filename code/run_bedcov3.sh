#!/bin/bash
set -uo pipefail
D=/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI
run_one () {
  ds=$1; bed=$2
  bams=(); smps=()
  while IFS=$'\t' read -r d s g sp b; do
    [ "$d" = "$ds" ] || continue
    bams+=("$b"); smps+=("$s")
  done < <(tail -n +2 $D/kod/ornekler.tsv)
  [ ${#bams[@]} -eq 0 ] && return
  { printf "chrom\tstart\tend\tisim"; printf "\t%s" "${smps[@]}"; printf "\n"; } > $D/sonuclar/apa_bedcov_${ds}.tsv
  samtools bedcov -Q 30 -j "$bed" "${bams[@]}" >> $D/sonuclar/apa_bedcov_${ds}.tsv
  echo "$ds: $(( $(wc -l < $D/sonuclar/apa_bedcov_${ds}.tsv) - 1 )) pencere x ${#bams[@]} örnek"
}
run_one iPSC_MN $D/kod/apa_pencereleri_human.bed
run_one NSC34   $D/kod/apa_pencereleri_mouse_chr.bed
run_one C2C12   $D/kod/apa_pencereleri_mouse_chr.bed
echo "BEDCOV3 TAMAM"

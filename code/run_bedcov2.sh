#!/bin/bash
set -uo pipefail
D=/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI
run_one () {           # $1=dataset $2=bed $3=tur(human/mouse)
  ds=$1; bed=$2; tur=$3
  bams=(); smps=()
  while IFS=$'\t' read -r d s g sp b; do
    [ "$d" = "$ds" ] || continue
    bams+=("$b"); smps+=("$s")
  done < <(tail -n +2 $D/kod/ornekler.tsv)
  [ ${#bams[@]} -eq 0 ] && { echo "$ds: örnek yok"; return; }
  { printf "chrom\tstart\tend\tisim"; printf "\t%s" "${smps[@]}"; printf "\n"; } > $D/sonuclar/apa_bedcov_${ds}.tsv
  samtools bedcov -Q 30 -j "$bed" "${bams[@]}" >> $D/sonuclar/apa_bedcov_${ds}.tsv
  echo "$ds: $(( $(wc -l < $D/sonuclar/apa_bedcov_${ds}.tsv) - 1 )) pencere x ${#bams[@]} örnek"
}
for ds in SH_SY5Y iPSC_koloni iPSC_MN K562_mRNA K562_totalRNA; do
  run_one $ds $D/kod/apa_pencereleri_human.bed human
done
for ds in C2C12 NSC34; do
  run_one $ds $D/kod/apa_pencereleri_mouse.bed mouse
done
echo "BEDCOV TAMAM"

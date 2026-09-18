#!/bin/bash
# v3 - APA coverage with the mandatory -j flag for the three datasets that were never
# recomputed after the correction: iPSC-MN (human windows) and C2C12 / NSC34 (mouse windows).
# Output goes to *_corrected_full.tsv, matching the SH-SY5Y naming; the pre-correction
# files are left untouched.
set -uo pipefail
D=/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI
SAM=/opt/homebrew/bin/samtools

run_one () {
  ds=$1; bed=$2
  bams=(); smps=()
  while IFS=$'\t' read -r d s g sp b; do
    [ "$d" = "$ds" ] || continue
    bams+=("$b"); smps+=("$s")
  done < <(tail -n +2 $D/kod/ornekler.tsv)
  [ ${#bams[@]} -eq 0 ] && { echo "$ds: ornek yok"; return; }
  out=$D/sonuclar/apa_bedcov_${ds}_corrected_full.tsv
  { printf "chrom\tstart\tend\tisim"; printf "\t%s" "${smps[@]}"; printf "\n"; } > "$out"
  echo "$ds: ${#bams[@]} BAM, basladi $(date +%H:%M:%S)"
  $SAM bedcov -Q 30 -j "$bed" "${bams[@]}" >> "$out"
  echo "$ds: $(( $(wc -l < "$out") - 1 )) pencere x ${#bams[@]} ornek, bitti $(date +%H:%M:%S)"
}

run_one iPSC_MN $D/kod/apa_pencereleri_human.bed
run_one C2C12   $D/kod/apa_pencereleri_mouse_chr.bed
run_one NSC34   $D/kod/apa_pencereleri_mouse_chr.bed
echo "BEDCOV v3 TAMAM"

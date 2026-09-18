#!/bin/bash
D=/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI
while IFS=$'\t' read -r ds smp grp sp bam; do
  [ "$ds" = "dataset" ] && continue
  out="$D/junctions/${smp}.junc"
  if [ -s "$out" ]; then echo "[SKIP] $smp"; continue; fi
  $D/kod/run_junctions.sh "$bam" "$out" 9
done < $D/kod/ornekler_oncelik.tsv
echo "TAMAMLANDI $(date)"

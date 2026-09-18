#!/bin/bash
# Bir BAM icin kromozom-paralel junction cikarimi
set -uo pipefail
BAM="$1"; OUT="$2"; NPROC="${3:-8}"
EX=/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/kod/extract_junctions.sh
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
CHRS=$(samtools idxstats "$BAM" | awk '$1!="*" && $3>0 {print $1}' | grep -Ev "_|^GL|^KI|^MT$|^chrM$" )
i=0
for c in $CHRS; do
  "$EX" "$BAM" "$TMP/$c.j" "$c" &
  i=$((i+1))
  if (( i % NPROC == 0 )); then wait; fi
done
wait
cat "$TMP"/*.j 2>/dev/null | sort -k1,1 -k2,2n -k3,3n > "$OUT"
echo "$(basename $BAM): $(wc -l < "$OUT") junction"

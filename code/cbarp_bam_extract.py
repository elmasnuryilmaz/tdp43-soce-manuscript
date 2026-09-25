"""CBARP lokusunda BAM'den kapsam ve birleşim (junction) sayımı.

Kurallar makaledeki kendi çıkarıcımızla aynıdır (Methods 2.5; code/extract_junctions.sh):
MAPQ >= 30, birincil ve eşlenmiş, kopya/QC-başarısız/ek hizalama değil (-F 3844),
her N öncesindeki M bloğu >= 8 bp, son M bloğu >= 8 bp, intron 50..500.000 bp.
Koordinatlar 1-tabanlı, intron kapsayıcıdır (start = ilk intronik baz).

Kapsam: aynı okuma süzgeciyle, baz başına derinlik; N (referans atlama) ve D
kapsam saymaz. Ayrıca kütüphane şeridine göre yalnız CBARP (eksi) şeridinden gelen
okumalarla ikinci bir kapsam izi üretilir (antisens CBARP-DT kontrolü).

Çalıştırma: ~/anaconda3/envs/tdp43_pipe/bin/python 01_cbarp_bam_extract.py
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd
import pysam

OUT = Path(__file__).resolve().parents[1] / "veri"
OUT.mkdir(exist_ok=True)

CHROM = "19"
START, END = 1234900, 1235650          # 1-tabanlı, kapsayıcı çizim penceresi
JWIN = (1234000, 1236200)              # birleşim arama penceresi

SH = "/Volumes/10TBElmas/tdp43_pipeline_SHSY5Y/bam/{}.sorted.bam"
IP = "/Volumes/10TBElmas/tdp43_pipeline_iPSC_MN/bam/{}.sorted.bam"
SAMPLES = [
    # veri seti, grup, örnek, yol, kütüphane şeridi
    *[("SH-SY5Y", "Control", s, SH.format(s), "fr-secondstrand")
      for s in ("SRR33374996", "SRR33375001", "SRR33374995")],
    *[("SH-SY5Y", "TDP-43 KD", s, SH.format(s), "fr-secondstrand")
      for s in ("SRR33374999", "SRR33374997", "SRR33375000")],
    *[("iPSC colonies", "Control", s, IP.format(s), "fr-firststrand")
      for s in ("SRR24314006", "SRR24314016", "SRR24314019", "SRR24314022")],
    *[("iPSC colonies", "TDP-43 KD", s, IP.format(s), "fr-firststrand")
      for s in ("SRR24314005", "SRR24314015", "SRR24314018", "SRR24314021")],
]
BAD_FLAGS = 3844  # unmapped, secondary, QC fail, duplicate, supplementary
M_OPS = {0, 7, 8}  # M, =, X


def passes(read):
    return read.mapping_quality >= 30 and not (read.flag & BAD_FLAGS)


def junctions(read):
    """extract_junctions.sh ile aynı mantık; geçersizse boş liste."""
    ref = read.reference_start + 1  # 1-tabanlı
    anchor, out, last_m = 0, [], 0
    for op, length in read.cigartuples:
        if op in M_OPS:
            ref += length; anchor = length; last_m = length
        elif op == 2:  # D
            ref += length
        elif op == 3:  # N
            if 50 <= length <= 500000 and anchor >= 8:
                out.append((ref, ref + length - 1))
            else:
                return []
            ref += length; anchor = 0
    if not out or last_m < 8:
        return []
    return out


def transcript_minus(read, libtype):
    """Okumanın temsil ettiği transkript eksi şeritte mi?"""
    r1 = read.is_read1 or not read.is_paired
    same = read.is_reverse if r1 else not read.is_reverse
    # fr-secondstrand: okuma1 transkriptle aynı şeritte
    if libtype == "fr-secondstrand":
        return same
    return not same


rows, cov = [], {}
n = END - START + 1
for dataset, group, sample, path, lib in SAMPLES:
    bam = pysam.AlignmentFile(path)
    depth_all = np.zeros(n, dtype=np.int32)
    depth_minus = np.zeros(n, dtype=np.int32)
    jc = {}
    mapped = bam.mapped
    for read in bam.fetch(CHROM, JWIN[0] - 1, JWIN[1]):
        if not passes(read):
            continue
        for s, e in junctions(read):
            if JWIN[0] <= s and e <= JWIN[1]:
                jc[(s, e)] = jc.get((s, e), 0) + 1
        minus = transcript_minus(read, lib)
        for bs, be in read.get_blocks():  # 0-tabanlı yarı açık; N ve D dışarıda
            a, b = max(bs + 1, START), min(be, END)
            if a <= b:
                depth_all[a - START:b - START + 1] += 1
                if minus:
                    depth_minus[a - START:b - START + 1] += 1
    bam.close()
    cov[sample] = {"all": depth_all.tolist(), "minus": depth_minus.tolist()}
    for (s, e), c in sorted(jc.items()):
        rows.append(dict(dataset=dataset, group=group, sample=sample, chrom=CHROM,
                         start=s, end=e, count=c, mapped_reads=mapped))
    print(sample, group, "junctions:", len(jc), "max depth:", depth_all.max(),
          "antisense share:", round(1 - depth_minus.sum() / max(depth_all.sum(), 1), 4))

pd.DataFrame(rows).to_csv(OUT / "CBARP_junction_counts_per_library.tsv", sep="\t", index=False)
meta = {"chrom": CHROM, "start": START, "end": END,
        "samples": [dict(dataset=d, group=g, sample=s, libtype=l) for d, g, s, _, l in SAMPLES]}
(OUT / "CBARP_coverage_per_library.json").write_text(json.dumps({"meta": meta, "coverage": cov}))
print("yazıldı:", OUT)

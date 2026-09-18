#!/usr/bin/env python3
"""A12 — Kontrol-karşı-kontrol boş test: kriptik çağrıların yanlış pozitif oranı.

Aynı gruptaki kontrol replikaları ikiye bölünüp analiz aynen tekrarlanır.
Gerçek bir biyolojik fark olmadığı için burada bulunan her "kriptik olay"
yanlış pozitiftir. Az replikalı veri setlerinde bu oran yüksektir.
"""
import os, sys, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import lib_junc as L
from A1_kriptik_lsv import lsv_frame
from A1b_lsv_regtools import matris

OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
TAB = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/tablolar"
MAN = pd.read_csv("/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/kod/ornekler_regtools.tsv", sep="\t")

rows = []
ann = {}
for ds in ["iPSC_koloni", "K562_totalRNA", "Fare_striatum", "SH_SY5Y"]:
    sub = MAN[MAN["dataset"] == ds]
    sp = sub["species"].iloc[0]
    ctrl = sub[sub["group"] == "CTRL"].drop_duplicates("replika")["replika"].tolist()
    if len(ctrl) < 4:
        print(f"[{ds}] kontrol replikası {len(ctrl)} — bölünemez, atlandı"); continue
    yari = len(ctrl) // 2
    sahte = dict(zip(ctrl, ["KD"] * yari + ["CTRL"] * (len(ctrl) - yari)))
    sub2 = sub[sub["replika"].isin(ctrl)].copy()
    jm = matris(sub2)
    jm = jm[jm[ctrl].max(axis=1) >= 5].reset_index(drop=True)
    if sp not in ann:
        ann[sp] = L.load_annotation(sp)
    ai, genes = ann[sp]
    jm = L.annotate_junctions(jm, ai, genes)
    D = lsv_frame(jm, ctrl, [sahte[c] for c in ctrl])
    kr = D[(D["sinif"] != "anotasyonlu") & (D["dPSI"] >= 0.05) & (D["q"] < 0.05)
           & (D["psi_CTRL"] <= 0.05) & (D["GA_alt"] > 0)]
    gercek = pd.read_csv(f"{OUT}/RT_KRIPTIK_{ds}.tsv", sep="\t")
    oran = len(kr) / max(len(gercek), 1)
    rows.append(dict(veri_seti=ds, n_kontrol=len(ctrl), bolme=f"{yari}v{len(ctrl)-yari}",
                     null_kriptik_olay=len(kr), null_kriptik_gen=kr["gene"].nunique(),
                     gercek_kriptik_olay=len(gercek), gercek_kriptik_gen=gercek["gene"].nunique(),
                     yanlis_pozitif_orani=round(oran, 3)))
    kr.to_csv(f"{OUT}/NULL_KRIPTIK_{ds}.tsv", sep="\t", index=False)
    print(f"[{ds}] boş test kriptik olay {len(kr)} ({kr['gene'].nunique()} gen) — "
          f"gerçek {len(gercek)} ({gercek['gene'].nunique()} gen), oran {oran:.3f}", flush=True)
N = pd.DataFrame(rows)
N.to_csv(f"{TAB}/S15_null_kontrol_yanlis_pozitif.tsv", sep="\t", index=False)
print("\n", N.to_string(index=False))

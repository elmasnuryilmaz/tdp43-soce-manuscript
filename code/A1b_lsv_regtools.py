#!/usr/bin/env python3
"""A1b — Aynı LSV analizi, diskteki regtools birleşim dosyaları üzerinde.

Kapsam A1'den geniştir: 11 karşılaştırma (SH-SY5Y iki doz, iPSC koloni,
iPSC-MN TDP-43/FUS/TAF15, K562 iki kütüphane, C2C12, NSC34, fare striatum).
Teknik dizileme parçaları biyolojik replika etiketine (replika sütunu) göre
havuzlanır.
"""
import os, sys, time, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import lib_junc as L
from A1_kriptik_lsv import lsv_frame

JD = "/Volumes/10TBElmas/thesis_addendum_2026/01_junctions"
OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
MAN = pd.read_csv("/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/kod/ornekler_regtools.tsv", sep="\t")

def matris(sub):
    """Replika düzeyinde birleşim matrisi (teknik parçalar toplanır)."""
    frames = {}
    for rep, g in sub.groupby("replika", sort=False):
        acc = None
        for _, r in g.iterrows():
            p = f"{JD}/{r['dizin']}/{r['sample']}.junc.bed"
            d = L.load_junc_regtools(p).set_index(["chrom", "start", "end"])["count"]
            acc = d if acc is None else acc.add(d, fill_value=0)
        frames[rep] = acc
    mat = pd.concat(frames, axis=1).fillna(0).astype(int)
    return mat.reset_index()

def main(only=None):
    ann = {}
    ozet = []
    for ds, sub in MAN.groupby("dataset", sort=False):
        if only and ds not in only:
            continue
        t0 = time.time()
        sp = sub["species"].iloc[0]
        reps = sub.drop_duplicates("replika")[["replika", "group"]]
        samples = list(reps["replika"]); groups = list(reps["group"])
        jm = matris(sub)
        jm = jm[jm[samples].max(axis=1) >= 5].reset_index(drop=True)
        if sp not in ann:
            ann[sp] = L.load_annotation(sp)
        ai, genes = ann[sp]
        jm = L.annotate_junctions(jm, ai, genes)
        vc = jm["sinif"].value_counts().to_dict()
        D = lsv_frame(jm, samples, groups)
        D.to_csv(f"{OUT}/RT_LSV_{ds}.tsv.gz", sep="\t", index=False)
        sig = D[(D["q"] < 0.05) & (D["dPSI"].abs() >= 0.10)]
        yeni = sig[sig["sinif"] != "anotasyonlu"]
        kr = D[(D["sinif"] != "anotasyonlu") & (D["dPSI"] >= 0.05) & (D["q"] < 0.05)
               & (D["psi_CTRL"] <= 0.05) & (D["GA_alt"] > 0)]
        kr.to_csv(f"{OUT}/RT_KRIPTIK_{ds}.tsv", sep="\t", index=False)
        print(f"[{ds}] n={len(samples)} replika | birleşim {len(jm)} "
              f"(anot {vc.get('anotasyonlu',0)}, yeni {len(jm)-vc.get('anotasyonlu',0)}) | "
              f"LSV {len(D)}, anlamlı {len(sig)}, yeni-anlamlı {len(yeni)}, "
              f"kriptik {len(kr)} ({kr['gene'].nunique()} gen) [{time.time()-t0:.0f}s]", flush=True)
        ozet.append(dict(veri_seti=ds, n_replika=len(samples), birlesim=len(jm),
                         anotasyonsuz=len(jm)-vc.get('anotasyonlu',0), lsv=len(D),
                         anlamli=len(sig), yeni_anlamli=len(yeni),
                         kriptik_olay=len(kr), kriptik_gen=kr["gene"].nunique()))
    pd.DataFrame(ozet).to_csv(f"{OUT}/RT_LSV_ozet.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main(only=sys.argv[1:] or None)

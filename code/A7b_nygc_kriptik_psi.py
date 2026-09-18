#!/usr/bin/env python3
"""A7b — NYGC/Target ALS kohortunda BİRLEŞİM DÜZEYİ kriptik STMN2/UNC13A ölçümü.

Makalenin en büyük kısıtı buydu: hasta dokusunda TDP-43 işlev kaybı gen düzeyi
STMN2 ile vekillenmişti; bu ölçü yığın dokuda nöronal içerikle karıştığı için
test bilgilendirici olmamıştı. Kriptik birleşim PSI'si aynı gen içindeki bir
ORAN olduğu için nöronal içerikten büyük ölçüde bağımsızdır.

Kriptik koordinatlar bu çalışmanın kendi SH-SY5Y verisinden de novo bulunmuştur:
  STMN2  chr8:79.611.215 -> 79.616.821  (kanonik karşılığı 79.636.801)
  UNC13A chr19:17.641.557 -> 17.642.413 (kanonik karşılığı 17.642.844)

Veri: recount3 SRP270799 (GSE153960) birleşim matrisi, 2.256 kütüphane.
"""
import os, sys, gzip, subprocess, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
import lib_junc as L

ND = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/nygc_junction"
OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
TRIPLE = f"{ND}/hedef_uclu.tsv"
KRIPTIK = {"STMN2": dict(start=79611215, kriptik=79616821, kanonik=79636801),
           "UNC13A": dict(start=17641557, kriptik=17642413, kanonik=17642844)}

# ------------------------------------------------------------ 1. MM çıkarımı
rr = pd.read_csv(f"{ND}/hedef_satirlar.tsv", sep="\t")
if not os.path.exists(TRIPLE) or os.path.getsize(TRIPLE) < 1000:
    print("MM taranıyor (616 M kayıt)...", flush=True)
    want = f"{ND}/_want_rows.txt"
    rr[["satir"]].to_csv(want, index=False, header=False)
    cmd = (f"gunzip -c {ND}/SRP270799.ALL.MM.gz | "
           f"awk 'NR==FNR{{w[$1];next}} /^%/{{next}} NF==3 && ($1 in w)"
           f"{{print $1\"\\t\"$2\"\\t\"$3}}' {want} - > {TRIPLE}")
    subprocess.run(["bash", "-c", cmd], check=True)
trip = pd.read_csv(TRIPLE, sep="\t", header=None, names=["satir", "sutun", "n"])
print(f"üçlü kayıt: {len(trip):,}", flush=True)

# ------------------------------------------------------------ 2. örnek eşlemesi
ids = [int(x) for x in gzip.open(f"{ND}/SRP270799.ALL.ID.gz", "rt").read().split()[1:]]
md = pd.read_csv(f"{ND}/SRP270799.MD.gz", sep="\t", dtype=str)
md["rail_id"] = md["rail_id"].astype(int)
attr = md["sample_attributes"].fillna("")
def alan(s, k):
    for p in s.split("|"):
        if p.startswith(k + ";;"): return p.split(";;", 1)[1]
    return ""
md["grup"] = [alan(s, "group") for s in attr]
md["doku"] = [alan(s, "tissue") for s in attr]
md["denek"] = [alan(s, "subject id") for s in attr]
md["cgnd"] = md["sample_title"].str.extract(r"\[(CGND-[^\]]+)\]")[0]
col2rail = pd.Series(ids, index=np.arange(1, len(ids) + 1))
trip["rail"] = trip["sutun"].map(col2rail)
info = rr.set_index("satir")[["chrom", "start", "end", "gene", "annotated"]]
trip = trip.join(info, on="satir")

# ------------------------------------------------------------ 3. kriptik PSI
psi = {}
for gen, k in KRIPTIK.items():
    sub = trip[(trip["gene"] == gen) & (trip["start"] == k["start"])]
    tot = sub.groupby("rail")["n"].sum()
    kr = sub[sub["end"] == k["kriptik"]].groupby("rail")["n"].sum().reindex(tot.index).fillna(0)
    kanon = sub[sub["end"] == k["kanonik"]].groupby("rail")["n"].sum().reindex(tot.index).fillna(0)
    psi[gen] = pd.DataFrame({f"{gen}_kriptik": kr, f"{gen}_kanonik": kanon,
                             f"{gen}_toplam": tot,
                             f"{gen}_PSI": kr / tot.replace(0, np.nan)})
    print(f"{gen}: {len(tot)} kütüphanede verici okundu; kriptik okuma toplamı {kr.sum():.0f}")
P = pd.concat(psi.values(), axis=1)
S = md.set_index("rail_id").join(P, how="inner")
S.to_csv(f"{OUT}/NYGC_kriptik_PSI_ornek_duzeyi.tsv", sep="\t")
print("\nörnek düzeyi tablo:", S.shape)

# ------------------------------------------------------------ 4. ALS vs kontrol
DERINLIK = 20     # verici üzerinde en az bu kadar okuma
def cliffs(a, b):
    a = np.asarray(a); b = np.asarray(b)
    return (np.sign(a[:, None] - b[None, :]).sum()) / (len(a) * len(b))

BOL = ["Cortex_Frontal", "Cerebellum", "Cortex_Motor", "Cortex_Temporal", "Hippocampus",
       "Spinal_Cord_Cervical", "Spinal_Cord_Lumbar", "Spinal_Cord_Thoracic"]
S["doku_k"] = S["doku"].str.replace(" ", "_")
rows = []
for gen in KRIPTIK:
    ok = S[S[f"{gen}_toplam"] >= DERINLIK]
    for doku, g in ok.groupby("doku_k"):
        als = g[g["grup"].str.strip() == "ALS Spectrum MND"][f"{gen}_PSI"].dropna()
        ctl = g[g["grup"].str.contains("Non-Neurological", na=False)][f"{gen}_PSI"].dropna()
        if len(als) < 10 or len(ctl) < 5: continue
        u = stats.mannwhitneyu(als, ctl, alternative="two-sided")
        rows.append(dict(gen=gen, doku=doku, n_ALS=len(als), n_kontrol=len(ctl),
                         medyan_ALS=als.median(), medyan_kontrol=ctl.median(),
                         cliffs_delta=cliffs(als.values, ctl.values), p=u.pvalue))
V = pd.DataFrame(rows)
if len(V):
    V["q"] = L.bh(V["p"].values)
    V = V.sort_values(["gen", "doku"])
    V.to_csv(f"{OUT}/NYGC_kriptik_PSI_ALS_vs_kontrol.tsv", sep="\t", index=False)
    print("\n=== Kriptik PSI: ALS vs nörolojik olmayan kontrol ===")
    print(V.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

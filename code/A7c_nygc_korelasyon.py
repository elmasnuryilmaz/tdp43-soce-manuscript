#!/usr/bin/env python3
"""A7c — Makalenin başarısız kalan mekanistik testinin doğru biçimi.

Eski test: TRPC1 ~ gen düzeyi STMN2 (dokuz bölgede de POZİTİF çıkmıştı; gen
düzeyi STMN2 yığın dokuda nöronal içeriği ölçtüğü için test bilgilendirici değil).
Yeni test: TRPC1 ~ kriptik STMN2 PSI (aynı gen içinde oran; nöronal içerikten
büyük ölçüde bağımsız).
"""
import os, sys, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
import lib_junc as L

OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
CACHE = "/Users/elmas/Desktop/TEZ/output/TRPC1_kisi_duzeyi_23.08.26"
SYM = "/Users/elmas/Desktop/TEZ/output/ek_analizler_2026-07-31/00_ham_veri_onbellek/encode_nmd/ensg_symbol.tsv"
DERINLIK = 20

S = pd.read_csv(f"{OUT}/NYGC_kriptik_PSI_ornek_duzeyi.tsv", sep="\t", index_col=0, low_memory=False)
print("kriptik PSI tablosu:", S.shape, "| grup değerleri:", S["grup"].value_counts().head(5).to_dict())

# --- gen düzeyi ifade (makalenin kullandığı yerel NYGC matrisi)
cnt = pickle.load(open(f"{CACHE}/_cnt.pkl", "rb"))
sym = pd.read_csv(SYM, sep="\t", index_col=0)
e2s = sym["gene_name"].to_dict(); s2e = {}
for e, s in e2s.items(): s2e.setdefault(s, e)
cnt = cnt.apply(pd.to_numeric, errors="coerce")
cnt = cnt.loc[:, cnt.notna().any(axis=0)].fillna(0)
cnt = cnt.loc[:, cnt.sum(0) > 0]
cpm = np.log2(cnt.divide(cnt.sum(0), axis=1) * 1e6 + 1)
print("gen düzeyi matris:", cpm.shape)

S = S[S["cgnd"].notna()].copy()
ortak = [c for c in cpm.columns if c in set(S["cgnd"])]
print("eşleşen örnek:", len(ortak))
cpm = cpm[ortak]
S = S.set_index("cgnd").loc[[c for c in ortak]]

ILGI = ["TRPC1", "SARAF", "CBARP", "STIMATE", "ORAI1", "ORAI2", "ATP2A2", "ATP2A3",
        "STIM1", "STIM2", "TARDBP", "STMN2", "UNC13A", "SNAP25", "GFAP", "RBFOX3"]
def vec(g):
    e = s2e.get(g)
    return cpm.loc[e] if (e and e in cpm.index) else None
G = pd.DataFrame({g: vec(g) for g in ILGI if vec(g) is not None}).loc[cpm.columns]
X = S.join(G)
X["kriptik_PSI"] = X["STMN2_PSI"]
X["kriptik_pozitif"] = (X["STMN2_kriptik"] > 0).astype(int)
X.to_csv(f"{OUT}/NYGC_kriptik_ve_ifade_birlesik.tsv", sep="\t")

als = X[(X["grup"].astype(str).str.strip() == "ALS Spectrum MND")
        & (X["STMN2_toplam"] >= DERINLIK)]
print(f"\nsaf ALS + yeterli kapsam: {len(als)} örnek")

HEDEF = ["TRPC1", "SARAF", "CBARP", "ATP2A2", "ORAI2", "STIM1", "STIMATE",
         "ORAI1", "ATP2A3", "SNAP25", "GFAP", "STMN2"]
rows = []
for doku, g in als.groupby("doku"):
    if len(g) < 25: continue
    for hedef in HEDEF:
        if hedef not in g.columns: continue
        y = g[hedef].astype(float)
        for vekil, ad in ((g["kriptik_PSI"].astype(float), "kriptik_STMN2_PSI"),
                          (g["STMN2"].astype(float), "gen_duzeyi_STMN2")):
            ok = y.notna() & vekil.notna() & np.isfinite(y) & np.isfinite(vekil)
            if ok.sum() < 25 or vekil[ok].nunique() < 5: continue
            rho, p = stats.spearmanr(vekil[ok], y[ok])
            rows.append(dict(doku=doku, hedef=hedef, vekil=ad, n=int(ok.sum()),
                             rho=round(rho, 3), p=p))
R = pd.DataFrame(rows)
R["q"] = L.bh(R["p"].values)
R.to_csv(f"{OUT}/NYGC_kriptikPSI_korelasyon.tsv", sep="\t", index=False)

print("\n=== TRPC1 / SARAF / CBARP: iki vekil karşılaştırması (ALS örnekleri) ===")
piv = R[R["hedef"].isin(["TRPC1", "SARAF", "CBARP", "STMN2", "SNAP25"])].pivot_table(
    index=["doku", "hedef"], columns="vekil", values="rho")
print(piv.to_string(float_format=lambda x: f"{x:+.3f}"))

print("\n=== Kriptik PSI vekiliyle anlamlı korelasyonlar (q<0.05) ===")
sig = R[(R["vekil"] == "kriptik_STMN2_PSI") & (R["q"] < 0.05)].sort_values("q")
print(sig.to_string(index=False, float_format=lambda x: f"{x:.4g}") if len(sig) else "yok")

# --- kriptik-pozitif vs kriptik-negatif ALS örnekleri
print("\n=== Kriptik-pozitif vs kriptik-negatif ALS örnekleri ===")
rows2 = []
for doku, g in als.groupby("doku"):
    a = g[g["kriptik_pozitif"] == 1]; b = g[g["kriptik_pozitif"] == 0]
    if len(a) < 10 or len(b) < 10: continue
    for hedef in ["TRPC1", "SARAF", "CBARP", "ATP2A2", "ORAI2", "SNAP25"]:
        if hedef not in g.columns: continue
        u = stats.mannwhitneyu(a[hedef].dropna(), b[hedef].dropna(), alternative="two-sided")
        d = (np.sign(a[hedef].dropna().values[:, None] - b[hedef].dropna().values[None, :]).sum()
             / (a[hedef].notna().sum() * b[hedef].notna().sum()))
        rows2.append(dict(doku=doku, hedef=hedef, n_poz=len(a), n_neg=len(b),
                          cliffs_delta=round(d, 3), p=u.pvalue))
R2 = pd.DataFrame(rows2)
if len(R2):
    R2["q"] = L.bh(R2["p"].values)
    R2.to_csv(f"{OUT}/NYGC_kriptik_pozitif_vs_negatif.tsv", sep="\t", index=False)
    print(R2.sort_values("p").head(20).to_string(index=False, float_format=lambda x: f"{x:.4g}"))

#!/usr/bin/env python3
"""
MS'teki TRPC1 azalmasi demiyelinizasyonun mu yansimasi?

TRPC1 oligodendrosit/OPC'lerde ifade edilir; MS'te miyelin kaybi (MBP, PLP1 azalmasi)
tek basina TRPC1 dususu uretebilir. Iki kontrol:
  (1) NAWM alt kumesi — miyelin korunmus dokuda TRPC1 hala dusuk mu?
  (2) MBP/PLP1/GFAP icin duzeltilmis kismi analiz.
"""
import os, gzip
import numpy as np
import pandas as pd
from scipy import stats

O = "/Users/elmas/Desktop/MAKALE/06_MS_ANALIZI"
RAW = os.path.join(O, "ham_veri")
SYM = ("/Users/elmas/Desktop/TEZ/output/ek_analizler_2026-07-31/"
       "00_ham_veri_onbellek/encode_nmd/ensg_symbol.tsv")
sym = pd.read_csv(SYM, sep="\t", index_col=0)["gene_name"].to_dict()

meta_rows = {}
with gzip.open(os.path.join(RAW, "GSE138614_series_matrix.txt.gz"), "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("!Sample_"):
            p = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            meta_rows.setdefault(p[0], []).append(p[1:])
titles = meta_rows["!Sample_title"][0]
ch = meta_rows["!Sample_characteristics_ch1"]
def pick(pref):
    for row in ch:
        if row[0].lower().startswith(pref):
            return [v.split(":", 1)[1].strip() if ":" in v else v for v in row]
diag, lesion, indiv = pick("diagnosis"), pick("lesion type"), pick("individual")
m = pd.DataFrame({"title": titles, "tani": diag, "lezyon": lesion, "denek": indiv})
m["col"] = "G58-" + m.title.str.extract(r"Sample_(\d+)")[0]
m = m.set_index("col")

C = pd.read_csv(os.path.join(RAW, "GSE138614_countMatrix.txt.gz"), sep="\t", index_col=0)
C.index = [sym.get(str(i).split(".")[0], None) for i in C.index]
C = C[[i is not None for i in C.index]]
C = C.groupby(level=0).sum()
cols = [c for c in C.columns if c in m.index]
C, m = C[cols], m.loc[cols]

cpm = C.divide(C.sum(0), axis=1) * 1e6
cpm = cpm[(cpm > 0.4).sum(axis=1) >= 20]
lg = np.log2(cpm + 1)
lg = lg.sub(lg.median(axis=0), axis=1)

def cliffs_test(x, y):
    U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 2.0 * U / (len(x) * len(y)) - 1.0, p

print("=" * 96)
print("1) NAWM ALT KUMESI — miyelin korunmus MS dokusunda TRPC1")
print("=" * 96)
nawm = [c for c in lg.columns if "NAWM" in str(m.loc[c, "lezyon"])]
ktr = [c for c in lg.columns if m.loc[c, "tani"] == "Control"]
print(f"MS NAWM n = {len(nawm)}, kontrol WM n = {len(ktr)}\n")
print(f"{'gen':10s} {'delta':>8s} {'p':>10s}   yorum")
for g in ["MBP", "PLP1", "MOG", "MAG", "GFAP", "TRPC1"]:
    if g not in lg.index:
        continue
    d, p = cliffs_test(lg.loc[g, nawm].astype(float), lg.loc[g, ktr].astype(float))
    yorum = ""
    if g in ("MBP", "PLP1", "MOG", "MAG"):
        yorum = "miyelin korunmus" if abs(d) < 0.3 else "miyelin kaybi var"
    print(f"{g:10s} {d:+8.3f} {p:10.5f}   {yorum}")

print("\n" + "=" * 96)
print("2) MBP + PLP1 + GFAP ICIN DUZELTILMIS TRPC1 (tum ornekler)")
print("=" * 96)
ms = [c for c in lg.columns if m.loc[c, "tani"] != "Control"]
X = pd.DataFrame({
    "TRPC1": lg.loc["TRPC1"],
    "MBP": lg.loc["MBP"], "PLP1": lg.loc["PLP1"], "GFAP": lg.loc["GFAP"],
    "hasta": [1 if m.loc[c, "tani"] != "Control" else 0 for c in lg.columns],
}).astype(float).dropna()

# artik yaklasimi: TRPC1'i miyelin/gliya belirteclerine gore regresle, artigi test et
A = np.column_stack([X.MBP, X.PLP1, X.GFAP, np.ones(len(X))])
beta = np.linalg.lstsq(A, X.TRPC1.values, rcond=None)[0]
resid = X.TRPC1.values - A @ beta
X["TRPC1_artik"] = resid

d0, p0 = cliffs_test(X.loc[X.hasta == 1, "TRPC1"], X.loc[X.hasta == 0, "TRPC1"])
d1, p1 = cliffs_test(X.loc[X.hasta == 1, "TRPC1_artik"], X.loc[X.hasta == 0, "TRPC1_artik"])
print(f"  ham TRPC1        : delta = {d0:+.3f}   p = {p0:.5f}")
print(f"  duzeltilmis TRPC1: delta = {d1:+.3f}   p = {p1:.5f}")

r_mbp = stats.spearmanr(X.TRPC1, X.MBP)
r_plp = stats.spearmanr(X.TRPC1, X.PLP1)
print(f"\n  TRPC1 ~ MBP  : rho = {r_mbp.statistic:+.3f}  (p = {r_mbp.pvalue:.2e})")
print(f"  TRPC1 ~ PLP1 : rho = {r_plp.statistic:+.3f}  (p = {r_plp.pvalue:.2e})")

print("\n" + "=" * 96)
print("3) LEZYON TIPINE GORE TRPC1 (kontrol WM'ye karsi)")
print("=" * 96)
for lz in ["Normal appearing white matter (NAWM)", "Remyelinating (RL)",
           "Active (AL)", "Chronic active (CA)", "Inactive (IL)"]:
    sub = [c for c in lg.columns if str(m.loc[c, "lezyon"]) == lz]
    if len(sub) < 4:
        continue
    dt, pt = cliffs_test(lg.loc["TRPC1", sub].astype(float), lg.loc["TRPC1", ktr].astype(float))
    dm, _ = cliffs_test(lg.loc["MBP", sub].astype(float), lg.loc["MBP", ktr].astype(float))
    print(f"  {lz:44s} n={len(sub):3d}  TRPC1 δ={dt:+.3f} (p={pt:.4f})   MBP δ={dm:+.3f}")

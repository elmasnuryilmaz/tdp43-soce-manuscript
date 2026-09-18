#!/usr/bin/env python3
"""
MS (multipl skleroz) dokusunda TRPC1 yonu.

Yontem, NYGC/AD/PD analizleriyle BIREBIR ayni tutuldu:
ham sayim -> CPM filtresi -> log2CPM -> ornek bazinda medyan-merkezleme ->
Mann-Whitney U -> BH-FDR -> Cliff's delta.

Iki bagimsiz kohort:
  GSE123496 — 5 MS + 5 kontrol donor, 5 beyin bolgesi (gri + ak madde)
  GSE138614 — 10 MS + 5 kontrol donor, ak madde; lezyon tipine gore
"""
import os, gzip
import numpy as np
import pandas as pd
from scipy import stats

O = "/Users/elmas/Desktop/MAKALE/06_MS_ANALIZI"
RAW = os.path.join(O, "ham_veri")
SYM = ("/Users/elmas/Desktop/TEZ/output/ek_analizler_2026-07-31/"
       "00_ham_veri_onbellek/encode_nmd/ensg_symbol.tsv")

GENLER = ["TRPC1", "SARAF", "CBARP", "STIMATE", "STIM1", "STIM2",
          "ORAI1", "ORAI2", "ORAI3", "ATP2A2", "ATP2A3",
          "STMN2", "UNC13A", "TARDBP",
          "GFAP", "AIF1", "SNAP25", "RBFOX3", "MBP", "PLP1"]
CPM_ESIK = 0.40   # AD/PD analizindeki DE ile ayni


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def cliffs(U, n1, n2):
    return 2.0 * U / (n1 * n2) - 1.0


sym = pd.read_csv(SYM, sep="\t", index_col=0)["gene_name"].to_dict()


def prep(counts):
    """CPM filtresi -> log2CPM -> ornek bazinda medyan merkezleme."""
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    counts = counts.loc[:, counts.sum(0) > 0]
    cpm = counts.divide(counts.sum(0), axis=1) * 1e6
    keep = (cpm > CPM_ESIK).sum(axis=1) >= max(3, int(0.2 * cpm.shape[1]))
    cpm = cpm[keep]
    lg = np.log2(cpm + 1)
    return lg.sub(lg.median(axis=0), axis=1)   # ornek bazinda medyan merkezleme


def test(lg, a, b, etiket, genler=GENLER):
    """a = hasta sutunlari, b = kontrol sutunlari."""
    out = []
    for g in genler:
        if g not in lg.index:
            continue
        x = lg.loc[g, a].astype(float).dropna()
        y = lg.loc[g, b].astype(float).dropna()
        if len(x) < 3 or len(y) < 3:
            continue
        U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        out.append(dict(karsilastirma=etiket, gen=g, n_hasta=len(x), n_kontrol=len(y),
                        medyan_fark=round(float(x.median() - y.median()), 4),
                        cliffs_delta=round(cliffs(U, len(x), len(y)), 3), p=p))
    d = pd.DataFrame(out)
    if not d.empty:
        d["q"] = bh(d.p.values)
        d["anlamli"] = np.where(d.q < 0.05, "*", "")
        d["yon"] = np.where(d.cliffs_delta > 0, "artmis", "azalmis")
    return d


sonuc = []

# =====================================================================
# GSE123496
# =====================================================================
print("=" * 104)
print("KOHORT 1 — GSE123496 · 5 MS + 5 kontrol donor · 5 beyin bolgesi")
print("=" * 104)

C1 = pd.read_csv(os.path.join(RAW, "GSE123496_Human_MSNL_counts.csv.gz"), index_col=0)
C1.index = [sym.get(str(i).split(".")[0], None) for i in C1.index]
C1 = C1[[i is not None for i in C1.index]]
C1 = C1.groupby(level=0).sum()
print(f"Matris: {C1.shape[0]} gen x {C1.shape[1]} ornek")

SUF = {"CC": "Corpus callosum", "Cp": "Parietal cortex", "H": "Hippocampus",
       "IC": "Internal capsule", "": "Frontal cortex"}
info = {}
for c in C1.columns:
    grp = "MS" if c.startswith("MS") else "Kontrol"
    suf = "".join(ch for ch in c if ch.isalpha())
    suf = suf.replace("MS", "").replace("NL", "")
    info[c] = (grp, SUF.get(suf, suf or "Frontal cortex"))
inf = pd.DataFrame(info, index=["grup", "bolge"]).T
print(inf.groupby(["bolge", "grup"]).size().unstack(fill_value=0).to_string())

lg1 = prep(C1)
for bol in sorted(inf.bolge.unique()):
    a = [c for c in lg1.columns if inf.loc[c, "bolge"] == bol and inf.loc[c, "grup"] == "MS"]
    b = [c for c in lg1.columns if inf.loc[c, "bolge"] == bol and inf.loc[c, "grup"] == "Kontrol"]
    d = test(lg1, a, b, f"GSE123496 · {bol}")
    if not d.empty:
        sonuc.append(d)

# tum bolgeler havuzlanmis (bolge icinde merkezlenmis)
lg1c = lg1.copy()
for bol in inf.bolge.unique():
    cols = [c for c in lg1c.columns if inf.loc[c, "bolge"] == bol]
    lg1c[cols] = lg1c[cols].sub(lg1c[cols].mean(axis=1), axis=0)
a = [c for c in lg1c.columns if inf.loc[c, "grup"] == "MS"]
b = [c for c in lg1c.columns if inf.loc[c, "grup"] == "Kontrol"]
d = test(lg1c, a, b, "GSE123496 · TUM BOLGELER (bolge icinde merkezlenmis)")
sonuc.append(d)

# =====================================================================
# GSE138614
# =====================================================================
print("\n" + "=" * 104)
print("KOHORT 2 — GSE138614 · 10 MS + 5 kontrol donor · ak madde, lezyon tipine gore")
print("=" * 104)

meta_rows = {}
with gzip.open(os.path.join(RAW, "GSE138614_series_matrix.txt.gz"), "rt", errors="replace") as fh:
    for line in fh:
        if not line.startswith("!Sample_"):
            continue
        p = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
        meta_rows.setdefault(p[0], []).append(p[1:])

titles = meta_rows["!Sample_title"][0]
ch = meta_rows["!Sample_characteristics_ch1"]
def pick(pref):
    for row in ch:
        if row[0].lower().startswith(pref):
            return [v.split(":", 1)[1].strip() if ":" in v else v for v in row]
    return None
diag = pick("diagnosis"); lesion = pick("lesion type"); indiv = pick("individual")

m2 = pd.DataFrame({"title": titles, "tani": diag, "lezyon": lesion, "denek": indiv})
m2["sample_no"] = m2.title.str.extract(r"Sample_(\d+)")[0]
m2["col"] = "G58-" + m2.sample_no
m2 = m2.set_index("col")

C2 = pd.read_csv(os.path.join(RAW, "GSE138614_countMatrix.txt.gz"), sep="\t", index_col=0)
C2.index = [sym.get(str(i).split(".")[0], None) for i in C2.index]
C2 = C2[[i is not None for i in C2.index]]
C2 = C2.groupby(level=0).sum()
ortak = [c for c in C2.columns if c in m2.index]
C2 = C2[ortak]; m2 = m2.loc[ortak]
print(f"Matris: {C2.shape[0]} gen x {C2.shape[1]} ornek")
print(m2.groupby(["tani", "lezyon"]).size().to_string())
print(f"\nDenek sayisi — MS: {m2[m2.tani!='Control'].denek.nunique()}, "
      f"Kontrol: {m2[m2.tani=='Control'].denek.nunique()}")

lg2 = prep(C2)
kontrol = [c for c in lg2.columns if m2.loc[c, "tani"] == "Control"]

# (a) MS NAWM vs kontrol ak madde  — en karsilastirilabilir olan
nawm = [c for c in lg2.columns if "NAWM" in str(m2.loc[c, "lezyon"])]
d = test(lg2, nawm, kontrol, "GSE138614 · MS NAWM vs kontrol WM")
sonuc.append(d)

# (b) tum MS lezyonlari vs kontrol
lez = [c for c in lg2.columns if m2.loc[c, "tani"] != "Control" and "NAWM" not in str(m2.loc[c, "lezyon"])]
d = test(lg2, lez, kontrol, "GSE138614 · MS lezyonlari vs kontrol WM")
sonuc.append(d)

# (c) DENEK duzeyinde (psodoreplikasyon kontrolu)
dn = lg2.T.copy(); dn["denek"] = m2.loc[dn.index, "denek"].values
dn["tani"] = m2.loc[dn.index, "tani"].values
agg = dn.groupby(["tani", "denek"]).mean(numeric_only=True).reset_index()
lgd = agg.drop(columns=["tani", "denek"]).T
lgd.columns = agg.denek.values
ms_d = list(agg[agg.tani != "Control"].denek)
kt_d = list(agg[agg.tani == "Control"].denek)
d = test(lgd, ms_d, kt_d, "GSE138614 · DENEK duzeyi (10 MS vs 5 kontrol)")
sonuc.append(d)

# =====================================================================
R = pd.concat(sonuc, ignore_index=True)
R.to_csv(os.path.join(O, "MS_TRPC1_sonuclar.tsv"), sep="\t", index=False)

pd.set_option("display.width", 210)
print("\n" + "=" * 104)
print("TRPC1 — TUM KARSILASTIRMALAR")
print("=" * 104)
t = R[R.gen == "TRPC1"]
print(t[["karsilastirma", "n_hasta", "n_kontrol", "medyan_fark",
         "cliffs_delta", "p", "q", "yon", "anlamli"]].to_string(index=False))

print("\n" + "=" * 104)
print("ODAK GENLER — her karsilastirmada Cliff's delta (yildiz = q < 0,05)")
print("=" * 104)
piv = R.pivot_table(index="gen", columns="karsilastirma", values="cliffs_delta")
qiv = R.pivot_table(index="gen", columns="karsilastirma", values="q")
show = piv.copy().astype(object)
for i in piv.index:
    for c in piv.columns:
        v, q = piv.loc[i, c], qiv.loc[i, c]
        show.loc[i, c] = "" if pd.isna(v) else f"{v:+.3f}" + ("*" if (pd.notna(q) and q < 0.05) else "")
ORDER = ["TRPC1", "SARAF", "CBARP", "STIMATE", "STIM1", "STIM2", "ORAI1", "ORAI2",
         "ORAI3", "ATP2A2", "ATP2A3", "STMN2", "UNC13A", "TARDBP",
         "GFAP", "AIF1", "SNAP25", "RBFOX3", "MBP", "PLP1"]
show = show.reindex([g for g in ORDER if g in show.index])
print(show.to_string())
print("\nSonuc ->", os.path.join(O, "MS_TRPC1_sonuclar.tsv"))

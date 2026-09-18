#!/usr/bin/env python3
"""
TRPC1/SARAF/CBARP artisi, TDP-43 kayip SIDDETIYLE olcekleniyor mu?

"ALS'te artmis" tanimlayici bir gozlemdir. Asil mekanistik iddia, degisimin
TDP-43 islev kaybiyla orantili olmasidir. NYGC kohortunda, YALNIZCA ALS
orneklerinde, her ornek icin bir TDP-43 kayip vekili hesaplanip ilgi genleriyle
korelasyonu olculur.

TDP-43 kayip vekili: STMN2 ve UNC13A gen duzeyi ifadesi. Her ikisi de TDP-43
kaybinda kriptik olaylar nedeniyle AZALIR (Melamed 2019, Klim 2019, Brown 2022;
NYGC'de gen duzeyinde Prudencio 2020 ayni sekilde kullanmistir).
Dolayisiyla: TRPC1 ~ STMN2 arasinda NEGATIF korelasyon beklenir.
"""
import os, pickle
import numpy as np
import pandas as pd
from scipy import stats

BASE = "/Users/elmas/Desktop/TEZ/output"
CACHE = os.path.join(BASE, "TRPC1_kisi_duzeyi_23.08.26")
SYM = os.path.join(BASE, "ek_analizler_2026-07-31/00_ham_veri_onbellek/encode_nmd/ensg_symbol.tsv")
OUT = "/Users/elmas/Desktop/MAKALE"
TAB = os.path.join(OUT, "03_TABLOLAR")

cnt = pickle.load(open(os.path.join(CACHE, "_cnt.pkl"), "rb"))
meta = pickle.load(open(os.path.join(CACHE, "_meta.pkl"), "rb"))
sym = pd.read_csv(SYM, sep="\t", index_col=0)
e2s = sym["gene_name"].to_dict()
s2e = {}
for e, s in e2s.items():
    s2e.setdefault(s, e)

ILGI = ["TRPC1", "SARAF", "CBARP", "STIMATE", "ORAI1", "ORAI2",
        "ATP2A2", "ATP2A3", "STIM1", "TARDBP"]
KAYIP = ["STMN2", "UNC13A"]
KONTROL_GEN = ["GFAP", "SNAP25", "RBFOX3"]   # hucre bilesimi vekilleri

# --- sayisal olmayan sutunlari at, CPM + log2
cnt = cnt.apply(pd.to_numeric, errors="coerce")
cnt = cnt.loc[:, cnt.notna().any(axis=0)]
cnt = cnt.fillna(0)
cnt = cnt.loc[:, cnt.sum(0) > 0]
cpm = np.log2(cnt.divide(cnt.sum(0), axis=1) * 1e6 + 1)

meta = meta.set_index("ornek_id")
ortak = [c for c in cpm.columns if c in meta.index]
cpm = cpm[ortak]
meta = meta.loc[ortak]

als = meta[meta.grup.astype(str).str.strip() == "ALS Spectrum MND"]
print(f"Toplam ornek: {len(meta)} | saf ALS ornegi: {len(als)}")

def vec(g):
    e = s2e.get(g)
    return cpm.loc[e] if (e and e in cpm.index) else None

BOLGELER = ["Cortex Frontal", "Cerebellum", "Cortex Motor Lateral",
            "Cortex Motor Medial", "Cortex Temporal", "Hippocampus",
            "Spinal Cord Cervical", "Spinal Cord Lumbar", "Spinal Cord Thoracic"]

rows = []
for bol in BOLGELER:
    idx = als[als.doku == bol].index
    if len(idx) < 20:
        continue
    for gen in ILGI:
        gv = vec(gen)
        if gv is None:
            continue
        for kg in KAYIP:
            kv = vec(kg)
            if kv is None:
                continue
            x = gv[idx].astype(float)
            y = kv[idx].astype(float)
            ok = x.notna() & y.notna() & np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 20:
                continue
            rho, p = stats.spearmanr(x[ok], y[ok])
            rows.append(dict(bolge=bol, gen=gen, kayip_vekili=kg, n=int(ok.sum()),
                             rho=round(rho, 3), p=p))

R = pd.DataFrame(rows)
R["q"] = np.nan
ok = R.p.notna()
pv = R.loc[ok, "p"].values
order = np.argsort(pv); ranked = pv[order]; n = len(pv)
q = np.minimum.accumulate((ranked * n / (np.arange(n) + 1))[::-1])[::-1]
qq = np.empty(n); qq[order] = np.minimum(q, 1.0)
R.loc[ok, "q"] = qq
R["anlamli"] = np.where(R.q < 0.05, "*", "")
R.to_csv(os.path.join(TAB, "NYGC_TDP43_kayip_korelasyonu.tsv"), sep="\t", index=False)

pd.set_option("display.width", 200)
print("\n" + "=" * 100)
print("ALS ORNEKLERINDE: ilgi geni ~ TDP-43 kayip vekili (Spearman rho)")
print("Beklenti: TDP-43 kaybi arttikca STMN2/UNC13A DUSER.")
print("TRPC1 TDP-43 kaybiyla artiyorsa rho NEGATIF olmalidir.")
print("=" * 100)

for gen in ILGI:
    sub = R[R.gen == gen]
    if sub.empty:
        continue
    stm = sub[sub.kayip_vekili == "STMN2"]
    neg = (stm.rho < 0).sum(); tot = len(stm)
    sig = (stm.q < 0.05).sum()
    print(f"\n### {gen}  —  STMN2 ile: {neg}/{tot} bolgede negatif, {sig} anlamli")
    print(sub.pivot_table(index="bolge", columns="kayip_vekili",
                          values="rho").to_string())

print("\n" + "=" * 100)
print("ANLAMLI KORELASYONLAR (q < 0,05)")
print("=" * 100)
S = R[R.q < 0.05].sort_values(["gen", "bolge"])
print(S[["bolge", "gen", "kayip_vekili", "n", "rho", "p", "q"]].to_string(index=False)
      if not S.empty else "  (yok)")

# ---- hucre bilesimi kontrolu: kismi korelasyon
print("\n" + "=" * 100)
print("HUCRE BILESIMI KONTROLU — GFAP ve SNAP25 icin duzeltilmis kismi korelasyon")
print("=" * 100)

def partial_spearman(x, y, covs):
    """Sirali donusum sonrasi kovaryatlara gore artiklarin korelasyonu."""
    from numpy.linalg import lstsq
    xr = stats.rankdata(x); yr = stats.rankdata(y)
    C = np.column_stack([stats.rankdata(c) for c in covs] + [np.ones(len(xr))])
    bx = lstsq(C, xr, rcond=None)[0]; by = lstsq(C, yr, rcond=None)[0]
    rx = xr - C @ bx; ry = yr - C @ by
    r, p = stats.pearsonr(rx, ry)
    return r, p

prows = []
for bol in BOLGELER:
    idx = als[als.doku == bol].index
    if len(idx) < 20:
        continue
    gf, sn = vec("GFAP"), vec("SNAP25")
    if gf is None or sn is None:
        continue
    for gen in ["TRPC1", "SARAF", "CBARP"]:
        gv, kv = vec(gen), vec("STMN2")
        if gv is None or kv is None:
            continue
        d = pd.DataFrame({"g": gv[idx], "k": kv[idx],
                          "gf": gf[idx], "sn": sn[idx]}).astype(float).dropna()
        if len(d) < 20:
            continue
        r0, p0 = stats.spearmanr(d.g, d.k)
        r1, p1 = partial_spearman(d.g.values, d.k.values,
                                  [d.gf.values, d.sn.values])
        prows.append(dict(bolge=bol, gen=gen, n=len(d),
                          rho_ham=round(r0, 3), p_ham=round(p0, 5),
                          rho_kismi=round(r1, 3), p_kismi=round(p1, 5)))
P = pd.DataFrame(prows)
P.to_csv(os.path.join(TAB, "NYGC_kismi_korelasyon_GFAP_SNAP25.tsv"), sep="\t", index=False)
for gen in ["TRPC1", "SARAF", "CBARP"]:
    print(f"\n### {gen} ~ STMN2 (GFAP + SNAP25 icin duzeltilmis)")
    print(P[P.gen == gen].drop(columns=["gen"]).to_string(index=False))

print("\nTablolar ->", TAB)

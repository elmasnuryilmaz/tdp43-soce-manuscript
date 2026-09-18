#!/usr/bin/env python3
"""
Oncelik 4.1 — TEZIN GERCEK gen panelleriyle eslenmis permutasyon null'u.
Paneller: output/calcium_panel_overlap_analysis_2026-04-25/01_gene_panels/
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy import stats
from importlib.machinery import SourceFileLoader

core = SourceFileLoader("core", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "01_rmats_core.py")).load_module()
TAB = os.path.join(core.OUT, "03_TABLOLAR")
PANEL = "/Users/elmas/Desktop/TEZ/output/calcium_panel_overlap_analysis_2026-04-25/01_gene_panels"
rng = np.random.default_rng(2026)

panels = {}
for fp in sorted(glob.glob(os.path.join(PANEL, "Tier*.csv"))):
    name = os.path.basename(fp).replace(".csv", "")
    d = pd.read_csv(fp)
    col = "gene" if "gene" in d.columns else d.columns[0]
    panels[name] = set(d[col].astype(str).str.upper().str.strip())
    print(f"{name}: {len(panels[name])} gen")

results = []
for DS in core.DATASETS:
    frames = [core.load_event(DS, ev, "JC") for ev in core.EVENTS]
    E = pd.concat([f for f in frames if f is not None], ignore_index=True)
    E["gene_up"] = E["geneSymbol"].str.upper()
    E["sig"] = (E["FDR"] < 0.05) & (E["IncLevelDifference"].abs() >= 0.10)
    g = (E.groupby("gene_up")
           .agg(n_olay=("ID", "size"), toplam_okuma=("total_reads", "sum"),
                anlamli=("sig", "any"))
           .reset_index())
    g = g[g.n_olay > 0]
    g["logev"] = np.log10(g.n_olay + 1)
    g["logreads"] = np.log10(g.toplam_okuma + 1)

    for pname, pset in panels.items():
        g["in_set"] = g.gene_up.isin(pset)
        n = int(g.in_set.sum()); k = int(g.loc[g.in_set, "anlamli"].sum())
        if n < 10:
            continue
        obs = k / n
        N = len(g); K = int(g.anlamli.sum())
        p_hyper = stats.hypergeom.sf(k - 1, N, K, n)

        cand = g[~g.in_set]
        C = cand[["logev", "logreads"]].values
        S = g[g.in_set][["logev", "logreads"]].values
        cs = cand["anlamli"].values.astype(float)
        pools = np.array([np.argsort(np.sqrt(((C - S[i]) ** 2).sum(1)))[:25]
                          for i in range(S.shape[0])])
        NP = 5000
        null = np.empty(NP)
        for b in range(NP):
            pick = pools[np.arange(pools.shape[0]),
                         rng.integers(0, pools.shape[1], pools.shape[0])]
            null[b] = cs[pick].mean()
        p_perm = (1 + (null >= obs).sum()) / (1 + NP)

        results.append(dict(
            veri_seti=DS, panel=pname, n_test_edilebilir=n, anlamli=k,
            gozlenen_oran=round(100 * obs, 1),
            arka_plan_ham=round(100 * cand.anlamli.mean(), 1),
            eslenmis_null=round(100 * null.mean(), 1),
            null_GA=f"{100*np.percentile(null,2.5):.1f}-{100*np.percentile(null,97.5):.1f}",
            p_hipergeometrik=round(p_hyper, 4),
            p_eslenmis_permutasyon=round(p_perm, 4),
            karar=("ZENGINLESME VAR" if p_perm < 0.05 else "zenginlesme yok"),
        ))

R = pd.DataFrame(results)
R.to_csv(os.path.join(TAB, "zenginlesme_gercek_paneller_permutasyon.tsv"),
         sep="\t", index=False)

pd.set_option("display.width", 240)
print("\n" + "=" * 120)
print("TEZIN GERCEK PANELLERIYLE ZENGINLESME — eslenmis permutasyon null'u")
print("=" * 120)
for DS in core.DATASETS:
    sub = R[R.veri_seti == DS]
    if sub.empty:
        continue
    print(f"\n### {DS}")
    print(sub.drop(columns=["veri_seti"]).to_string(index=False))

print("\n" + "=" * 120)
print("BIRINCIL MODEL (SH-SY5Y) OZETI")
print("=" * 120)
print(R[R.veri_seti == "GSE296712_SHSY5Y"].drop(columns=["veri_seti"]).to_string(index=False))
print("\nTablo ->", TAB)

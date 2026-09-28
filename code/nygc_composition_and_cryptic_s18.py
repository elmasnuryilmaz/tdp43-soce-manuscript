#!/usr/bin/env python3
"""Supplementary Tables S18, S18b and S18c: what the NYGC TRPC1 differences follow.

Three questions about the tissue comparison of Section 3.9, each answered on exactly the samples
of Table 5 (the GSE153960 count-matrix columns with metadata, single-label groups):

S18   Does cell composition explain the TRPC1 differences? Within each region TRPC1 is regressed
      on a neuronal marker (SNAP25 or RBFOX3), alone or with GFAP, by ordinary least squares with
      an intercept across all samples of the comparison, as for multiple sclerosis; the residuals
      are compared with Mann-Whitney U and Cliff's delta, Benjamini-Hochberg within group and model.
      The normalisation is that of Table 5: genes with CPM > 1 in at least max(10, 20% of samples),
      log2(CPM + 1), per-sample median centring.
S18b  How common is the cryptic STMN2 junction, the indicator of TDP-43 loss of function, in each
      group and region? The definition is that of Supplementary Table S8 (recount3 junctions, at
      least 20 reads at the shared exon-1 donor); a sample sequenced as more than one library has
      its reads pooled, so each sample counts once.
S18c  Within the comparison group, do the candidate transcripts follow the cryptic junction once
      neuronal content is held constant (partial Spearman correlation on SNAP25)? And how closely
      does TRPC1 follow SNAP25 in every group and region?

Run from 09_YAYIN_PAKETI: /usr/bin/python3 code/nygc_composition_and_cryptic_s18.py
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
BASE = Path("/Users/elmas/Desktop/TEZ/output/ek_analizler_2026-07-31")
JUNCTIONS = Path("/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar/NYGC_kriptik_ve_ifade_birlesik.tsv")
CONTROL, ALS, OND = "Non-Neurological Control", "ALS Spectrum MND", "Other Neurological Disorders"
LABEL = {ALS: "ALS", OND: "Other neurological disorders", CONTROL: "Non-neurological control"}
BRAIN = ["Cerebellum", "Cortex Frontal", "Cortex Motor Lateral", "Cortex Motor Medial",
         "Cortex Occipital", "Cortex Temporal", "Hippocampus"]
CORD = ["Spinal Cord Cervical", "Spinal Cord Lumbar", "Spinal Cord Thoracic"]
OND_REGIONS = ["Cerebellum", "Cortex Frontal", "Cortex Temporal"]
MODELS = {"none": [], "SNAP25": ["SNAP25"], "RBFOX3": ["RBFOX3"],
          "SNAP25 + GFAP": ["SNAP25", "GFAP"], "RBFOX3 + GFAP": ["RBFOX3", "GFAP"]}


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def cliff(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return float(np.sign(a[:, None] - b[None, :]).sum() / (len(a) * len(b)))


def partial_spearman(x, y, z):
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r, _ = stats.pearsonr(ex, ey)
    n = len(x); t = r * np.sqrt((n - 3) / (1 - r * r))
    return r, 2 * stats.t.sf(abs(t), n - 3)


# ------------------------------------------------------------------ the samples of Table 5
meta = pd.read_csv(BASE / "06_WGCNA_ALS_kohort/T60_NYGC_ornek_ustverisi.csv", low_memory=False)
meta = meta.dropna(subset=["ornek_id"]).drop_duplicates("ornek_id").set_index("ornek_id")
meta["grup"] = meta["grup"].astype(str).str.strip().str.replace("DIsorders", "Disorders")
sym = pd.read_csv(BASE / "00_ham_veri_onbellek/encode_nmd/ensg_symbol.tsv", sep="\t", index_col=0)
cnt = pd.read_csv(BASE / "00_ham_veri_onbellek/gse153960/GSE153960_counts.txt.gz", sep="\t", low_memory=False)
ens = "EnsemblID" if "EnsemblID" in cnt.columns else cnt.columns[0]
cnt = cnt.set_index(cnt[ens].astype(str)).drop(columns=[ens]); cnt.index = cnt.index.str.split(".").str[0]
cnt = cnt[[c for c in cnt.columns if str(c).startswith("CGND")]].apply(pd.to_numeric, errors="coerce").dropna(how="all")
cnt = cnt[cnt.index.isin(sym.index)]; cnt.index = sym.reindex(cnt.index)["gene_name"].values
cnt = cnt[~pd.isna(cnt.index)].groupby(level=0).sum()
meta = meta[meta.index.isin(cnt.columns)]
print(f"samples with metadata in the count matrix: {len(meta)}")


def comparison(region, group):
    s = meta[meta.doku == region]
    cols = [c for c in cnt.columns if c in s.index]
    case = [c for c in cols if s.loc[c, "grup"] == group]
    ctrl = [c for c in cols if s.loc[c, "grup"] == CONTROL]
    if len(case) < 10 or len(ctrl) < 10:
        return None
    sub = cnt[case + ctrl]; cpm = sub / sub.sum(axis=0) * 1e6
    sub = sub[(cpm > 1).sum(axis=1) >= max(10, int(0.2 * len(sub.columns)))]
    lg = np.log2(sub / sub.sum(axis=0) * 1e6 + 1)
    return lg.sub(lg.median(axis=0), axis=1), case, ctrl


def residual(y, markers):
    if not markers:
        return y
    Z = np.column_stack([np.ones(len(y))] + markers)
    beta, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return y - Z @ beta


# ------------------------------------------------------------------ S18
rows, coupling = [], []
for region in BRAIN:
    for group in (ALS, OND):
        if group == OND and region not in OND_REGIONS:
            continue
        r = comparison(region, group)
        if r is None:
            continue
        lg, case, ctrl = r
        y = lg.loc["TRPC1"].values
        for model, markers in MODELS.items():
            v = pd.Series(residual(y, [lg.loc[m].values for m in markers]), index=lg.columns)
            rows.append(dict(group=LABEL[group], region=region, variable="TRPC1", adjustment=model,
                             n_case=len(case), n_control=len(ctrl),
                             cliffs_delta=cliff(v[case], v[ctrl]),
                             p_value=stats.mannwhitneyu(v[case], v[ctrl], alternative="two-sided").pvalue))
        for marker in ("SNAP25", "RBFOX3", "GFAP"):
            rows.append(dict(group=LABEL[group], region=region, variable=marker, adjustment="none (marker itself)",
                             n_case=len(case), n_control=len(ctrl),
                             cliffs_delta=cliff(lg.loc[marker, case], lg.loc[marker, ctrl]),
                             p_value=stats.mannwhitneyu(lg.loc[marker, case], lg.loc[marker, ctrl],
                                                        alternative="two-sided").pvalue))
        for who, cols in ((group, case), (CONTROL, ctrl)):
            rho, p = stats.spearmanr(lg.loc["TRPC1", cols], lg.loc["SNAP25", cols])
            coupling.append(dict(group=LABEL[who], region=region, normalised_with=LABEL[group] + " comparison",
                                 n=len(cols), rho=rho, p=p))
S18 = pd.DataFrame(rows)
for _, g in S18.groupby(["group", "variable", "adjustment"]):
    S18.loc[g.index, "q_value"] = bh(g.p_value.values)

# ------------------------------------------------------------------ S18b
J = pd.read_csv(JUNCTIONS, sep="\t", low_memory=False)
J = J[J.cgnd.isin(meta.index)]
per_sample = J.groupby("cgnd").agg(cryptic=("STMN2_kriptik", "sum"), total=("STMN2_toplam", "sum"),
                                   donor=("denek", "first"),
                                   **{g: (g, "first") for g in ("TRPC1", "SARAF", "CBARP", "SNAP25")})
S = meta[["grup", "doku"]].join(per_sample, how="inner")
S = S[S.total >= 20].copy()
S["PSI"] = S.cryptic / S.total
rows = []
for region in BRAIN + CORD:
    s = S[S.doku == region]
    ctrl = s[s.grup == CONTROL]
    for group in (ALS, OND, CONTROL):
        g = s[s.grup == group]
        if group == OND and region not in OND_REGIONS:
            continue
        if len(g) == 0:
            continue
        row = dict(group=LABEL[group], region=region, n_samples=len(g), n_donors=g.donor.nunique(),
                   n_detected=int((g.cryptic > 0).sum()), n_PSI_above_0_01=int((g.PSI > 0.01).sum()),
                   fraction_PSI_above_0_01=float((g.PSI > 0.01).mean()), mean_PSI=float(g.PSI.mean()),
                   cliffs_delta_vs_control=np.nan, p_value=np.nan)
        if group == OND:
            row["cliffs_delta_vs_control"] = cliff(g.PSI, ctrl.PSI)
            row["p_value"] = stats.mannwhitneyu(g.PSI, ctrl.PSI, alternative="two-sided").pvalue
        rows.append(row)
S18b = pd.DataFrame(rows)
m = S18b.group == LABEL[OND]
S18b.loc[m, "q_value"] = bh(S18b.loc[m, "p_value"].values)

# ------------------------------------------------------------------ S18c
rows = []
for region in OND_REGIONS:
    g = S[(S.grup == OND) & (S.doku == region)]
    if len(g) < 25 or g.PSI.nunique() < 5:        # the S9 criteria
        continue
    for target in ("TRPC1", "SARAF", "CBARP", "SNAP25"):
        rho, p = stats.spearmanr(g.PSI, g[target])
        pr, pp = (np.nan, np.nan) if target == "SNAP25" else partial_spearman(g.PSI.values, g[target].values,
                                                                              g.SNAP25.values)
        rows.append(dict(family="cryptic STMN2 PSI within the comparison group", group=LABEL[OND],
                         region=region, x="cryptic STMN2 PSI", y=target, n=len(g), spearman_rho=rho,
                         p_value=p, partial_rho_given_SNAP25=pr, partial_p_value=pp))
C = pd.DataFrame(rows); C["q_value"] = bh(C.p_value.values)
K = pd.DataFrame(coupling).drop_duplicates(["group", "region"])
K = K.rename(columns={"rho": "spearman_rho", "p": "p_value"})
K["family"] = "TRPC1 with the neuronal marker SNAP25"; K["x"] = "SNAP25"; K["y"] = "TRPC1"
K["q_value"] = bh(K.p_value.values)
S18c = pd.concat([C, K[["family", "group", "region", "x", "y", "n", "spearman_rho", "p_value", "q_value",
                        "normalised_with"]]], ignore_index=True)

out = ROOT / "supplementary"
S18.to_csv(out / "S18_TRPC1_cell_composition_adjustment.csv", index=False)
S18b.to_csv(out / "S18b_cryptic_STMN2_by_group_and_region.csv", index=False)
S18c.to_csv(out / "S18c_cryptic_STMN2_within_comparison_group.csv", index=False)
print("written: S18, S18b, S18c")

# ------------------------------------------------------------------ the numbers quoted in the text
t = S18[(S18.variable == "TRPC1")]
sig6 = ["Cerebellum", "Cortex Frontal", "Cortex Motor Lateral", "Cortex Motor Medial", "Cortex Temporal", "Hippocampus"]
a = t[(t.group == "ALS") & t.region.isin(sig6) & (t.adjustment != "none")]
print(f"ALS adjusted delta over six regions and four models: {a.cliffs_delta.min():+.3f} to {a.cliffs_delta.max():+.3f}")
robust = [r for r in sig6 if (a[a.region == r].q_value < 0.05).all()]
print("significant under every model:", robust)
print("significant after RBFOX3 alone:", list(a[(a.adjustment == "RBFOX3") & (a.q_value < 0.05)].region))
print(f"TRPC1-SNAP25 rho range: {K.spearman_rho.min():.2f} to {K.spearman_rho.max():.2f}")

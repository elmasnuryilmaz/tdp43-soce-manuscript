#!/usr/bin/env python3
"""MS: donor-level re-analysis of the TRPC1 result (v3 correction).

run_ms_confound.py treated several samples from the same donor as independent.
Here every test is repeated with the donor as the unit of inference.
Output: 03_TABLOLAR/v3/MS_donor_level.csv
"""
import os, gzip
import numpy as np
import pandas as pd
from scipy import stats

O = "/Users/elmas/Desktop/MAKALE/06_MS_ANALIZI"
RAW = os.path.join(O, "ham_veri")
OUT = "/Users/elmas/Desktop/MAKALE/03_TABLOLAR/v3"
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
# same low-expression filter as run_ms.py, so the two tables share the per-sample
# median used for centring and the donor-level delta is directly comparable
cpm = cpm[(cpm > 0.4).sum(axis=1) >= max(3, int(0.2 * cpm.shape[1]))]
lg = np.log2(cpm + 1)
lg = lg.sub(lg.median(axis=0), axis=1)

def cliffs_test(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 2.0 * U / (len(x) * len(y)) - 1.0, p

rows = []
def add(label, gene, x, y, unit, nx, ny):
    d, p = cliffs_test(x, y)
    rows.append(dict(comparison=label, gene=gene, unit=unit, n_case=nx, n_control=ny,
                     cliffs_delta=round(d, 3), p=float(f"{p:.5g}")))
    print(f"{label:46s} {gene:7s} {unit:7s} n={nx:3d}/{ny:3d}  delta={d:+.3f}  p={p:.5f}")

ktr_cols = [c for c in lg.columns if m.loc[c, "tani"] == "Control"]
nawm_cols = [c for c in lg.columns if "NAWM" in str(m.loc[c, "lezyon"])]
ms_cols = [c for c in lg.columns if m.loc[c, "tani"] != "Control"]

print("=== 1. NAWM vs control white matter ===")
for g in ["MBP", "PLP1", "MOG", "MAG", "GFAP", "TRPC1"]:
    if g not in lg.index:
        continue
    # sample level (as published)
    add("NAWM vs control WM", g, lg.loc[g, nawm_cols], lg.loc[g, ktr_cols], "sample",
        len(nawm_cols), len(ktr_cols))
    # donor level
    dn_case = lg.loc[g, nawm_cols].groupby(m.loc[nawm_cols, "denek"].values).mean()
    dn_ctrl = lg.loc[g, ktr_cols].groupby(m.loc[ktr_cols, "denek"].values).mean()
    add("NAWM vs control WM", g, dn_case, dn_ctrl, "donor", len(dn_case), len(dn_ctrl))

print("\n=== 2. myelin/glia-adjusted TRPC1 (residuals) ===")
X = pd.DataFrame({"TRPC1": lg.loc["TRPC1"], "MBP": lg.loc["MBP"],
                  "PLP1": lg.loc["PLP1"], "GFAP": lg.loc["GFAP"]}).astype(float).dropna()
A = np.column_stack([X.MBP, X.PLP1, X.GFAP, np.ones(len(X))])
beta = np.linalg.lstsq(A, X.TRPC1.values, rcond=None)[0]
X["resid"] = X.TRPC1.values - A @ beta
X["case"] = [1 if m.loc[c, "tani"] != "Control" else 0 for c in X.index]
X["denek"] = [m.loc[c, "denek"] for c in X.index]
for col, lab in [("TRPC1", "raw TRPC1"), ("resid", "myelin+glia-adjusted TRPC1")]:
    add(f"MS vs control, {lab}", "TRPC1", X.loc[X.case == 1, col], X.loc[X.case == 0, col],
        "sample", int((X.case == 1).sum()), int((X.case == 0).sum()))
    dm = X.groupby(["denek", "case"])[col].mean().reset_index()
    add(f"MS vs control, {lab}", "TRPC1", dm.loc[dm.case == 1, col], dm.loc[dm.case == 0, col],
        "donor", int((dm.case == 1).sum()), int((dm.case == 0).sum()))

print("\n=== 3. lesion types, donor level ===")
for lz in ["Normal appearing white matter (NAWM)", "Remyelinating (RL)",
           "Active (AL)", "Chronic active (CA)", "Inactive (IL)"]:
    sub = [c for c in lg.columns if str(m.loc[c, "lezyon"]) == lz]
    if len(sub) < 4:
        continue
    dn_case = lg.loc["TRPC1", sub].groupby(m.loc[sub, "denek"].values).mean()
    dn_ctrl = lg.loc["TRPC1", ktr_cols].groupby(m.loc[ktr_cols, "denek"].values).mean()
    add(f"{lz} vs control WM", "TRPC1", dn_case, dn_ctrl, "donor", len(dn_case), len(dn_ctrl))

os.makedirs(OUT, exist_ok=True)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "MS_donor_level.csv"), index=False)
print("\nwritten:", os.path.join(OUT, "MS_donor_level.csv"))

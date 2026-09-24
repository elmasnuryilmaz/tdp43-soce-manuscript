#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v3 — S1: laboratory source data workbook for the submission package.

Sheets
  README            what each sheet contains and how the values were produced
  TARDBP_qPCR       raw Ct, dCt, ddCt and 2^-ddCt for the three groups (n = 4)
  Target_qPCR_Ct    raw Ct for the four SOCE-associated targets (n = 4 per group)
  Target_qPCR_rel   relative expression per replicate (the values plotted in Figure 1B)
  Fura2             per-replicate ER release and SOCE amplitudes (n = 3)
  WST1              per-well metabolic signal at 48 h (n = 4)
  Summary_stats     group means, SEM and tests used for Figures 1 and 2
"""
import os
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats

TEZ = "/Users/elmas/Desktop/TEZ"
# the comparison group of the target-gene RT-qPCR, Fura-2 and WST-1 experiments
NT = "Non-targeting shRNA control"
# the four targets reported in the manuscript; the lab workbook also holds apoptosis markers,
# which are not part of this study and are not published
TARGETS = ["TRPC1", "STIM1", "ORAI1", "ATP2A3"]
PZ = os.path.join(TEZ, "10_GRAPHPAD_PRISM_DOSYALARI", "01_tez_sekil_kaynaklari")
SUP = "/Users/elmas/Desktop/MAKALE/09_YAYIN_PAKETI/supplementary"
os.makedirs(SUP, exist_ok=True)


def pzfx(path):
    out = {}
    for tab in ET.parse(path).getroot().iter("Table"):
        out[tab.findtext("Title")] = {c.findtext("Title"): [float(d.text) for d in c.iter("d")]
                                      for c in tab.findall("YColumn")}
    return out


def sem(v):
    return float(np.std(v, ddof=1) / np.sqrt(len(v)))


qp = pzfx(os.path.join(PZ, "sekil_4.17_4.18_qPCR.pzfx"))
fu = pzfx(os.path.join(PZ, "sekil_4.20_fura2.pzfx"))
ws = pzfx(os.path.join(PZ, "sekil_4.21_wst1.pzfx"))

# ------------------------------------------------------------------ TARDBP Ct
td = pd.read_excel(os.path.join(TEZ, "SHSY5Y_TDP43_qPCR_Ct_Data.xlsx"),
                   sheet_name="TARDBP_Knockdown_Ct", skiprows=4, nrows=12,
                   names=["group", "replicate", "date", "TARDBP_Ct", "GAPDH_Ct",
                          "dCt", "ddCt", "rel_expression"])
td["group"] = td.group.map({"Kontrol": "Untransduced control",
                            "Hedef dışı kontrol": "Non-targeting shRNA control",
                            "shTDP-43": "shTDP-43"})

# ------------------------------------------------------------- target gene Ct
tg = pd.read_excel(os.path.join(TEZ, "SHSY5Y_TDP43_qPCR_Ct_Data.xlsx"),
                   sheet_name="SOCE_Apoptosis_Genes_Ct", skiprows=4,
                   names=["gene", "group", "replicate", "date", "target_Ct", "GAPDH_Ct",
                          "dCt", "ddCt", "rel_expression", "note"])
# Keep the replicate rows only. The lab workbook ends with a summary block in Turkish (group
# means and an interpretation column) and carries per-replicate direction notes (UP/DOWN);
# neither is source data, so neither is published.
tg = tg[tg.group.isin(["Kontrol", "shTDP-43"]) & tg.gene.isin(TARGETS)].drop(columns="note").copy()
tg["group"] = tg.group.map({"Kontrol": NT, "shTDP-43": "shTDP-43"})

rel = []
for g, tab in qp.items():
    if g not in TARGETS:
        continue
    for grp, vals in tab.items():
        for i, v in enumerate(vals, 1):
            rel.append(dict(gene=g, group=(NT if grp == "Control" else "shTDP-43"),
                            replicate=i, rel_expression=v))
rel = pd.DataFrame(rel)

fura = []
for tab, lab in [("ER_Ca2_release", "ER Ca2+ release"), ("SOCE", "SOCE")]:
    for grp, vals in fu[tab].items():
        for i, v in enumerate(vals, 1):
            fura.append(dict(measurement=lab,
                             group=(NT if grp == "Control" else "shTDP-43"),
                             replicate=i, delta_F340_F380=v))
fura = pd.DataFrame(fura)

wst = []
for tab, lab in [("WST_1_48h", "48 h")]:
    for grp, vals in ws[tab].items():
        for i, v in enumerate(vals, 1):
            wst.append(dict(time=lab, group=(NT if grp == "Control" else "shTDP-43"),
                            well=i, signal_pct_of_control=v))
wst = pd.DataFrame(wst)

# ------------------------------------------------------------------ summaries
rows = []
for g in ["Untransduced control", "Non-targeting shRNA control", "shTDP-43"]:
    v = td.loc[td.group == g, "rel_expression"].astype(float).values
    rows.append(dict(panel="1A", measurement="TARDBP relative expression", group=g,
                     n=len(v), mean=round(v.mean(), 4), SEM=round(sem(v), 4), test="", p_value=""))
un = td.loc[td.group == "Untransduced control", "rel_expression"].astype(float).values
nt = td.loc[td.group == "Non-targeting shRNA control", "rel_expression"].astype(float).values
kd = td.loc[td.group == "shTDP-43", "rel_expression"].astype(float).values
F, pa = stats.f_oneway(np.log2(un), np.log2(nt), np.log2(kd))
rows.append(dict(panel="1A", measurement="TARDBP silencing", group="shTDP-43 vs controls",
                 n=4, mean=round(100 * (1 - kd.mean() / nt.mean()), 1), SEM="",
                 test="one-way ANOVA on log2 values, Tukey post hoc", p_value=f"{pa:.3g}"))
genes = ["TRPC1", "STIM1", "ORAI1", "ATP2A3"]
raw_p = {}
for g in genes:
    c = np.array(qp[g]["Control"]); k = np.array(qp[g]["shTDP-43"])
    raw_p[g] = stats.ttest_ind(c, k).pvalue
order = sorted(raw_p, key=raw_p.get)
holm = {}
running = 0.0
for rank, g in enumerate(order):
    running = max(running, raw_p[g] * (len(genes) - rank))
    holm[g] = min(1.0, running)
for g in genes:
    c = np.array(qp[g]["Control"]); k = np.array(qp[g]["shTDP-43"])
    rows.append(dict(panel="1B", measurement=f"{g} relative mRNA", group=NT,
                     n=len(c), mean=round(c.mean(), 3), SEM=round(sem(c), 3), test="", p_value=""))
    rows.append(dict(panel="1B", measurement=f"{g} relative mRNA", group="shTDP-43",
                     n=len(k), mean=round(k.mean(), 3), SEM=round(sem(k), 3),
                     test="two-tailed Student's t-test; Holm-adjusted across four targets",
                     p_value=f"{holm[g]:.6g}"))
# Figure 1C is WST-1; Figure 2C-D are Fura-2 group amplitudes.
for tab, lab in [("WST_1_48h", "48 h")]:
    c = np.array(ws[tab]["Control"]); k = np.array(ws[tab]["shTDP-43"])
    rows.append(dict(panel="1C", measurement=f"WST-1 signal {lab}", group=NT,
                     n=len(c), mean=round(c.mean(), 1), SEM=round(sem(c), 1), test="", p_value=""))
    rows.append(dict(panel="1C", measurement=f"WST-1 signal {lab}", group="shTDP-43",
                     n=len(k), mean=round(k.mean(), 1), SEM=round(sem(k), 1),
                     test="descriptive only; four wells from one experiment", p_value=""))
for tab, lab, panel in [("ER_Ca2_release", "ER Ca2+ release", "2C"),
                        ("SOCE", "SOCE", "2D")]:
    c = np.array(fu[tab]["Control"]); k = np.array(fu[tab]["TDP-43 KD"])
    p = stats.ttest_ind(c, k).pvalue
    rows.append(dict(panel=panel, measurement=f"{lab} delta F340/F380", group=NT,
                     n=len(c), mean=round(c.mean(), 3), SEM=round(sem(c), 3), test="", p_value=""))
    rows.append(dict(panel=panel, measurement=f"{lab} delta F340/F380", group="shTDP-43",
                     n=len(k), mean=round(k.mean(), 3), SEM=round(sem(k), 3),
                     test="two-tailed Student's t-test; n = 3",
                     p_value=f"{p:.8g}"))
summ = pd.DataFrame(rows)

readme = pd.DataFrame({"sheet": ["TARDBP_qPCR", "Target_qPCR_Ct", "Target_qPCR_rel",
                                 "Fura2", "WST1", "Summary_stats", "Primers", "Thermal_profile"],
    "content": [
        "RT-qPCR of TARDBP in three groups (n = 4 each); raw Ct for "
        "TARDBP and GAPDH, dCt, ddCt and 2^-ddCt. Experiment dates 14.04-21.05.2026.",
        "Raw Ct values for the four SOCE-associated targets and GAPDH in shTDP-43 cells and "
        "the non-targeting (scrambled) shRNA control, n = 4, "
        "same RNA set; experiment window 03.06-10.07.2026.",
        "Relative expression per replicate (2^-ddCt) for the four targets plotted in Figure 1B.",
        "Fura-2/AM measurements. ER Ca2+ release is the rise in F340/F380 after 10 uM "
        "cyclopiazonic acid in Ca2+-free HBS with EGTA; SOCE is the rise after re-addition "
        "of 1.5 mM CaCl2. Both as delta(F340/F380) versus the preceding baseline (n = 3). The control group is the "
        "non-targeting (scrambled) shRNA control.",
        "WST-1 metabolic signal at 48 h (n = 4; not a direct cell count or viability measure), normalised to the mean of the non-targeting (scrambled) shRNA control.",
        "Group means, SEM and the statistical test behind Figures 1 and 2.",
        "Primer sequences, product sizes and annealing temperatures for the RT-qPCR targets "
        "and the GAPDH reference.",
        "Thermal cycling profile of the RT-qPCR reactions."]})

primers = pd.DataFrame([
    ["NM_007375.4", "TARDBP", "GATGGTGTGACTGCAAACTTC", "CAGCTCATCCTCAGTCATGTC", 110, 60],
    ["NM_003156.4", "STIM1", "AGCAGAGTTTTGCCGAATTG", "ATCACTTTCTTCCACATCCACAT", 133, 60],
    ["NM_032790.4", "ORAI1", "CAGAGTTACTCCGAGGTGATGAG", "GAGAGCAGAGCCGAGGTCC", 119, 60],
    ["NM_003304.5", "TRPC1", "TGCGACAAGGGTGACTATTA", "TCCATTAGTTTCTGACAACCG", 176, 57],
    ["NM_005173.4", "ATP2A3 (SERCA3)", "GTCATCAACATCGGCCACTT", "GCCAGGCATGTAGTGATGAC", 143, 60],
    ["NM_002046.7", "GAPDH", "CAGTCAGCCGCATCTTCTTT", "GCCCAATACGACCAAATCC", 192, 60]],
    columns=["accession", "gene", "forward_primer_5_3", "reverse_primer_5_3",
             "product_size_bp", "annealing_temperature_C"])
cycling = pd.DataFrame([
    ["Initial denaturation", "95 °C, 5 min", 1],
    ["Amplification", "95 °C 10 s; annealing 57 or 60 °C 20 s (per gene); 72 °C 10 s", 40],
    ["Extension", "72 °C, 10 min", 1],
    ["Melt curve", "95 °C 10 s; 63 °C 15 s; continuous ramp to 95 °C", 1]],
    columns=["step", "thermal_condition", "cycles"])

p = os.path.join(SUP, "S1_laboratory_source_data.xlsx")
with pd.ExcelWriter(p, engine="openpyxl") as xw:
    readme.to_excel(xw, sheet_name="README", index=False)
    td.to_excel(xw, sheet_name="TARDBP_qPCR", index=False)
    tg.to_excel(xw, sheet_name="Target_qPCR_Ct", index=False)
    rel.to_excel(xw, sheet_name="Target_qPCR_rel", index=False)
    fura.to_excel(xw, sheet_name="Fura2", index=False)
    wst.to_excel(xw, sheet_name="WST1", index=False)
    summ.to_excel(xw, sheet_name="Summary_stats", index=False)
    primers.to_excel(xw, sheet_name="Primers", index=False)
    cycling.to_excel(xw, sheet_name="Thermal_profile", index=False)
print("written:", p)
for n, d in [("TARDBP_qPCR", td), ("Target_qPCR_Ct", tg), ("Target_qPCR_rel", rel),
             ("Fura2", fura), ("WST1", wst), ("Summary_stats", summ)]:
    print(f"  {n:16s} {len(d):4d} rows")

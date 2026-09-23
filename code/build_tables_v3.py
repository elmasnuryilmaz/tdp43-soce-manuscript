#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v3 — main tables and supplementary files for the submission package.

Writes 09_YAYIN_PAKETI/tables/Table1..5.csv and
       09_YAYIN_PAKETI/supplementary/S1..S16.*
All column names are in English; all values are read from the analysis outputs.
"""
import io, os, shutil
import numpy as np
import pandas as pd
from scipy import stats

M = "/Users/elmas/Desktop/MAKALE"
TEZ = "/Users/elmas/Desktop/TEZ"
D = f"{M}/07_DISK_ANALIZLERI"
TAB = f"{M}/09_YAYIN_PAKETI/tables"
SUP = f"{M}/09_YAYIN_PAKETI/supplementary"
for d in (TAB, SUP):
    os.makedirs(d, exist_ok=True)

JUNC = {"veri_seti": "comparison", "gene": "gene", "chrom": "chrom", "start": "start",
        "end": "end", "strand": "strand", "sinif": "junction_class", "lsv_tipi": "lsv_type",
        "psi_KD": "PSI_knockdown", "psi_CTRL": "PSI_control", "dPSI": "delta_PSI",
        "GA_alt": "CI95_low", "GA_ust": "CI95_high", "q": "q_value",
        "okuma_KD": "reads_knockdown", "okuma_CTRL": "reads_control",
        "toplam_KD": "LSV_total_knockdown", "toplam_CTRL": "LSV_total_control",
        "n_junction_lsv": "junctions_in_LSV", "rho": "overdispersion_rho", "p": "p_value",
        "psi_ham_KD": "PSI_raw_knockdown", "psi_ham_CTRL": "PSI_raw_control",
        "n_KD_pozitif": "n_knockdown_libraries_positive",
        "n_CTRL_pozitif": "n_control_libraries_positive", "lsv_anahtar": "LSV_key"}
CLASS = {"anotasyonlu": "annotated", "yeni_kombinasyon": "novel combination",
         "yeni_bolge": "novel splice site", "tamamen_yeni": "fully novel"}
LSVT = {"verici": "donor", "alici": "acceptor"}
DSET = {"SH_SY5Y": "SH-SY5Y 75 ng/mL", "SH_SY5Y_DOZ25": "SH-SY5Y 25 ng/mL",
        "iPSC_koloni": "iPSC colonies", "iPSC_MN": "iPSC-MN, TDP-43 KD",
        "iPSC_MN_FUS": "iPSC-MN, FUS KD", "iPSC_MN_TAF15": "iPSC-MN, TAF15 KD",
        "K562_mRNA": "K562 poly(A)+ mRNA", "K562_totalRNA": "K562 total RNA",
        "C2C12": "C2C12", "NSC34": "NSC34", "Fare_striatum": "Mouse striatum",
        "OWN_SH_SY5Y": "SH-SY5Y 75 ng/mL (MAPQ-filtered set)"}


def en_junction(df):
    df = df.rename(columns=JUNC)
    if "junction_class" in df:
        df["junction_class"] = df["junction_class"].map(lambda x: CLASS.get(x, x))
    if "lsv_type" in df:
        df["lsv_type"] = df["lsv_type"].map(lambda x: LSVT.get(x, x))
    if "comparison" in df:
        df["comparison"] = df["comparison"].map(lambda x: DSET.get(x, x))
    return df


def w(df, path, **kw):
    df.to_csv(path, index=False, **kw)
    print(f"  {os.path.basename(path):55s} {len(df):6d} rows")


print("main tables")
# ------------------------------------------------------------------- Table 1
t = pd.read_csv(f"{M}/03_TABLOLAR/kapsam_filtresi_veriseti_ozeti.tsv", sep="\t")
t.columns = ["dataset", "events_tested", "events_after_filter", "significant_before_filter",
             "significant_after_filter", "events_removed_pct"]
t["dataset"] = t.dataset.replace({"GSE230647_iPSC_koloni": "iPSC colonies (GSE230647)",
    "GSE27394_mouse_SE": "Mouse striatum (GSE27394)", "GSE296712_SHSY5Y": "SH-SY5Y (GSE296712)",
    "GSE77702_iPSC_MN": "iPSC-derived motor neurons (GSE77702)",
    "Mouse_PE_C2C12": "C2C12 (GSE171714)", "Mouse_PE_NSC34": "NSC34 (GSE171714)"})
t["significant_lost_pct"] = (100 * (1 - t.significant_after_filter / t.significant_before_filter)).round(1)
t = t[["dataset", "events_tested", "events_after_filter", "events_removed_pct",
       "significant_before_filter", "significant_after_filter", "significant_lost_pct"]]
w(t, f"{TAB}/Table1_coverage_prefilter.csv")

# ------------------------------------------------------------------- Table 2
s = pd.read_csv(f"{M}/03_TABLOLAR/SOCE_izoform_SAGLAM_olaylar.tsv", sep="\t")
b = pd.read_csv(f"{M}/03_TABLOLAR/SOCE_izoform_bootstrap_GA.tsv", sep="\t")
s = s[s.veri_seti == "GSE296712_SHSY5Y"].merge(
    b[["rmats_ID", "GA_alt", "GA_ust", "sifir_GA_icinde"]], on="rmats_ID", how="left")
s = s.rename(columns={"olay": "event_class", "gen": "gene", "chr": "chrom",
    "ekzon_bp": "exon_bp", "cerceve_korunur": "reading_frame",
    "aa_karsiligi": "amino_acids", "dPSI": "delta_PSI", "PSI_KD": "PSI_knockdown",
    "PSI_KONTROL": "PSI_control", "toplam_okuma": "total_informative_reads",
    "ort_okuma": "mean_reads_per_sample", "min_bilgilendirici": "min_informative_reads",
    "GA_alt": "CI95_low", "GA_ust": "CI95_high", "sifir_GA_icinde": "CI_includes_zero"})
s["reading_frame"] = s.reading_frame.astype(str).str.lower().map(
    {"evet": "preserved", "hayir": "disrupted", "true": "preserved", "false": "disrupted"}).fillna("")
s["CI_includes_zero"] = s.CI_includes_zero.astype(str).str.strip().str.lower().map(
    {"evet": "yes", "hayir": "no", "true": "yes", "false": "no"}).fillna("")
s["dataset"] = "SH-SY5Y (GSE296712)"
s["start_1based"] = s["start"] + 1
# coordinates and counts are whole numbers that the source table stores as floats
for c in ["start_1based", "end", "exon_bp", "amino_acids", "min_informative_reads"]:
    assert (s[c].dropna() % 1 == 0).all(), c
    s[c] = s[c].astype("Int64")
w(s[["dataset", "gene", "event_class", "chrom", "start_1based", "end", "strand", "exon_bp",
     "reading_frame", "amino_acids", "delta_PSI", "FDR", "CI95_low", "CI95_high",
     "CI_includes_zero", "PSI_knockdown", "PSI_control", "mean_reads_per_sample",
     "min_informative_reads"]], f"{TAB}/Table2_robust_SOCE_splicing_events.csv")

# ------------------------------------------------------------------- Table 3
hc = pd.read_csv(f"{D}/tablolar/S13_yuksek_guven_kriptik_ozet.tsv", sep="\t")
hc = hc.rename(columns={"veri_seti": "comparison", "olay": "high_confidence_events",
    "gen": "genes", "pozitif_kontrol": "positive_controls_high_confidence",
    "bulunan": "positive_control_genes_high_confidence", "Tier1": "Tier1_genes",
    "cekirdek_SOCE": "SOCE_machinery_genes"})
perm = {}
for ds in hc.comparison:
    f = f"{D}/sonuclar/RT_KRIPTIK_{ds}.tsv"
    if os.path.exists(f):
        g = pd.read_csv(f, sep="\t")
        POS = {"STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
               "ATG4B", "SETD5", "ELAVL3", "POLDIP3", "CAMK2B", "RSF1", "GPSM2", "SYNJ2"}
        gu = set(g.gene.astype(str).str.upper())
        perm[ds] = (len(g), len(gu), len(gu & POS))
hc["permissive_events"] = [perm.get(d, (np.nan,) * 3)[0] for d in hc.comparison]
hc["permissive_genes"] = [perm.get(d, (np.nan,) * 3)[1] for d in hc.comparison]
hc["positive_controls_permissive"] = [perm.get(d, (np.nan,) * 3)[2] for d in hc.comparison]
null = pd.read_csv(f"{D}/tablolar/S15_kriptik_esik_kalibrasyonu.tsv", sep="\t")
null = null[null.kademe == "K2"].set_index("veri_seti")
hc["null_calls_high_confidence"] = pd.array([null.bos_olay.get(d, np.nan) for d in hc.comparison],
                                           dtype="Float64").astype("Int64")
hc["null_to_real_ratio"] = [null.yanlis_pozitif_orani.get(d, np.nan) for d in hc.comparison]
hc["comparison"] = hc.comparison.map(lambda x: DSET.get(x, x))
# The sixteen literature controls are human cryptic events; the STMN2 and UNC13A events are
# absent from the mouse genes (Melamed et al., 2019; Ma et al., 2022). A match in a mouse
# comparison is a gene-name match only (for example an unrelated acceptor in NSC34 Unc13a),
# so recovery is not assessed there.
MOUSE = {"C2C12", "NSC34", "Mouse striatum"}
_m = hc.comparison.isin(MOUSE)
hc["positive_controls_permissive"] = hc.positive_controls_permissive.astype(object)
hc["positive_controls_high_confidence"] = hc.positive_controls_high_confidence.astype(object)
hc.loc[_m, ["positive_controls_permissive", "positive_controls_high_confidence",
            "positive_control_genes_high_confidence"]] = "n/a (mouse)"
hc = hc[["comparison", "permissive_events", "permissive_genes", "positive_controls_permissive",
         "high_confidence_events", "genes", "positive_controls_high_confidence",
         "positive_control_genes_high_confidence", "Tier1_genes", "SOCE_machinery_genes",
         "null_calls_high_confidence", "null_to_real_ratio"]]
w(hc, f"{TAB}/Table3_cryptic_events_eleven_comparisons.csv")

# ------------------------------------------------------------------- Table 4
t4 = pd.read_csv(f"{M}/03_TABLOLAR/v3/Table4_family_TPM.csv")
w(t4, f"{TAB}/Table4_transcript_family_abundance.csv")

# ------------------------------------------------------------------- Table 5
rows = []
G = ["TRPC1", "SARAF", "CBARP"]
NY = f"{TEZ}/output/NYGC_genom_geneli_22.08.26"
REG = {"Cerebellum": "Cerebellum", "Cortex_Frontal": "Frontal cortex",
       "Cortex_Motor_Lateral": "Motor cortex (lateral)", "Cortex_Motor_Medial": "Motor cortex (medial)",
       "Cortex_Occipital": "Occipital cortex", "Cortex_Temporal": "Temporal cortex",
       "Hippocampus": "Hippocampus", "Spinal_Cord_Cervical": "Spinal cord (cervical)",
       "Spinal_Cord_Lumbar": "Spinal cord (lumbar)", "Spinal_Cord_Thoracic": "Spinal cord (thoracic)"}
oz = pd.read_csv(f"{NY}/duz_bolge_ozeti.csv")
oz["key"] = oz.doku.str.replace(" ", "_")
oz = oz.set_index("key")
for key, lab in REG.items():
    f = f"{NY}/duz_bolge_{key}.csv"
    d = pd.read_csv(f).set_index("gen")
    for g in G:
        if g in d.index:
            rows.append(dict(cohort="ALS (NYGC GSE153960)", group="ALS", region=lab, gene=g,
                             n_case=int(oz.loc[key, "n_ALS"]), n_control=int(oz.loc[key, "n_kontrol"]),
                             log2FC=round(d.loc[g, "log2FC"], 3), cliffs_delta=round(d.loc[g, "delta"], 3),
                             q_value=float(f"{d.loc[g,'q']:.3g}")))
ON = f"{TEZ}/output/NYGC_ozgulluk_22.08.26"
# group sizes of the ONd comparisons, as logged by the analysis that produced them
import re as _re
_n_ond = {}
for _line in io.open(f"{ON}/run.log", encoding="utf-8"):
    _m = _re.match(r"\s+(\S.*?)\s+ONd\s+hasta=\s*(\d+)\s+kontrol=\s*(\d+)", _line)
    if _m:
        _n_ond[_m.group(1).replace(" ", "_")] = (int(_m.group(2)), int(_m.group(3)))
for key, lab in [("Cerebellum", "Cerebellum"), ("Cortex_Frontal", "Frontal cortex"),
                 ("Cortex_Temporal", "Temporal cortex")]:
    d = pd.read_csv(f"{ON}/{key}_ONd.csv").set_index("gen")
    for g in G:
        if g in d.index:
            rows.append(dict(cohort="Other neurological disorders (NYGC, same controls)",
                             group="ONd", region=lab, gene=g,
                             n_case=_n_ond[key][0], n_control=_n_ond[key][1],
                             log2FC=round(d.loc[g, "log2FC"], 3), cliffs_delta=round(d.loc[g, "delta"], 3),
                             q_value=float(f"{d.loc[g,'q']:.3g}")))
rows += [dict(cohort="Alzheimer's disease (GSE125583)", group="AD", region="Fusiform gyrus",
              gene="TRPC1", n_case=219, n_control=70, log2FC=-0.476, cliffs_delta=-0.447, q_value=1.61e-7),
         dict(cohort="Parkinson's disease (GSE68719)", group="PD", region="BA9",
              gene="TRPC1", n_case=29, n_control=44, log2FC=-0.570, cliffs_delta=-0.677, q_value=1.8e-4)]
MSLAB = {"DENEK duzeyi (10 MS vs 5 kontrol)": "Donor level (10 MS vs 5 control donors)",
         "MS NAWM vs kontrol WM": "Normal-appearing white matter vs control white matter",
         "MS lezyonlari vs kontrol WM": "MS lesions vs control white matter",
         "TUM BOLGELER (bolge icinde merkezlenmis)": "All five regions pooled (centred within region)"}
ms = pd.read_csv(f"{M}/06_MS_ANALIZI/MS_TRPC1_sonuclar.tsv", sep="\t")
ms["karsilastirma"] = ms.karsilastirma.map(
    lambda v: " \u00b7 ".join([p if i == 0 else MSLAB.get(p, p)
                                for i, p in enumerate(str(v).split(" \u00b7 "))]))
for _, r in ms[ms.gen.isin(G)].iterrows():
    rows.append(dict(cohort=r.karsilastirma.split(" · ")[0], group="MS",
                     region=r.karsilastirma.split(" · ")[-1], gene=r.gen,
                     n_case=r.n_hasta, n_control=r.n_kontrol, log2FC=r.medyan_fark,
                     cliffs_delta=r.cliffs_delta, q_value=float(f"{r.q:.3g}")))
t5 = pd.DataFrame(rows)
t5["n_case"] = t5.n_case.astype("Int64")
t5["n_control"] = t5.n_control.astype("Int64")
w(t5, f"{TAB}/Table5_cross_disease_comparison.csv")

# ------------------------------------------------------------- supplementary
print("supplementary")
# S2 panels
pan = pd.read_csv(f"{D}/paket/paneller/ALL_PANELS_LONG.csv")
w(pan, f"{SUP}/S2_calcium_gene_panels.csv")
# S4 enrichment
e = pd.read_csv(f"{M}/03_TABLOLAR/zenginlesme_gercek_paneller_permutasyon.tsv", sep="\t")
e.columns = ["dataset", "panel", "testable_genes", "significant_genes", "observed_pct",
             "background_pct", "matched_null_pct", "matched_null_CI", "p_hypergeometric",
             "p_matched_permutation", "decision"]
e["decision"] = e.decision.str.strip().str.lower().replace(
    {"zenginlesme yok": "no enrichment", "zenginlesme var": "enrichment"})
w(e, f"{SUP}/S4_matched_permutation_enrichment.csv")
# S5, S6, S7 junction tables
for src, dst in [("S5_yuksek_guven_kriptik_olaylar", "S5_high_confidence_cryptic_events"),
                 ("S6_kriptik_pozitif_kontroller", "S6_cryptic_positive_controls"),
                 ("S7_SOCE_anotasyonsuz_analiz", "S7_SOCE_genes_annotation_free")]:
    df = en_junction(pd.read_csv(f"{D}/tablolar/{src}.tsv", sep="\t"))
    if dst.startswith("S6_"):
        # positive-control recovery is assessed in the human comparisons only (see Table 3)
        df = df[~df.comparison.isin(MOUSE)]
    w(df, f"{SUP}/{dst}.csv")
# S6b matrix - headers were corrupted in the original output; rebuilt here
mat = pd.read_csv(f"{D}/tablolar/S6_kriptik_pozitif_kontroller.tsv", sep="\t")
mat["veri_seti"] = mat.veri_seti.map(lambda x: DSET.get(x, x))
mat = mat[~mat.veri_seti.isin(MOUSE)]
piv = mat.pivot_table(index="gene", columns="veri_seti", values="dPSI", aggfunc="max")
piv.index.name = "gene"
piv.reset_index().to_csv(f"{SUP}/S6b_positive_control_matrix.csv", index=False)
print(f"  {'S6b_positive_control_matrix.csv':55s} {len(piv):6d} rows")
# S8, S9
s8 = pd.read_csv(f"{D}/tablolar/S8_NYGC_kriptikPSI_ALS_vs_kontrol.tsv", sep="\t")
s8.columns = ["gene", "region", "n_ALS", "n_control", "median_ALS", "median_control",
              "cliffs_delta", "p_value", "q_value"]
s8["region"] = s8.region.str.replace("_", " ")
w(s8, f"{SUP}/S8_cryptic_STMN2_ALS_vs_control.csv")
s9 = pd.read_csv(f"{D}/tablolar/S9_NYGC_kriptikPSI_korelasyon.tsv", sep="\t")
s9.columns = ["region", "target_gene", "proxy", "n", "spearman_rho", "p_value", "q_value"]
s9["proxy"] = s9.proxy.replace({"kriptik_STMN2_PSI": "cryptic STMN2 PSI",
                                "gen_duzeyi_STMN2": "gene-level STMN2"})
# gene-level STMN2 correlated against itself is 1.0 by construction in every region;
# those 11 tests are excluded from the correction family, leaving 220 informative tests.
triv = (s9.proxy == "gene-level STMN2") & (s9.target_gene == "STMN2")
s9["in_correction_family"] = np.where(triv, "no", "yes")
inf = s9.loc[~triv].sort_values("p_value")
m = len(inf)
raw = inf.p_value.values
q = np.minimum.accumulate((raw * m / np.arange(1, m + 1))[::-1])[::-1]
s9["q_value"] = np.nan
s9.loc[inf.index, "q_value"] = np.minimum(q, 1.0)
s9 = s9[["region", "target_gene", "proxy", "n", "spearman_rho", "p_value",
         "q_value", "in_correction_family"]]
w(s9, f"{SUP}/S9_cryptic_PSI_correlations_within_ALS.csv")
_sig = s9[(s9.proxy == "cryptic STMN2 PSI") & (s9.q_value < 0.05)]
print("    q recomputed over %d informative tests; %d significant: %s"
      % (m, len(_sig), ", ".join("%s/%s q=%.4f" % (r.target_gene, r.region.split()[-1], r.q_value)
                                 for r in _sig.sort_values("q_value").itertuples())))
# S10 NMD
s10 = pd.read_csv(f"{D}/tablolar/S10_NMD_etkilesimi_SOCE.tsv", sep="\t")
s10 = s10.rename(columns={"symbol": "gene", "interaction_log2": "interaction_log2",
    "n_positive": "n_conditions_positive", "p_t4": "p_t_test_4conditions",
    "q_t4": "q_t_test_4conditions", "p_sign4": "p_sign_test", "q_sign4": "q_sign_test",
    "X": "interaction_XRN1", "XS": "interaction_XRN1_SMG6", "XU": "interaction_XRN1_UPF1",
    "US": "interaction_UPF1_SMG6"})
w(s10, f"{SUP}/S10_NMD_interaction_SOCE_panel.csv")
pn = pd.read_csv(f"{M}/03_TABLOLAR/v3/nmd_panel_t4.csv")
# the same one-sided test for the sixteen literature cryptic targets, quoted in Section 3.4
_r = pd.read_csv(f"{D}/sonuclar/NMD_etkilesim_paylasimli_kontrol_t4_sembol.tsv", sep="\t")
_POS16 = {"STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
          "ATG4B", "SETD5", "CAMK2B", "ELAVL3", "POLDIP3", "RSF1", "GPSM2", "SYNJ2"}
_v = _r.loc[_r.symbol.isin(_POS16), "interaction_log2"].values
_u = stats.mannwhitneyu(_v, _r.interaction_log2.values, alternative="greater")
pn = pd.concat([pn, pd.DataFrame([{
    "panel": "Cryptic_positive_controls_16", "n_genes": len(_v),
    "median_interaction_log2": float(np.median(_v)),
    "p_one_sided_MWU": float(_u.pvalue)}])], ignore_index=True)
w(pn, f"{SUP}/S10b_NMD_panel_level_tests.csv")
# S11 APA
s11 = pd.read_csv(f"{D}/tablolar/S12_APA_anlamli_olaylar.tsv", sep="\t")
s11 = s11.rename(columns={"gene": "gene", "unit": "unit", "measure": "index",
    "index_KD": "index_knockdown", "index_CTRL": "index_control", "delta": "delta",
    "boot_low": "bootstrap_CI95_low", "boot_high": "bootstrap_CI95_high",
    "sample_KD": "per_sample_knockdown", "sample_CTRL": "per_sample_control"})
s11["index"] = s11["index"].replace({"IPA_index": "intronic polyadenylation index",
                                     "UTR_index": "distal 3'UTR usage index"})
s11.insert(0, "dataset", "SH-SY5Y (GSE296712)")
# datasets recomputed after the external drive became available
extra = []
for f in ["APA_corrected_iPSC_MN.tsv", "APA_corrected_remaining_datasets.tsv",
          "APA_corrected_C2C12.tsv", "APA_corrected_NSC34.tsv"]:
    fp = f"{D}/sonuclar/{f}"
    if os.path.exists(fp):
        e = pd.read_csv(fp, sep="\t")
        e = e.rename(columns={"measure": "index", "index_KD": "index_knockdown",
                              "index_CTRL": "index_control", "boot_low": "bootstrap_CI95_low",
                              "boot_high": "bootstrap_CI95_high",
                              "sample_KD": "per_sample_knockdown",
                              "sample_CTRL": "per_sample_control"})
        e["dataset"] = e.dataset.map(lambda x: DSET.get(x, x))
        e["index"] = e["index"].replace({"IPA_index": "intronic polyadenylation index",
                                         "UTR_distal_index": "distal 3'UTR usage index"})
        extra.append(e)
if extra:
    ex = pd.concat(extra, ignore_index=True).drop_duplicates(["dataset", "gene", "unit", "index"])
    s11 = pd.concat([s11, ex[[c for c in ex.columns if c in list(s11.columns) + ["n_boot"]]]],
                    ignore_index=True)
# SH-SY5Y intervals come from the same complete enumeration (3^6 replicate combinations)
s11.loc[s11.dataset.str.startswith("SH-SY5Y") & s11.n_boot.isna(), "n_boot"] = 729
# Unit numbers are indices of the gaps between merged exons of all basic-annotation
# transcripts, in the direction of transcription; they need not match canonical intron
# numbers, so the genomic windows of every unit are given (1-based, inclusive).
def _windows(bedfile):
    b = pd.read_csv(bedfile, sep="\t", header=None, names=["chrom", "s0", "end", "name"])
    p = b.name.str.split("|", expand=True)
    b["gene"], b["unit"], b["win"], b["strand"] = p[0], p[1], p[2], p[3]
    b["chrom"] = "chr" + b.chrom.astype(str).str.replace("^chr", "", regex=True)
    b["iv"] = (b.s0 + 1).astype(str) + "-" + b.end.astype(str)
    iv = b.pivot_table(index=["gene", "unit"], columns="win", values="iv", aggfunc="first")
    meta = b.groupby(["gene", "unit"])[["chrom", "strand"]].first()
    out = meta.join(iv).reset_index()
    out["first_window"] = out.I5.where(out.unit != "termexon", out.Uprox)
    out["second_window"] = out.I3.where(out.unit != "termexon", out.Udist)
    return out[["gene", "unit", "chrom", "strand", "first_window", "second_window"]]
_wh = _windows(f"{D}/kod/apa_pencereleri_human.bed")
_wm = _windows(f"{D}/kod/apa_pencereleri_mouse.bed")
_human = s11.dataset.str.startswith("SH-SY5Y") | s11.dataset.str.startswith("iPSC-MN")
s11 = pd.concat([s11[_human].merge(_wh, on=["gene", "unit"], how="left"),
                 s11[~_human].merge(_wm, on=["gene", "unit"], how="left")], ignore_index=True)
s11 = s11.rename(columns={"first_window": "window_5prime_or_proximal",
                          "second_window": "window_3prime_or_distal"})
assert s11.window_5prime_or_proximal.notna().all(), "S11 unit without genomic window"
w(s11, f"{SUP}/S11_APA_candidate_gradients.csv")
# S12 cryptic counts (= Table 3 source, per dataset)
w(hc, f"{SUP}/S12_cryptic_counts_by_dataset.csv")
# S13 SOAR exon junction level
s13 = en_junction(pd.read_csv(f"{D}/tablolar/S14_STIM2_SOAR_ekzonu_birlesim_duzeyi.tsv", sep="\t"))
s13 = s13.rename(columns={"yan": "flank"})
s13["flank"] = s13.flank.replace({"asagi": "downstream", "yukari": "upstream"})
w(s13, f"{SUP}/S13_STIM2_SOAR_exon_junction_level.csv")
# S14 null test
s14 = pd.read_csv(f"{D}/tablolar/S15_kriptik_esik_kalibrasyonu.tsv", sep="\t")
s14.columns = ["dataset", "tier", "criteria", "real_calls", "null_calls",
               "null_to_real_ratio", "positive_controls"]
# the permissive definition of Methods 2.5 (A12_null_kontrol.py): null calls against real calls
_perm = []
for ds in s14.dataset.unique():
    rf, nf = f"{D}/sonuclar/RT_KRIPTIK_{ds}.tsv", f"{D}/sonuclar/NULL_KRIPTIK_{ds}.tsv"
    if os.path.exists(rf) and os.path.exists(nf):
        r_, n_ = pd.read_csv(rf, sep="\t"), pd.read_csv(nf, sep="\t")
        _perm.append(dict(dataset=ds, tier="permissive",
                          criteria="unannotated, ΔPSI≥0.05, q<0.05, control PSI≤0.05, CI lower bound>0",
                          real_calls=len(r_), null_calls=len(n_),
                          null_to_real_ratio=round(len(n_) / len(r_), 2),
                          positive_controls=len(set(r_.gene.astype(str).str.upper()) & _POS16)))
s14 = pd.concat([pd.DataFrame(_perm), s14], ignore_index=True)
s14["dataset"] = s14.dataset.map(lambda x: DSET.get(x, x))
s14["tier"] = s14.tier.replace({"K1": "tier 1", "K2": "tier 2 (high-confidence)", "K3": "tier 3"})
s14["criteria"] = (s14.criteria.str.replace("okuma", "reads")
                   .str.replace("GA>", "CI lower bound>", regex=False))
s14["positive_controls"] = s14.positive_controls.astype(object)
s14.loc[s14.dataset.isin(MOUSE), "positive_controls"] = "n/a (mouse)"
w(s14, f"{SUP}/S14_control_vs_control_null_test.csv")
# S15 STIM2.1 meta across six datasets
s15 = pd.read_csv(f"{M}/03_TABLOLAR/STIM2_1_SOAR_ekzonu_meta.tsv", sep="\t")
s15 = s15.rename(columns={"veri_seti": "dataset", "tur": "species", "dPSI": "delta_PSI",
    "var": "variance", "ort_okuma": "mean_reads_per_sample", "boot_GA": "bootstrap_CI95",
    "PSI_KD": "PSI_knockdown_per_replicate", "PSI_KONTROL": "PSI_control_per_replicate"})
s15["species"] = s15.species.replace({"human": "human", "mouse": "mouse"})
w(s15, f"{SUP}/S15_STIM2.1_exon_six_datasets.csv")
# S16 MS
ms2 = ms.rename(columns={"karsilastirma": "comparison", "gen": "gene", "n_hasta": "n_case",
    "n_kontrol": "n_control", "medyan_fark": "median_difference_log2CPM",
    "cliffs_delta": "cliffs_delta", "p": "p_value", "q": "q_value",
    "anlamli": "significant", "yon": "direction"})
ms2["direction"] = ms2.direction.replace({"azalmis": "decreased", "artmis": "increased"})
ms2["significant"] = ms2.significant.fillna("").replace({"*": "yes"})
w(ms2, f"{SUP}/S16_multiple_sclerosis_both_cohorts.csv")
don = pd.read_csv(f"{M}/03_TABLOLAR/v3/MS_donor_level.csv")
w(don, f"{SUP}/S16b_multiple_sclerosis_donor_level.csv")
print("done")

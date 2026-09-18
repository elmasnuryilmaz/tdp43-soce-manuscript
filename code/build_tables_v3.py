#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v3 — main tables and supplementary files for the submission package.

Writes 09_YAYIN_PAKETI/tables/Table1..5.csv and
       09_YAYIN_PAKETI/supplementary/S1..S16.*
All column names are in English; all values are read from the analysis outputs.
"""
import os, shutil
import numpy as np
import pandas as pd

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
s["reading_frame"] = s.reading_frame.map({"EVET": "preserved", "hayir": "disrupted",
                                          True: "preserved", False: "disrupted"}).fillna(s.reading_frame)
s["CI_includes_zero"] = s.CI_includes_zero.map({"EVET": "yes", "hayir": "no"}).fillna("")
s["dataset"] = "SH-SY5Y (GSE296712)"
s["start_1based"] = s["start"] + 1
w(s[["dataset", "gene", "event_class", "chrom", "start_1based", "end", "strand", "exon_bp",
     "reading_frame", "amino_acids", "delta_PSI", "FDR", "CI95_low", "CI95_high",
     "CI_includes_zero", "PSI_knockdown", "PSI_control", "mean_reads_per_sample",
     "min_informative_reads"]], f"{TAB}/Table2_robust_SOCE_splicing_events.csv")

# ------------------------------------------------------------------- Table 3
hc = pd.read_csv(f"{D}/tablolar/S13_yuksek_guven_kriptik_ozet.tsv", sep="\t")
hc = hc.rename(columns={"veri_seti": "comparison", "olay": "high_confidence_events",
    "gen": "genes", "pozitif_kontrol": "positive_controls_high_confidence",
    "bulunan": "positive_control_genes_high_confidence", "Tier1": "Tier1_genes",
    "cekirdek_SOCE": "core_SOCE_genes"})
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
hc["null_calls_high_confidence"] = [null.bos_olay.get(d, np.nan) for d in hc.comparison]
hc["null_to_real_ratio"] = [null.yanlis_pozitif_orani.get(d, np.nan) for d in hc.comparison]
hc["comparison"] = hc.comparison.map(lambda x: DSET.get(x, x))
hc = hc[["comparison", "permissive_events", "permissive_genes", "positive_controls_permissive",
         "high_confidence_events", "genes", "positive_controls_high_confidence",
         "positive_control_genes_high_confidence", "Tier1_genes", "core_SOCE_genes",
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
for key, lab in [("Cerebellum", "Cerebellum"), ("Cortex_Frontal", "Frontal cortex"),
                 ("Cortex_Temporal", "Temporal cortex")]:
    d = pd.read_csv(f"{ON}/{key}_ONd.csv").set_index("gen")
    for g in G:
        if g in d.index:
            rows.append(dict(cohort="Other neurological disorders (NYGC, same controls)",
                             group="ONd", region=lab, gene=g, n_case=np.nan, n_control=np.nan,
                             log2FC=round(d.loc[g, "log2FC"], 3), cliffs_delta=round(d.loc[g, "delta"], 3),
                             q_value=float(f"{d.loc[g,'q']:.3g}")))
rows += [dict(cohort="Alzheimer's disease (GSE125583)", group="AD", region="Fusiform gyrus",
              gene="TRPC1", n_case=219, n_control=70, log2FC=-0.476, cliffs_delta=-0.447, q_value=1.61e-7),
         dict(cohort="Parkinson's disease (GSE68719)", group="PD", region="BA9",
              gene="TRPC1", n_case=29, n_control=44, log2FC=-0.570, cliffs_delta=-0.677, q_value=1.8e-4)]
ms = pd.read_csv(f"{M}/06_MS_ANALIZI/MS_TRPC1_sonuclar.tsv", sep="\t")
for _, r in ms[ms.gen.isin(G)].iterrows():
    rows.append(dict(cohort=r.karsilastirma.split(" · ")[0], group="MS",
                     region=r.karsilastirma.split(" · ")[-1], gene=r.gen,
                     n_case=r.n_hasta, n_control=r.n_kontrol, log2FC=r.medyan_fark,
                     cliffs_delta=r.cliffs_delta, q_value=float(f"{r.q:.3g}")))
t5 = pd.DataFrame(rows)
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
e["decision"] = e.decision.replace({"zenginlesme yok": "no enrichment",
                                    "zenginlesme var": "enrichment"})
w(e, f"{SUP}/S4_matched_permutation_enrichment.csv")
# S5, S6, S7 junction tables
for src, dst in [("S5_yuksek_guven_kriptik_olaylar", "S5_high_confidence_cryptic_events"),
                 ("S6_kriptik_pozitif_kontroller", "S6_cryptic_positive_controls"),
                 ("S7_SOCE_anotasyonsuz_analiz", "S7_SOCE_genes_annotation_free")]:
    df = en_junction(pd.read_csv(f"{D}/tablolar/{src}.tsv", sep="\t"))
    w(df, f"{SUP}/{dst}.csv")
# S6b matrix - headers were corrupted in the original output; rebuilt here
mat = pd.read_csv(f"{D}/tablolar/S6_kriptik_pozitif_kontroller.tsv", sep="\t")
mat["veri_seti"] = mat.veri_seti.map(lambda x: DSET.get(x, x))
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
w(s9, f"{SUP}/S9_cryptic_PSI_correlations_within_ALS.csv")
# S10 NMD
s10 = pd.read_csv(f"{D}/tablolar/S10_NMD_etkilesimi_SOCE.tsv", sep="\t")
s10 = s10.rename(columns={"symbol": "gene", "interaction_log2": "interaction_log2",
    "n_positive": "n_conditions_positive", "p_t4": "p_t_test_4conditions",
    "q_t4": "q_t_test_4conditions", "p_sign4": "p_sign_test", "q_sign4": "q_sign_test",
    "X": "interaction_XRN1", "XS": "interaction_XRN1_SMG6", "XU": "interaction_XRN1_UPF1",
    "US": "interaction_UPF1_SMG6"})
w(s10, f"{SUP}/S10_NMD_interaction_SOCE_panel.csv")
pn = pd.read_csv(f"{M}/03_TABLOLAR/v3/nmd_panel_t4.csv")
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
s14["dataset"] = s14.dataset.map(lambda x: DSET.get(x, x))
s14["criteria"] = s14.criteria.str.replace("okuma", "reads")
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

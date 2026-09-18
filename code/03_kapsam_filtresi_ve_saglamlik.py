#!/usr/bin/env python3
"""
Oncelik 2.2 + 4.x — Kapsam on filtresi ile FDR'nin yeniden hesaplanmasi
ve Ca2+ adaylarinin saglamlik siralamasi.

Mantik: rMATS tum olaylara FDR uygular; dusuk kapsamli on binlerce olay
coklu test yukunu sisirir. Alan standardi, FDR'den ONCE okuma destegi
filtresi uygulamaktir.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from importlib.machinery import SourceFileLoader

core = SourceFileLoader("core", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "01_rmats_core.py")).load_module()
TAB = os.path.join(core.OUT, "03_TABLOLAR")
os.makedirs(TAB, exist_ok=True)

GROUPS = core.load_gene_groups()
MIN_MEAN_READS = 10      # ornek basina ortalama bilgilendirici okuma
MIN_INFORMATIVE = 5      # her ornekte en az bu kadar bilgilendirici okuma

summary = []
allrows = []

for ds in core.DATASETS:
    for ev in core.EVENTS:
        df = core.load_event(ds, ev, "JC")
        if df is None:
            continue
        df["gene_up"] = df["geneSymbol"].str.upper()
        n_all = len(df)

        # --- filtresiz (tezdeki yaklasim)
        sig_raw = ((df["FDR"] < 0.05) & (df["IncLevelDifference"].abs() >= 0.10)).sum()

        # --- kapsam on filtresi + FDR'nin YENIDEN hesaplanmasi
        keep = ((df["mean_reads_per_sample"] >= MIN_MEAN_READS)
                & (df["min_informative_reads"] >= MIN_INFORMATIVE)
                & df["PValue"].notna())
        sub = df[keep].copy()
        sub["FDR_yeni"] = core.bh_fdr(sub["PValue"].values)
        sig_new = ((sub["FDR_yeni"] < 0.05) & (sub["IncLevelDifference"].abs() >= 0.10)).sum()

        summary.append(dict(veri_seti=ds, olay=ev, toplam_olay=n_all,
                            filtre_sonrasi=len(sub),
                            anlamli_filtresiz=int(sig_raw),
                            anlamli_filtreli=int(sig_new)))
        sub["dataset"] = ds
        sub["eventType"] = ev
        allrows.append(sub)

sm = pd.DataFrame(summary)
tot = sm.groupby("veri_seti")[["toplam_olay", "filtre_sonrasi",
                               "anlamli_filtresiz", "anlamli_filtreli"]].sum().reset_index()
tot["test_yuku_azalmasi_%"] = (100 * (1 - tot.filtre_sonrasi / tot.toplam_olay)).round(1)
sm.to_csv(os.path.join(TAB, "kapsam_filtresi_olay_bazinda.tsv"), sep="\t", index=False)
tot.to_csv(os.path.join(TAB, "kapsam_filtresi_veriseti_ozeti.tsv"), sep="\t", index=False)

print("=" * 92)
print("KAPSAM ON FILTRESININ ETKISI (JC, tum olay turleri birlestirilmis)")
print(f"olcut: ornek basina ort. >= {MIN_MEAN_READS} okuma VE her orneklte >= {MIN_INFORMATIVE}")
print("=" * 92)
print(tot.to_string(index=False))

# ------------------------------------------------------------------ birlesik tablo
A = pd.concat(allrows, ignore_index=True)
A.to_pickle(os.path.join(TAB, "_filtreli_olaylar.pkl"))

# ------------------------------------------------------------------ Ca2+ adaylari
print("\n" + "=" * 92)
print("Ca2+ GRUP 3 ADAYLARI — kapsam filtresinden GECEN ve esikleri karsilayan olaylar")
print("=" * 92)

ca = A[A["gene_up"].isin(GROUPS["Grup3"])].copy()
ca_sig = ca[(ca["FDR_yeni"] < 0.05) & (ca["IncLevelDifference"].abs() >= 0.10)].copy()

# gen basina en iyi olay + kac veri setinde gorulmus
ca_sig["abs_dpsi"] = ca_sig["IncLevelDifference"].abs()
best = (ca_sig.sort_values("abs_dpsi", ascending=False)
             .groupby(["gene_up", "dataset"], as_index=False).first())
rank = (best.groupby("gene_up")
            .agg(veri_seti_sayisi=("dataset", "nunique"),
                 datasetler=("dataset", lambda s: ";".join(sorted(set(s)))),
                 max_abs_dPSI=("abs_dpsi", "max"),
                 medyan_abs_dPSI=("abs_dpsi", "median"),
                 min_FDR=("FDR_yeni", "min"),
                 toplam_okuma_max=("total_reads", "max"),
                 yonler=("IncLevelDifference",
                         lambda s: "artis" if (s > 0).all() else ("azalis" if (s < 0).all() else "karisik")))
            .reset_index()
            .sort_values(["veri_seti_sayisi", "max_abs_dPSI"], ascending=[False, False]))
rank.to_csv(os.path.join(TAB, "Ca_adaylari_saglamlik_siralamasi.tsv"), sep="\t", index=False)
print(rank.head(25).to_string(index=False))

# ------------------------------------------------------------------ TRPC1 filtre sonrasi
print("\n" + "-" * 92)
print("Kritik kontrol: TRPC1 ana olayi kapsam filtresinden geciyor mu?")
t = A[(A.dataset == "GSE296712_SHSY5Y") & (A.gene_up == "TRPC1") & (A.eventType == "SE")]
if t.empty:
    print("  HAYIR — TRPC1 SH-SY5Y SE olaylarinin hicbiri kapsam filtresini gecemedi.")
    df0 = core.load_event("GSE296712_SHSY5Y", "SE", "JC")
    df0["gene_up"] = df0["geneSymbol"].str.upper()
    m = df0[(df0.gene_up == "TRPC1") & (df0.ID == 33987)]
    if not m.empty:
        r = m.iloc[0]
        print(f"  ana olay (ID 33987): ornek basina ort. okuma = "
              f"{r['mean_reads_per_sample']:.1f} (esik {MIN_MEAN_READS}), "
              f"en dusuk bilgilendirici = {int(r['min_informative_reads'])} (esik {MIN_INFORMATIVE})")
else:
    print(t[["ID", "gene_up", "IncLevelDifference", "FDR", "FDR_yeni",
             "mean_reads_per_sample", "min_informative_reads"]].to_string(index=False))

print("\nTablolar ->", TAB)

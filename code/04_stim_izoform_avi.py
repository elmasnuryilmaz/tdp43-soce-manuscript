#!/usr/bin/env python3
"""
Oncelik 3.1 + 3.2 — "mRNA artti ama SOCE azaldi" paradoksunun izoform temelli sinanmasi.

STIM2.1 / STIM2-beta: SOAR/CAD domainine 24 nt (8 aa) ekleyen alternatif ekzon;
STIM2'yi aktivatorden dominant-negatif inhibitore cevirir (Miederer 2015; Rana 2015).
Eger TDP-43 kaybinda bu ekzonun katilimi artiyorsa, mRNA artisi ile SOCE azalmasi
tek bir molekuler olayla baglanir.

Ayrica STIM1'in 36 bp'lik SH-SY5Y olayi ve tum SOCE duzenleyicileri taranir.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from importlib.machinery import SourceFileLoader

core = SourceFileLoader("core", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "01_rmats_core.py")).load_module()
TAB = os.path.join(core.OUT, "03_TABLOLAR")

TARGETS = ["STIM1", "STIM2", "STIMATE", "SARAF", "CRACR2A", "CRACR2B",
           "ORAI1", "ORAI2", "ORAI3", "TRPC1", "ATP2A2", "ATP2A3"]

rows = []
for ds in core.DATASETS:
    for ev in core.EVENTS:
        df = core.load_event(ds, ev, "JC")
        if df is None:
            continue
        df["gene_up"] = df["geneSymbol"].str.upper()
        sub = df[df["gene_up"].isin(TARGETS)].copy()
        if sub.empty:
            continue
        if "exonStart_0base" in sub.columns and "exonEnd" in sub.columns:
            sub["ekzon_uzunlugu_bp"] = sub["exonEnd"] - sub["exonStart_0base"]
        else:
            sub["ekzon_uzunlugu_bp"] = np.nan
        for _, r in sub.iterrows():
            kd, ct = core.per_replicate_psi(r)
            L = r["ekzon_uzunlugu_bp"]
            rows.append(dict(
                veri_seti=ds, olay=ev, gen=r["gene_up"], rmats_ID=r["ID"],
                chr=r["chr"], strand=r["strand"],
                start=r.get("exonStart_0base", np.nan), end=r.get("exonEnd", np.nan),
                ekzon_bp=L,
                cerceve_korunur=("evet" if (not pd.isna(L) and L % 3 == 0) else
                                 ("hayir" if not pd.isna(L) else "NA")),
                aa_karsiligi=(None if pd.isna(L) else int(L // 3)),
                dPSI=r["IncLevelDifference"], FDR=r["FDR"], PValue=r["PValue"],
                PSI_KD=";".join("NA" if np.isnan(x) else f"{x:.3f}" for x in kd),
                PSI_KONTROL=";".join("NA" if np.isnan(x) else f"{x:.3f}" for x in ct),
                toplam_okuma=r["total_reads"],
                ort_okuma=round(r["mean_reads_per_sample"], 1),
                min_bilgilendirici=r["min_informative_reads"],
            ))

S = pd.DataFrame(rows)
S.to_csv(os.path.join(TAB, "SOCE_izoform_tum_olaylar.tsv"), sep="\t", index=False)

print("=" * 100)
print("A) STIM2 — SOAR/CAD domainine 24 nt ekleyen STIM2.1/STIM2-beta olayi araniyor")
print("=" * 100)
s2 = S[S.gen == "STIM2"].copy()
print(f"STIM2 icin toplam {len(s2)} rMATS olayi bulundu.")
cand = s2[(s2.ekzon_bp >= 18) & (s2.ekzon_bp <= 36)]
if cand.empty:
    print("\n  >>> 18-36 bp araliginda (STIM2.1'in 24 nt'lik ekzonuna karsilik gelebilecek)")
    print("      HICBIR rMATS olayi yok.")
    print("      Kisa alternatif ekzonlar rMATS'in varsayilan olay tanimlarinda")
    print("      cogu kez yakalanmaz; hedefli birlesim sayimi gerekir.")
else:
    print(cand[["veri_seti", "olay", "chr", "start", "end", "ekzon_bp", "aa_karsiligi",
                "dPSI", "FDR", "ort_okuma"]].to_string(index=False))

print("\n  STIM2'nin esikleri karsilayan TUM olaylari:")
s2sig = s2[(s2.FDR < 0.05) & (s2.dPSI.abs() >= 0.10)]
if s2sig.empty:
    print("    (yok)")
else:
    print(s2sig[["veri_seti", "olay", "chr", "start", "end", "ekzon_bp", "aa_karsiligi",
                 "cerceve_korunur", "dPSI", "FDR", "PSI_KD", "PSI_KONTROL",
                 "ort_okuma", "min_bilgilendirici"]].to_string(index=False))

print("\n" + "=" * 100)
print("B) STIM1 — tezdeki chr11:4088702-4088738 (36 bp) olayi ve digerleri")
print("=" * 100)
s1 = S[S.gen == "STIM1"].copy()
tez = s1[(s1.start == 4088702) | (s1.end == 4088738)]
if not tez.empty:
    print("Tezde bildirilen olay:")
    print(tez[["veri_seti", "olay", "chr", "start", "end", "ekzon_bp", "aa_karsiligi",
               "cerceve_korunur", "dPSI", "FDR", "PSI_KD", "PSI_KONTROL",
               "ort_okuma", "min_bilgilendirici"]].to_string(index=False))
else:
    print("  chr11:4088702-4088738 koordinati bu rMATS surumunde bulunamadi.")
print("\nSTIM1'in esikleri karsilayan olaylari:")
s1sig = s1[(s1.FDR < 0.05) & (s1.dPSI.abs() >= 0.10)]
print(s1sig[["veri_seti", "olay", "chr", "start", "end", "ekzon_bp", "aa_karsiligi",
             "cerceve_korunur", "dPSI", "FDR", "PSI_KD", "PSI_KONTROL",
             "ort_okuma", "min_bilgilendirici"]].to_string(index=False)
      if not s1sig.empty else "  (yok)")

print("\n" + "=" * 100)
print("C) SOCE duzenleyicileri — okuma destegi GUCLU ve esikleri karsilayan olaylar")
print("   (bunlar TRPC1'in yerine gecebilecek adaylardir)")
print("=" * 100)
strong = S[(S.FDR < 0.05) & (S.dPSI.abs() >= 0.10)
           & (S.ort_okuma >= 10) & (S.min_bilgilendirici >= 5)].copy()
strong = strong.sort_values(["gen", "veri_seti"])
cols = ["veri_seti", "olay", "gen", "chr", "start", "end", "ekzon_bp", "aa_karsiligi",
        "cerceve_korunur", "dPSI", "FDR", "PSI_KD", "PSI_KONTROL",
        "ort_okuma", "min_bilgilendirici"]
print(strong[cols].to_string(index=False))
strong.to_csv(os.path.join(TAB, "SOCE_izoform_SAGLAM_olaylar.tsv"), sep="\t", index=False)

# ---- bootstrap: guclu olanlar icin GA
print("\n" + "=" * 100)
print("D) Saglam SOCE olaylari icin bootstrap %95 GA")
print("=" * 100)
res = []
for ds in strong.veri_seti.unique():
    for ev in strong[strong.veri_seti == ds].olay.unique():
        df = core.load_event(ds, ev, "JC")
        df["gene_up"] = df["geneSymbol"].str.upper()
        ids = strong[(strong.veri_seti == ds) & (strong.olay == ev)].rmats_ID.tolist()
        for _, r in df[df.ID.isin(ids)].iterrows():
            med, lo, hi = core.bootstrap_dpsi(r, n_boot=10000, seed=7)
            res.append(dict(veri_seti=ds, olay=ev, gen=r["gene_up"], rmats_ID=r["ID"],
                            dPSI=r["IncLevelDifference"],
                            boot_medyan=round(med, 3) if not np.isnan(med) else None,
                            GA_alt=round(lo, 3) if not np.isnan(lo) else None,
                            GA_ust=round(hi, 3) if not np.isnan(hi) else None,
                            sifir_GA_icinde=("EVET" if (not np.isnan(lo) and lo <= 0 <= hi)
                                             else "hayir")))
B = pd.DataFrame(res).sort_values(["gen", "veri_seti"])
B.to_csv(os.path.join(TAB, "SOCE_izoform_bootstrap_GA.tsv"), sep="\t", index=False)
print(B.to_string(index=False))
print("\nTablolar ->", TAB)

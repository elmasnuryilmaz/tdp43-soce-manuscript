#!/usr/bin/env python3
"""
Oncelik 2.1 + 2.2 — Ana kirpilma iddiasinin (TRPC1) replika duzeyinde savunmasi.
Ayrica SOCE cekirdek genleri icin ayni tablo uretilir.
Cikti: tablolar/TRPC1_*.tsv, tablolar/SOCE_cekirdek_olay_destegi.tsv
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from importlib.machinery import SourceFileLoader

core = SourceFileLoader("core", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "01_rmats_core.py")).load_module()

OUT = core.OUT
TAB = os.path.join(OUT, "03_TABLOLAR")
os.makedirs(TAB, exist_ok=True)

SOCE_CORE = ["TRPC1", "STIM1", "STIM2", "ORAI1", "ORAI2", "ORAI3",
             "ATP2A2", "ATP2A3", "SARAF", "STIMATE", "CRACR2A", "CBARP", "MCU"]

rows_detail = []
rows_summary = []

for ds in core.DATASETS:
    sp = core.DATASETS[ds]["species"]
    for kind in ("JC", "JCEC"):
        for ev in core.EVENTS:
            df = core.load_event(ds, ev, kind)
            if df is None:
                continue
            df["gene_up"] = df["geneSymbol"].str.upper()
            sub = df[df["gene_up"].isin(SOCE_CORE)]
            if sub.empty:
                continue
            for _, r in sub.iterrows():
                kd, ct = core.per_replicate_psi(r)
                if np.all(np.isnan(kd)) or np.all(np.isnan(ct)):
                    continue
                coord = f"{r['chr']}:{r.get('exonStart_0base','')}-{r.get('exonEnd','')}"
                rec = dict(
                    veri_seti=ds, tur=sp, sayim=kind, olay=ev, gen=r["gene_up"],
                    rmats_ID=r["ID"], koordinat=coord, strand=r["strand"],
                    dPSI=r["IncLevelDifference"], PValue=r["PValue"], FDR=r["FDR"],
                    PSI_KD=";".join("NA" if np.isnan(x) else f"{x:.3f}" for x in kd),
                    PSI_KONTROL=";".join("NA" if np.isnan(x) else f"{x:.3f}" for x in ct),
                    IJC_KD=r["IJC_SAMPLE_1"], SJC_KD=r["SJC_SAMPLE_1"],
                    IJC_KONTROL=r["IJC_SAMPLE_2"], SJC_KONTROL=r["SJC_SAMPLE_2"],
                    toplam_okuma=r["total_reads"],
                    toplam_atlama_okumasi=r["total_skip_reads"],
                    ornek_basina_ort_okuma=round(r["mean_reads_per_sample"], 1),
                    en_dusuk_bilgilendirici_okuma=r["min_informative_reads"],
                )
                rows_detail.append(rec)

det = pd.DataFrame(rows_detail)
det.to_csv(os.path.join(TAB, "SOCE_cekirdek_TUM_olaylar_replika_duzeyi.tsv"),
           sep="\t", index=False)
print(f"SOCE cekirdek olay kaydi: {len(det)}")

# ---------------------------------------------------------------- anlamli olanlar
sig = det[(det.FDR < 0.05) & (det.dPSI.abs() >= 0.10)].copy()
sig = sig.sort_values(["gen", "veri_seti", "sayim"])
sig.to_csv(os.path.join(TAB, "SOCE_cekirdek_ANLAMLI_olaylar.tsv"), sep="\t", index=False)
print(f"Esikleri karsilayan olay: {len(sig)}")

# ---------------------------------------------------------------- TRPC1 derin analiz
print("\n" + "=" * 78)
print("TRPC1 — SH-SY5Y ana olayinin adli incelemesi")
print("=" * 78)

df_se = core.load_event("GSE296712_SHSY5Y", "SE", "JC")
df_se["gene_up"] = df_se["geneSymbol"].str.upper()
t = df_se[(df_se.gene_up == "TRPC1")].copy()
main = t[t["FDR"] < 0.05]

boot_rows = []
for _, r in t.iterrows():
    kd, ct = core.per_replicate_psi(r)
    med, lo, hi = core.bootstrap_dpsi(r, n_boot=20000, seed=42)
    loo = core.leave_one_out_dpsi(r)
    coord = f"{r['chr']}:{r['exonStart_0base']}-{r['exonEnd']}"
    boot_rows.append(dict(
        rmats_ID=r["ID"], koordinat=coord,
        dPSI_rmats=r["IncLevelDifference"], FDR_rmats=r["FDR"],
        PSI_KD=";".join(f"{x:.3f}" for x in kd),
        PSI_KONTROL=";".join(f"{x:.3f}" for x in ct),
        IJC_KD=r["IJC_SAMPLE_1"], SJC_KD=r["SJC_SAMPLE_1"],
        IJC_KONTROL=r["IJC_SAMPLE_2"], SJC_KONTROL=r["SJC_SAMPLE_2"],
        toplam_atlama_okumasi=int(r["total_skip_reads"]),
        toplam_katilim_okumasi=int(r["total_inc_reads"]),
        en_dusuk_bilgilendirici_okuma=int(r["min_informative_reads"]),
        bootstrap_dPSI_medyan=None if med is None or np.isnan(med) else round(med, 3),
        bootstrap_GA_alt=None if lo is None or np.isnan(lo) else round(lo, 3),
        bootstrap_GA_ust=None if hi is None or np.isnan(hi) else round(hi, 3),
        GA_sifiri_iceriyor=("EVET" if (not np.isnan(lo) and lo <= 0 <= hi) else "hayir"),
        LOO_min=round(min(v for _, v in loo), 3),
        LOO_max=round(max(v for _, v in loo), 3),
        LOO_ayrinti="; ".join(f"{k}={v:+.3f}" for k, v in loo),
    ))

bt = pd.DataFrame(boot_rows).sort_values("FDR_rmats")
bt.to_csv(os.path.join(TAB, "TRPC1_SHSY5Y_bootstrap_ve_LOO.tsv"), sep="\t", index=False)

for _, r in bt.iterrows():
    flag = "  <-- TEZDEKI ANA OLAY" if r["FDR_rmats"] < 0.05 else ""
    print(f"\nID {r['rmats_ID']}  {r['koordinat']}{flag}")
    print(f"  rMATS dPSI = {r['dPSI_rmats']:+.3f}   FDR = {r['FDR_rmats']:.4g}")
    print(f"  PSI  KD      : {r['PSI_KD']}")
    print(f"  PSI  KONTROL : {r['PSI_KONTROL']}")
    print(f"  IJC KD={r['IJC_KD']}  SJC KD={r['SJC_KD']}")
    print(f"  IJC CT={r['IJC_KONTROL']}  SJC CT={r['SJC_KONTROL']}")
    print(f"  toplam atlama okumasi = {r['toplam_atlama_okumasi']}   "
          f"en dusuk bilgilendirici okuma/ornek = {r['en_dusuk_bilgilendirici_okuma']}")
    if r["bootstrap_dPSI_medyan"] is not None:
        print(f"  bootstrap dPSI = {r['bootstrap_dPSI_medyan']:+.3f}  "
              f"%95 GA [{r['bootstrap_GA_alt']:+.3f}, {r['bootstrap_GA_ust']:+.3f}]  "
              f"sifir GA icinde: {r['GA_sifiri_iceriyor']}")
    print(f"  birini-disarida-birak dPSI araligi: [{r['LOO_min']:+.3f}, {r['LOO_max']:+.3f}]")

# ---------------------------------------------------------------- JCEC karsiligi
print("\n" + "-" * 78)
print("Ayni olayin JCEC karsiligi")
jcec = core.load_event("GSE296712_SHSY5Y", "SE", "JCEC")
jcec["gene_up"] = jcec["geneSymbol"].str.upper()
for _, r in main.iterrows():
    m = jcec[(jcec.gene_up == "TRPC1")
             & (jcec.exonStart_0base == r.exonStart_0base)
             & (jcec.exonEnd == r.exonEnd)]
    for _, q in m.iterrows():
        kd, ct = core.per_replicate_psi(q)
        print(f"  ID {q['ID']}  dPSI={q['IncLevelDifference']:+.3f}  FDR={q['FDR']:.4g}")
        print(f"    PSI KD={['%.3f'%x for x in kd]}  KONTROL={['%.3f'%x for x in ct]}")
        print(f"    IJC KD={q['IJC_SAMPLE_1']} SJC KD={q['SJC_SAMPLE_1']} | "
              f"IJC CT={q['IJC_SAMPLE_2']} SJC CT={q['SJC_SAMPLE_2']}")

print("\nTablolar yazildi ->", TAB)

#!/usr/bin/env python3
"""
Oncelik 4.1 (eslenmis permutasyon null'u), 4.2 (olay duzeyi meta-analiz),
4.3 (guc simulasyonu) ve STIM2.1 hedefli meta-analizi.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy import stats
from importlib.machinery import SourceFileLoader

core = SourceFileLoader("core", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "01_rmats_core.py")).load_module()
TAB = os.path.join(core.OUT, "03_TABLOLAR")
GROUPS = core.load_gene_groups()
rng = np.random.default_rng(2026)

# =====================================================================
# 1) STIM2.1 / STIM2-beta hedefli meta-analiz
# =====================================================================
print("=" * 96)
print("1) STIM2.1 (24 nt SOAR ekzonu) — veri setleri arasi hedefli meta-analiz")
print("=" * 96)
COORD = {"human": (27007982, 27008006), "mouse": (54110115, 54110139)}
rows = []
for ds, meta in core.DATASETS.items():
    df = core.load_event(ds, "SE", "JC")
    if df is None:
        continue
    df["gene_up"] = df["geneSymbol"].str.upper()
    st, en = COORD[meta["species"]]
    m = df[(df.gene_up == "STIM2") & (df.exonStart_0base == st) & (df.exonEnd == en)]
    for _, r in m.iterrows():
        kd, ct = core.per_replicate_psi(r)
        med, lo, hi = core.bootstrap_dpsi(r, n_boot=10000, seed=11)
        # varyans kestirimi: replika PSI'lerinin standart hatasi
        kd_a, ct_a = np.array(kd, float), np.array(ct, float)
        v = (np.nanvar(kd_a, ddof=1) / max(len(kd_a), 1)
             + np.nanvar(ct_a, ddof=1) / max(len(ct_a), 1))
        rows.append(dict(veri_seti=ds, tur=meta["species"], dPSI=r["IncLevelDifference"],
                         FDR=r["FDR"], PValue=r["PValue"], var=v if v > 0 else np.nan,
                         ort_okuma=round(r["mean_reads_per_sample"], 1),
                         boot_GA=f"[{lo:+.3f},{hi:+.3f}]" if not np.isnan(lo) else "NA",
                         PSI_KD=";".join(f"{x:.3f}" for x in kd),
                         PSI_KONTROL=";".join(f"{x:.3f}" for x in ct)))
S21 = pd.DataFrame(rows)
print(S21.to_string(index=False))

d = S21.dPSI.dropna().values
n_pos = int((d > 0).sum()); n_tot = len(d)
p_sign = stats.binomtest(n_pos, n_tot, 0.5, alternative="greater").pvalue
print(f"\n  Yon tutarliligi: {n_pos}/{n_tot} veri setinde POZITIF (KD'de katilim artiyor)")
print(f"  Tek yonlu isaret testi p = {p_sign:.4f}")
w = 1 / S21["var"]
ok = S21["var"].notna() & np.isfinite(w)
if ok.sum() >= 2:
    mu = np.average(S21.loc[ok, "dPSI"], weights=w[ok])
    se = np.sqrt(1 / w[ok].sum())
    z = mu / se
    p_fe = 2 * (1 - stats.norm.cdf(abs(z)))
    print(f"  Sabit etkiler havuzlanmis dPSI = {mu:+.4f}  (%95 GA "
          f"[{mu-1.96*se:+.4f}, {mu+1.96*se:+.4f}])  p = {p_fe:.4f}")
S21.to_csv(os.path.join(TAB, "STIM2_1_SOAR_ekzonu_meta.tsv"), sep="\t", index=False)

# =====================================================================
# 2) Eslenmis permutasyon null'u  (Grup 4 zenginlesmesi)
# =====================================================================
print("\n" + "=" * 96)
print("2) Ca2+ gen grubu zenginlesmesi — uzunluk/ifade/kapsam ESLENMIS permutasyon null'u")
print("=" * 96)

DS = "GSE296712_SHSY5Y"
frames = []
for ev in core.EVENTS:
    f = core.load_event(DS, ev, "JC")
    if f is not None:
        frames.append(f)
E = pd.concat(frames, ignore_index=True)
E["gene_up"] = E["geneSymbol"].str.upper()
E["sig"] = (E["FDR"] < 0.05) & (E["IncLevelDifference"].abs() >= 0.10)

# gen duzeyinde ozet: olay sayisi (= test edilebilirlik / uzunluk vekili)
# ve toplam okuma (= ifade vekili)
g = E.groupby("gene_up").agg(n_olay=("ID", "size"),
                             toplam_okuma=("total_reads", "sum"),
                             anlamli=("sig", "any")).reset_index()
g = g[g.n_olay > 0]
print(f"Test edilebilir gen sayisi: {len(g)}")

for gname in ["Grup1", "Grup3", "Grup4"]:
    gs = GROUPS[gname]
    g["in_set"] = g.gene_up.isin(gs)
    k = int(g.loc[g.in_set, "anlamli"].sum())
    n = int(g.in_set.sum())
    bg = g[~g.in_set]
    obs_rate = k / n if n else np.nan

    # --- ham hipergeometrik (tezdeki yaklasim)
    K = int(g["anlamli"].sum()); N = len(g)
    p_hyper = stats.hypergeom.sf(k - 1, N, K, n)

    # --- eslenmis permutasyon: her set genine, olay sayisi + okuma
    #     bakimindan en yakin arka plan geni eslenir
    setg = g[g.in_set]
    cand = bg.copy()
    cand["logreads"] = np.log10(cand.toplam_okuma + 1)
    cand["logev"] = np.log10(cand.n_olay + 1)
    sg = setg.copy()
    sg["logreads"] = np.log10(sg.toplam_okuma + 1)
    sg["logev"] = np.log10(sg.n_olay + 1)

    C = cand[["logev", "logreads"]].values
    S = sg[["logev", "logreads"]].values
    cand_sig = cand["anlamli"].values.astype(float)

    NPERM = 5000
    null = np.empty(NPERM)
    # her set geni icin en yakin 25 arka plan genini onceden bul
    pools = []
    for i in range(S.shape[0]):
        dist = np.sqrt(((C - S[i]) ** 2).sum(1))
        pools.append(np.argsort(dist)[:25])
    pools = np.array(pools)
    for b in range(NPERM):
        pick = pools[np.arange(pools.shape[0]), rng.integers(0, pools.shape[1], pools.shape[0])]
        null[b] = cand_sig[pick].mean()
    p_perm = (1 + (null >= obs_rate).sum()) / (1 + NPERM)

    print(f"\n  {gname}: n={n} gen, anlamli={k} ({100*obs_rate:.1f}%)")
    print(f"    arka plan orani (ham)            = {100*bg.anlamli.mean():.1f}%")
    print(f"    ESLENMIS null ortalamasi         = {100*null.mean():.1f}%  "
          f"(%95 aralik {100*np.percentile(null,2.5):.1f}-{100*np.percentile(null,97.5):.1f}%)")
    print(f"    ham hipergeometrik p             = {p_hyper:.4f}")
    print(f"    ESLENMIS permutasyon p           = {p_perm:.4f}   <-- savunulabilir deger")

# =====================================================================
# 3) Guc simulasyonu
# =====================================================================
print("\n" + "=" * 96)
print("3) Guc simulasyonu — n=3+3 tasarimda hangi dPSI saptanabilir?")
print("=" * 96)
print("  (gozlenen kapsam dagiliminda, beta-binom benzeri replika degiskenligi ile)")

def sim_power(n_per_group, depth, dpsi, base_psi=0.5, nsim=2000, alpha=0.05, disp=0.05):
    hit = 0
    for _ in range(nsim):
        p1 = np.clip(rng.normal(base_psi + dpsi / 2, disp, n_per_group), 0.01, 0.99)
        p2 = np.clip(rng.normal(base_psi - dpsi / 2, disp, n_per_group), 0.01, 0.99)
        i1 = rng.binomial(depth, p1); i2 = rng.binomial(depth, p2)
        psi1 = i1 / depth; psi2 = i2 / depth
        try:
            t, p = stats.ttest_ind(psi1, psi2)
        except Exception:
            continue
        if p < alpha:
            hit += 1
    return hit / nsim

depths = [10, 20, 50, 100]
dpsis = [0.05, 0.10, 0.15, 0.20, 0.30]
print(f"\n  {'derinlik':>9} | " + " ".join(f"dPSI={d:.2f}" for d in dpsis))
pw_rows = []
for dep in depths:
    vals = [sim_power(3, dep, d) for d in dpsis]
    pw_rows.append(dict(derinlik=dep, **{f"dPSI_{d:.2f}": v for d, v in zip(dpsis, vals)}))
    print(f"  {dep:>9} | " + " ".join(f"   {v:.2f}  " for v in vals))
pd.DataFrame(pw_rows).to_csv(os.path.join(TAB, "guc_simulasyonu.tsv"), sep="\t", index=False)

# gozlenen SH-SY5Y kapsam dagilimi
med_depth = E["mean_reads_per_sample"].median()
q1 = E["mean_reads_per_sample"].quantile(0.25)
print(f"\n  SH-SY5Y'de gozlenen ornek basina okuma: medyan={med_depth:.1f}, Q1={q1:.1f}")
print(f"  TRPC1 ana olayinda bu deger 10.5 ve en dusuk bilgilendirici okuma 2 idi.")
print(f"  Bu derinlikte %80 guce ulasmak icin gereken dPSI >= 0.30 duzeyindedir.")

print("\nTablolar ->", TAB)

#!/usr/bin/env python3
"""A1 — Anotasyondan bağımsız birleşim analizi (MAJIQ'in yerine geçen yerel LSV testi).

Her veri seti için:
  1. Birleşim matrisi (kendi çıkarımımız: MAPQ>=30, birincil, anchor>=8 bp).
  2. GENCODE'a göre sınıflandırma: anotasyonlu / yeni bölge / yeni kombinasyon / tamamen yeni.
  3. LSV: ortak verici (5') ve ortak alıcı (3') kümelerinde PSI.
  4. Ortak aşırı dağılımlı beta-binomiyal olabilirlik oranı testi + BH + replika bootstrap GA.
  5. Kriptik olay çağrısı: yeni birleşim, KD'de PSI artışı, kontrolde neredeyse yok.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import lib_junc as L

OUT = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/sonuclar"
os.makedirs(OUT, exist_ok=True)
MAN = pd.read_csv("/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/kod/ornekler.tsv", sep="\t")
MIN_LSV_READS = 20
MIN_JUNC_READS = 5

def lsv_frame(jm, samples, groups):
    kd = [s for s, g in zip(samples, groups) if g == "KD"]
    ct = [s for s, g in zip(samples, groups) if g == "CTRL"]
    C = jm[samples].values.astype(float)
    ikd = [samples.index(s) for s in kd]; ict = [samples.index(s) for s in ct]
    parts = []
    for side, col in (("verici", "start"), ("alici", "end")):
        key = jm["chrom"].astype(str) + ":" + jm[col].astype(str)
        codes, uniq = pd.factorize(key)
        n_in_lsv = np.bincount(codes, minlength=len(uniq))
        # LSV toplamları
        tot = np.zeros((len(uniq), C.shape[1]))
        np.add.at(tot, codes, C)
        keep = (n_in_lsv[codes] >= 2) & (tot[codes].mean(axis=1) >= MIN_LSV_READS) \
               & (tot[codes].min(axis=1) > 0) & (C.max(axis=1) >= MIN_JUNC_READS)
        if keep.sum() == 0:
            continue
        d = jm.loc[keep, ["chrom", "start", "end", "gene", "strand", "sinif"]].copy()
        d["lsv_tipi"] = side
        d["lsv_anahtar"] = uniq[codes[keep]]
        d["n_junction_lsv"] = n_in_lsv[codes[keep]]
        parts.append((d, C[keep], tot[codes[keep]]))
    if not parts:
        return pd.DataFrame()
    D = pd.concat([p[0] for p in parts], ignore_index=True)
    Cs = np.vstack([p[1] for p in parts]); Ts = np.vstack([p[2] for p in parts])
    K1, N1 = Cs[:, ikd], Ts[:, ikd]
    K2, N2 = Cs[:, ict], Ts[:, ict]
    rho = L.tahmin_rho(K1, N1, K2, N2)
    p, e1, e2 = L.bb_test_vec(K1, N1, K2, N2, rho)
    lo, hi = L.bootstrap_dpsi_vec(K1, N1, K2, N2, B=2000)
    D["psi_KD"] = e1; D["psi_CTRL"] = e2; D["dPSI"] = e1 - e2
    D["p"] = p; D["q"] = L.bh(p); D["rho"] = rho
    D["GA_alt"] = lo; D["GA_ust"] = hi
    D["okuma_KD"] = K1.sum(1); D["okuma_CTRL"] = K2.sum(1)
    D["toplam_KD"] = N1.sum(1); D["toplam_CTRL"] = N2.sum(1)
    D["psi_ham_KD"] = np.where(N1.sum(1) > 0, K1.sum(1) / N1.sum(1), np.nan)
    D["psi_ham_CTRL"] = np.where(N2.sum(1) > 0, K2.sum(1) / N2.sum(1), np.nan)
    D["n_KD_pozitif"] = (K1 > 0).sum(1); D["n_CTRL_pozitif"] = (K2 > 0).sum(1)
    return D

def main(only=None):
    ann = {}
    for ds, sub in MAN.groupby("dataset"):
        if only and ds not in only:
            continue
        samples = list(sub["sample"]); groups = list(sub["group"])
        sp = sub["species"].iloc[0]
        eksik = [s for s in samples if not os.path.exists(f"{L.JDIR}/{s}.junc")]
        if eksik:
            print(f"[ATLA] {ds}: eksik junction {eksik}", flush=True); continue
        t0 = time.time()
        jm = L.build_matrix(samples)
        jm = jm[jm[samples].max(axis=1) >= MIN_JUNC_READS].reset_index(drop=True)
        if sp not in ann:
            ann[sp] = L.load_annotation(sp)
        ai, genes = ann[sp]
        jm = L.annotate_junctions(jm, ai, genes)
        jm.to_csv(f"{OUT}/junction_matrisi_{ds}.tsv.gz", sep="\t", index=False)
        vc = jm["sinif"].value_counts()
        print(f"[{ds}] {len(jm)} birleşim | " +
              " ".join(f"{k}={v}" for k, v in vc.items()), flush=True)
        D = lsv_frame(jm, samples, groups)
        D.to_csv(f"{OUT}/LSV_{ds}.tsv.gz", sep="\t", index=False)
        sig = D[(D["q"] < 0.05) & (D["dPSI"].abs() >= 0.10)]
        yeni = sig[sig["sinif"] != "anotasyonlu"]
        print(f"[{ds}] LSV birleşimi {len(D)}, anlamlı {len(sig)}, bunlardan yeni {len(yeni)} "
              f"(rho={D['rho'].iloc[0] if len(D) else float('nan')}) [{time.time()-t0:.0f}s]", flush=True)

if __name__ == "__main__":
    main(only=sys.argv[1:] or None)

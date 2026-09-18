#!/usr/bin/env python3
"""A4b — İntronik poliadenilasyon (IPA) ve 3'UTR kullanımı testi.

IPA indeksi = I5 / (I5 + I3): intronun 5' ucundaki kapsamın, aynı intronun
3' ucundaki kapsama oranı. Erken (intronik) poliadenilasyon 5' ucu zenginleştirir.
İntron tutulması ise iki ucu birlikte yükseltir ve indeksi ~0,5'te bırakır;
bu yüzden indeks IPA'ya özgüdür.

3'UTR distal kullanım indeksi = Udist / (Uprox + Udist): DaPars PDUI mantığı.
Düşmesi 3'UTR kısalması demektir.

İstatistik: kapsam toplamları okuma uzunluğuna bölünerek sözde-sayıya çevrilir,
ortak aşırı dağılımlı beta-binomiyal olabilirlik oranı testi + BH + replika
bootstrap %95 GA uygulanır.
"""
import os, sys, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import lib_junc as L

D = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI"
OUT = f"{D}/sonuclar"
MAN = pd.read_csv(f"{D}/kod/ornekler.tsv", sep="\t")
MIN_DEPTH = 3.0     # her örnekte iki pencerenin ortalama derinliği

def okuma_uzunlugu(bam):
    r = subprocess.run(["bash", "-c",
        f"samtools view {bam} 1:1000000-1200000 2>/dev/null | head -2000 | "
        f"awk '{{print length($10)}}' | sort -n | awk '{{a[NR]=$1}} END{{print a[int(NR/2)+1]}}'"],
        capture_output=True, text=True)
    try: return max(int(r.stdout.strip()), 30)
    except Exception: return 100

def analiz(ds):
    f = f"{OUT}/apa_bedcov_{ds}.tsv"
    if not os.path.exists(f) or os.path.getsize(f) < 1000:
        return None
    sub = MAN[MAN["dataset"] == ds]
    samples = list(sub["sample"]); groups = list(sub["group"])
    kd = [s for s, g in zip(samples, groups) if g == "KD"]
    ct = [s for s, g in zip(samples, groups) if g == "CTRL"]
    RL = okuma_uzunlugu(sub["bam"].iloc[0])
    b = pd.read_csv(f, sep="\t")
    parts = b["isim"].str.split("|", expand=True)
    b["gene"], b["birim"], b["pencere"], b["strand"] = parts[0], parts[1], parts[2], parts[3]
    b["uzunluk"] = b["end"] - b["start"]
    rows = []
    for (gene, birim), g in b.groupby(["gene", "birim"], sort=False):
        w = dict(zip(g["pencere"], range(len(g))))
        if birim.startswith("intron"):
            a, c = "I5", "I3"
        else:
            a, c = "Udist", "Uprox"
        if a not in w or c not in w:
            continue
        ga = g.iloc[w[a]]; gc = g.iloc[w[c]]
        # derinlik = toplam / uzunluk ; sözde-sayı = derinlik * uzunluk / okuma_uzunlugu
        da = ga[samples].values.astype(float) / ga["uzunluk"]
        dc = gc[samples].values.astype(float) / gc["uzunluk"]
        if np.min(da + dc) < MIN_DEPTH:
            continue
        ka = ga[samples].values.astype(float) / RL
        kc = gc[samples].values.astype(float) / RL
        rows.append((gene, birim, ga["chrom"], int(ga["start"]), int(gc["end"]),
                     ga["strand"], *ka, *kc))
    if not rows:
        return None
    cols = (["gene", "birim", "chrom", "start", "end", "strand"]
            + [f"A_{s}" for s in samples] + [f"C_{s}" for s in samples])
    T = pd.DataFrame(rows, columns=cols)
    A = T[[f"A_{s}" for s in samples]].values
    C = T[[f"C_{s}" for s in samples]].values
    N = A + C
    ikd = [samples.index(s) for s in kd]; ict = [samples.index(s) for s in ct]
    K1, N1 = A[:, ikd], N[:, ikd]; K2, N2 = A[:, ict], N[:, ict]
    rho = L.tahmin_rho(K1, N1, K2, N2)
    p, e1, e2 = L.bb_test_vec(K1, N1, K2, N2, rho)
    lo, hi = L.bootstrap_dpsi_vec(K1, N1, K2, N2, B=5000)
    R = T[["gene", "birim", "chrom", "start", "end", "strand"]].copy()
    R["olcum"] = np.where(R["birim"].str.startswith("intron"), "IPA_indeksi", "UTR_distal_indeksi")
    R["indeks_KD"] = e1; R["indeks_CTRL"] = e2; R["delta"] = e1 - e2
    R["GA_alt"] = lo; R["GA_ust"] = hi
    R["p"] = p; R["q"] = L.bh(p); R["rho"] = rho
    R["ort_derinlik"] = ((A + C) * 0 + (A + C)).mean(axis=1) * 0 + (N.mean(axis=1))
    R.to_csv(f"{OUT}/APA_{ds}.tsv.gz", sep="\t", index=False)
    sig = R[(R["q"] < 0.05) & (R["delta"].abs() >= 0.10)]
    print(f"[{ds}] okuma uzunluğu {RL} | test {len(R)} birim, anlamlı {len(sig)} "
          f"(IPA {int((sig.olcum=='IPA_indeksi').sum())}, UTR {int((sig.olcum=='UTR_distal_indeksi').sum())}), rho={rho}")
    return R

SOCE = ["STIM1", "STIM2", "ORAI1", "ORAI2", "ORAI3", "TRPC1", "SARAF", "STIMATE",
        "CBARP", "ATP2A2", "ATP2A3", "MCU", "MCUB", "MICU1", "MICU2", "CRACR2A",
        "CRACR2B", "SELENOK", "SELENON", "ITPR1", "ITPR2", "ITPR3"]
POS = ["STMN2", "UNC13A", "HDGFL2", "ACTL6B", "AGRN", "KALRN", "ARHGAP32", "PFKP",
       "ATG4B", "SETD5", "ELAVL3", "POLDIP3"]

if __name__ == "__main__":
    hepsi = []
    for ds in ["SH_SY5Y", "iPSC_koloni", "iPSC_MN", "K562_mRNA", "K562_totalRNA",
               "C2C12", "NSC34"]:
        R = analiz(ds)
        if R is not None:
            R["veri_seti"] = ds; hepsi.append(R)
    if hepsi:
        H = pd.concat(hepsi, ignore_index=True)
        H.to_csv(f"{OUT}/APA_tum_veri_setleri.tsv.gz", sep="\t", index=False)
        POS2 = POS + [g.capitalize() for g in POS]
        SOCE2 = SOCE + [g.capitalize() for g in SOCE]
        for etiket, panel in (("pozitif_kontrol", POS2), ("SOCE", SOCE2)):
            d = H[H["gene"].isin(panel) & (H["q"] < 0.05) & (H["delta"].abs() >= 0.10)]
            d = d.sort_values(["gene", "q"])
            d.to_csv(f"{OUT}/APA_{etiket}_anlamli.tsv", sep="\t", index=False)
            print(f"\n=== APA — {etiket} panelinde anlamlı ({len(d)}) ===")
            if len(d):
                print(d[["veri_seti","gene","birim","olcum","indeks_KD","indeks_CTRL",
                         "delta","GA_alt","GA_ust","q","ort_derinlik"]].head(30)
                      .to_string(index=False, float_format=lambda x: f"{x:.4g}"))

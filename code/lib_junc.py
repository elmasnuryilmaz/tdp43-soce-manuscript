"""Junction-düzeyi analiz kütüphanesi.

Anotasyondan bağımsız birleşim (junction) analizi: yükleme, anotasyon,
LSV (yerel kırpılma varyasyonu) kurma, beta-binomiyal diferansiyel test.
Tüm koordinatlar 1-tabanlı, intron kapsayıcı (start = ilk intronik baz).
"""
import gzip, os, sys, math
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.special import betaln, gammaln
from scipy.optimize import minimize_scalar, minimize
from scipy import stats

REF = "/Volumes/10TBElmas/thesis_addendum_2026/ref"
JDIR = "/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI/junctions"
RNG = np.random.default_rng(42)

# ---------------------------------------------------------------- yükleme
def load_junc(sample):
    """Kendi çıkardığımız .junc dosyası: chrom start end strand count."""
    p = os.path.join(JDIR, f"{sample}.junc")
    df = pd.read_csv(p, sep="\t", header=None,
                     names=["chrom", "start", "end", "strand", "count"],
                     dtype={"chrom": str})
    # aynı intron farklı XS etiketleriyle gelebilir; şeritten bağımsız topla
    df = df.groupby(["chrom", "start", "end"], as_index=False)["count"].sum()
    return df

def load_junc_regtools(path):
    """regtools BED12 -> 1-tabanlı intron koordinatı."""
    df = pd.read_csv(path, sep="\t", header=None, dtype={0: str})
    bs = df[10].str.split(",", expand=True).astype(int)
    bt = df[11].str.split(",", expand=True).astype(int)
    out = pd.DataFrame({"chrom": df[0].astype(str),
                        "start": df[1] + bs[0] + 1,
                        "end": df[1] + bt[1],
                        "count": df[4]})
    return out.groupby(["chrom", "start", "end"], as_index=False)["count"].sum()

def build_matrix(samples, loader=load_junc):
    """Örnekler arası birleşim x örnek sayım matrisi (uzun -> geniş)."""
    frames = []
    for s in samples:
        d = loader(s) if loader is load_junc else loader(s)
        d = d.rename(columns={"count": s})
        frames.append(d.set_index(["chrom", "start", "end"])[s])
    mat = pd.concat(frames, axis=1).fillna(0).astype(int)
    return mat.reset_index()

# ------------------------------------------------------------- anotasyon
def load_annotation(species="human"):
    ai = pd.read_csv(f"{REF}/annotated_introns_{species}.bed", sep="\t", header=None,
                     names=["chrom", "start0", "end", "info", "n", "strand"],
                     dtype={"chrom": str})
    ai["start"] = ai["start0"] + 1
    ai["gene"] = ai["info"].str.split("|").str[0]
    genes = pd.read_csv(f"{REF}/genes_{species}.bed", sep="\t", header=None,
                        names=["chrom", "start0", "end", "info", "n", "strand"],
                        dtype={"chrom": str})
    genes["start"] = genes["start0"] + 1
    parts = genes["info"].str.split("|", expand=True)
    genes["gene"] = parts[0]; genes["gid"] = parts[1]; genes["biotype"] = parts[2]
    return ai, genes

def annotate_junctions(jm, ai, genes):
    """Her birleşime gen ata ve verici/alıcı bölgelerin bilinip bilinmediğini işaretle."""
    known_intron = set(zip(ai["chrom"], ai["start"], ai["end"]))
    known_don = defaultdict(set); known_acc = defaultdict(set)
    for c, s, e in zip(ai["chrom"], ai["start"], ai["end"]):
        known_don[c].add(s); known_acc[c].add(e)

    # gen ataması: kromozoma göre sıralı aralık araması
    gidx = defaultdict(list)
    for c, s, e, g, st in zip(genes["chrom"], genes["start"], genes["end"],
                              genes["gene"], genes["strand"]):
        gidx[c].append((s, e, g, st))
    for c in gidx:
        gidx[c].sort()
    starts = {c: np.array([x[0] for x in v]) for c, v in gidx.items()}

    gene_col, strand_col, don_col, acc_col, ann_col = [], [], [], [], []
    for c, s, e in zip(jm["chrom"], jm["start"], jm["end"]):
        # gen
        g, st = ".", "."
        if c in gidx:
            i = np.searchsorted(starts[c], s, side="right")
            best = None
            for j in range(max(0, i - 400), min(len(gidx[c]), i + 2)):
                gs, ge, gn, gst = gidx[c][j]
                if gs <= s and e <= ge:
                    ln = ge - gs
                    if best is None or ln < best[0]:
                        best = (ln, gn, gst)
            if best:
                g, st = best[1], best[2]
        gene_col.append(g); strand_col.append(st)
        don_col.append(s in known_don.get(c, ()))
        acc_col.append(e in known_acc.get(c, ()))
        ann_col.append((c, s, e) in known_intron)
    jm = jm.copy()
    jm["gene"] = gene_col; jm["strand"] = strand_col
    jm["donor_known"] = don_col; jm["acceptor_known"] = acc_col
    jm["annotated"] = ann_col
    jm["sinif"] = np.where(jm["annotated"], "anotasyonlu",
                  np.where(jm["donor_known"] & jm["acceptor_known"], "yeni_kombinasyon",
                  np.where(jm["donor_known"] | jm["acceptor_known"], "yeni_bolge", "tamamen_yeni")))
    return jm

# -------------------------------------------------- beta-binomiyal test
def _bb_nll(counts, totals, p, rho):
    p = min(max(p, 1e-9), 1 - 1e-9); rho = min(max(rho, 1e-9), 1 - 1e-9)
    s = (1 - rho) / rho
    a, b = p * s, (1 - p) * s
    k = np.asarray(counts, float); n = np.asarray(totals, float)
    return -np.sum(betaln(k + a, n - k + b) - betaln(a, b))

def bb_test(k1, n1, k2, n2, rho_grid=(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2)):
    """İki grup arasında beta-binomiyal olabilirlik oranı testi (ortak rho).

    k,n replika düzeyi dizileri. Döndürür: (p_deger, psi1, psi2, rho).
    """
    k1 = np.asarray(k1, float); n1 = np.asarray(n1, float)
    k2 = np.asarray(k2, float); n2 = np.asarray(n2, float)
    if n1.sum() == 0 or n2.sum() == 0:
        return np.nan, np.nan, np.nan, np.nan
    kk = np.concatenate([k1, k2]); nn = np.concatenate([n1, n2])
    best = None
    for rho in rho_grid:
        # H1: iki ayrı p
        p1 = k1.sum() / max(n1.sum(), 1); p2 = k2.sum() / max(n2.sum(), 1)
        f1 = minimize(lambda x: _bb_nll(k1, n1, x[0], rho), [max(min(p1,.99),.01)],
                      bounds=[(1e-6, 1 - 1e-6)], method="L-BFGS-B")
        f2 = minimize(lambda x: _bb_nll(k2, n2, x[0], rho), [max(min(p2,.99),.01)],
                      bounds=[(1e-6, 1 - 1e-6)], method="L-BFGS-B")
        p0 = kk.sum() / max(nn.sum(), 1)
        f0 = minimize(lambda x: _bb_nll(kk, nn, x[0], rho), [max(min(p0,.99),.01)],
                      bounds=[(1e-6, 1 - 1e-6)], method="L-BFGS-B")
        ll1 = -(f1.fun + f2.fun); ll0 = -f0.fun
        if best is None or ll1 > best[0]:
            best = (ll1, ll0, f1.x[0], f2.x[0], rho)
    ll1, ll0, e1, e2, rho = best
    lr = 2 * (ll1 - ll0)
    p = stats.chi2.sf(max(lr, 0), 1)
    return p, e1, e2, rho

def bh(pvals):
    p = np.asarray(pvals, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pp = p[ok]; n = len(pp)
    if n == 0:
        return q
    order = np.argsort(pp)
    ranked = pp[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q

def bootstrap_dpsi(k1, n1, k2, n2, B=10000, rng=None):
    """Replika düzeyi bootstrap ile ΔPSI %95 güven aralığı."""
    rng = rng or RNG
    k1 = np.asarray(k1, float); n1 = np.asarray(n1, float)
    k2 = np.asarray(k2, float); n2 = np.asarray(n2, float)
    a = len(k1); b = len(k2)
    out = np.empty(B)
    for i in range(B):
        i1 = rng.integers(0, a, a); i2 = rng.integers(0, b, b)
        t1 = n1[i1].sum(); t2 = n2[i2].sum()
        if t1 == 0 or t2 == 0:
            out[i] = np.nan; continue
        out[i] = k1[i1].sum() / t1 - k2[i2].sum() / t2
    out = out[~np.isnan(out)]
    if len(out) < 100:
        return np.nan, np.nan
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))

# ---------------------------------------------- vektörleştirilmiş beta-binomiyal
def _bb_ll_grid(K, N, pgrid, rho):
    """K,N: (J, R) sayım/toplam. Dönüş: (J, G) her p ızgara noktası için log-olabilirlik."""
    s = (1.0 - rho) / rho
    a = pgrid * s                     # (G,)
    b = (1.0 - pgrid) * s
    Kx = K[:, :, None]; Nx = N[:, :, None]
    ll = (betaln(Kx + a[None, None, :], Nx - Kx + b[None, None, :])
          - betaln(a[None, None, :], b[None, None, :]))
    return ll.sum(axis=1)

def bb_test_vec(K1, N1, K2, N2, rho, grid=801, chunk=1500):
    """Gruplar arası beta-binomiyal olabilirlik oranı testi, ızgara MLE ile.

    K1/N1: (J, R1), K2/N2: (J, R2). Ortak rho. Dönüş: p, psi1, psi2 dizileri.
    """
    J = K1.shape[0]
    pg = np.linspace(1e-4, 1 - 1e-4, grid)
    p_out = np.empty(J); e1 = np.empty(J); e2 = np.empty(J)
    Kc = np.concatenate([K1, K2], axis=1); Nc = np.concatenate([N1, N2], axis=1)
    for i in range(0, J, chunk):
        sl = slice(i, min(i + chunk, J))
        l1 = _bb_ll_grid(K1[sl], N1[sl], pg, rho)
        l2 = _bb_ll_grid(K2[sl], N2[sl], pg, rho)
        l0 = _bb_ll_grid(Kc[sl], Nc[sl], pg, rho)
        i1 = l1.argmax(1); i2 = l2.argmax(1)
        lr = 2 * (l1.max(1) + l2.max(1) - l0.max(1))
        p_out[sl] = stats.chi2.sf(np.maximum(lr, 0), 1)
        e1[sl] = pg[i1]; e2[sl] = pg[i2]
    return p_out, e1, e2

def tahmin_rho(K1, N1, K2, N2, adaylar=(0.0005, 0.001, 0.0025, 0.005, 0.01, 0.02,
                                        0.05, 0.1, 0.2), n_ornek=4000, rng=None):
    """Veri seti geneli ortak aşırı dağılım (rho) profil olabilirlikle seçilir."""
    rng = rng or RNG
    J = K1.shape[0]
    idx = rng.choice(J, size=min(n_ornek, J), replace=False)
    pg = np.linspace(1e-4, 1 - 1e-4, 401)
    en_iyi = (None, -np.inf)
    for r in adaylar:
        ll = float(np.nansum(_bb_ll_grid(K1[idx], N1[idx], pg, r).max(1))
                   + np.nansum(_bb_ll_grid(K2[idx], N2[idx], pg, r).max(1)))
        if np.isfinite(ll) and ll > en_iyi[1]:
            en_iyi = (r, ll)
    return en_iyi[0] if en_iyi[0] is not None else 0.01

def bootstrap_dpsi_vec(K1, N1, K2, N2, B=2000, rng=None):
    """Tüm birleşimler için eşzamanlı replika bootstrap ΔPSI %95 GA."""
    rng = rng or RNG
    J, R1 = K1.shape; R2 = K2.shape[1]
    lo = np.empty(J); hi = np.empty(J)
    i1 = rng.integers(0, R1, size=(B, R1)); i2 = rng.integers(0, R2, size=(B, R2))
    for i in range(0, J, 500):
        sl = slice(i, min(i + 500, J))
        k1 = K1[sl][:, i1].sum(axis=2); n1 = N1[sl][:, i1].sum(axis=2)
        k2 = K2[sl][:, i2].sum(axis=2); n2 = N2[sl][:, i2].sum(axis=2)
        with np.errstate(invalid="ignore", divide="ignore"):
            d = np.where(n1 > 0, k1 / n1, np.nan) - np.where(n2 > 0, k2 / n2, np.nan)
        lo[sl] = np.nanpercentile(d, 2.5, axis=1)
        hi[sl] = np.nanpercentile(d, 97.5, axis=1)
    return lo, hi

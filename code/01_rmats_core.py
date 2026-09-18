#!/usr/bin/env python3
"""
Ortak rMATS yardimcilari.
b1 = KD, b2 = KONTROL  (GSE27394 b1=kd_rep*, b2=ctrl_rep* ile dogrulandi)
IncLevelDifference = IncLevel1(KD) - IncLevel2(Kontrol)
"""
import os, glob
import numpy as np
import pandas as pd

RM = "/Users/elmas/Desktop/TEZ/output/rmats_libtype_corrected_2026-08-21"
OUT = "/Users/elmas/Desktop/MAKALE"

DATASETS = {
    "GSE296712_SHSY5Y":   dict(label="SH-SY5Y (birincil)", species="human", n1=3, n2=3),
    "GSE230647_iPSC_koloni": dict(label="iPSC koloni",     species="human", n1=4, n2=4),
    "GSE77702_iPSC_MN":   dict(label="iPSC-MN",            species="human", n1=2, n2=2),
    "GSE27394_mouse_SE":  dict(label="Fare striatum SE",   species="mouse", n1=4, n2=4),
    "Mouse_PE_C2C12":     dict(label="C2C12",              species="mouse", n1=3, n2=3),
    "Mouse_PE_NSC34":     dict(label="NSC34",              species="mouse", n1=3, n2=3),
}
EVENTS = ["SE", "A5SS", "A3SS", "MXE", "RI"]


def _ints(s):
    if pd.isna(s):
        return []
    out = []
    for x in str(s).split(","):
        x = x.strip()
        out.append(np.nan if x in ("", "NA") else float(x))
    return out


def psi_from_counts(ijc, sjc, inc_len, skip_len):
    """rMATS PSI = (IJC/IncFormLen) / (IJC/IncFormLen + SJC/SkipFormLen)."""
    if inc_len in (0, None) or skip_len in (0, None):
        return np.nan
    a = ijc / inc_len
    b = sjc / skip_len
    return np.nan if (a + b) == 0 else a / (a + b)


def load_event(ds, ev, kind="JC"):
    """Bir veri seti + olay turu icin rMATS tablosunu okur, sayimlari acar."""
    fp = os.path.join(RM, ds, f"{ev}.MATS.{kind}.txt")
    if not os.path.exists(fp):
        return None
    df = pd.read_csv(fp, sep="\t", low_memory=False)
    for c in ("geneSymbol", "GeneID", "chr", "strand"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.strip('"')
    df["eventType"] = ev
    df["dataset"] = ds
    # sayim listeleri
    for col in ["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2"]:
        df[col + "_list"] = df[col].map(_ints)
    df["n_kd"] = df["IJC_SAMPLE_1_list"].map(len)
    df["n_ctrl"] = df["IJC_SAMPLE_2_list"].map(len)
    # toplam okuma destegi
    df["sum_IJC_KD"] = df["IJC_SAMPLE_1_list"].map(lambda v: np.nansum(v) if v else 0)
    df["sum_SJC_KD"] = df["SJC_SAMPLE_1_list"].map(lambda v: np.nansum(v) if v else 0)
    df["sum_IJC_CTRL"] = df["IJC_SAMPLE_2_list"].map(lambda v: np.nansum(v) if v else 0)
    df["sum_SJC_CTRL"] = df["SJC_SAMPLE_2_list"].map(lambda v: np.nansum(v) if v else 0)
    df["total_reads"] = (df["sum_IJC_KD"] + df["sum_SJC_KD"]
                         + df["sum_IJC_CTRL"] + df["sum_SJC_CTRL"])
    n_tot = df["n_kd"] + df["n_ctrl"]
    df["mean_reads_per_sample"] = df["total_reads"] / n_tot.replace(0, np.nan)
    # ornek basina en dusuk bilgilendirici okuma
    def _min_inf(r):
        vals = []
        for i in range(len(r["IJC_SAMPLE_1_list"])):
            vals.append(np.nansum([r["IJC_SAMPLE_1_list"][i], r["SJC_SAMPLE_1_list"][i]]))
        for i in range(len(r["IJC_SAMPLE_2_list"])):
            vals.append(np.nansum([r["IJC_SAMPLE_2_list"][i], r["SJC_SAMPLE_2_list"][i]]))
        return min(vals) if vals else 0
    df["min_informative_reads"] = df.apply(_min_inf, axis=1)
    # atlama formu toplam destegi (SE icin kritik)
    df["total_skip_reads"] = df["sum_SJC_KD"] + df["sum_SJC_CTRL"]
    df["total_inc_reads"] = df["sum_IJC_KD"] + df["sum_IJC_CTRL"]
    df["event_key"] = (df["dataset"] + "|" + df["eventType"] + "|" + df["chr"] + ":"
                       + df.get("exonStart_0base", pd.Series([""] * len(df))).astype(str)
                       + "-" + df.get("exonEnd", pd.Series([""] * len(df))).astype(str)
                       + "|" + df["ID"].astype(str))
    return df


def load_all(kind="JC", events=None, datasets=None):
    events = events or EVENTS
    datasets = datasets or list(DATASETS)
    frames = []
    for ds in datasets:
        for ev in events:
            d = load_event(ds, ev, kind)
            if d is not None:
                frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else None


def bh_fdr(p):
    """Benjamini-Hochberg. NaN'lar korunur."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    res = np.empty(n)
    res[order] = q
    out[ok] = res
    return out


def per_replicate_psi(row):
    """Her replika icin PSI dondurur: (kd_list, ctrl_list)."""
    il, sl = row["IncFormLen"], row["SkipFormLen"]
    kd = [psi_from_counts(i, s, il, sl)
          for i, s in zip(row["IJC_SAMPLE_1_list"], row["SJC_SAMPLE_1_list"])]
    ct = [psi_from_counts(i, s, il, sl)
          for i, s in zip(row["IJC_SAMPLE_2_list"], row["SJC_SAMPLE_2_list"])]
    return kd, ct


def bootstrap_dpsi(row, n_boot=20000, seed=0):
    """
    Replika duzeyinde bootstrap: her grupta replikalar yerine koyarak yeniden
    orneklenir, ardindan havuzlanmis sayimlardan dPSI hesaplanir.
    Replika sayisi az oldugunda GA genis cikar - amac da budur.
    """
    rng = np.random.default_rng(seed)
    il, sl = row["IncFormLen"], row["SkipFormLen"]
    kd_i = np.array(row["IJC_SAMPLE_1_list"], dtype=float)
    kd_s = np.array(row["SJC_SAMPLE_1_list"], dtype=float)
    ct_i = np.array(row["IJC_SAMPLE_2_list"], dtype=float)
    ct_s = np.array(row["SJC_SAMPLE_2_list"], dtype=float)
    n1, n2 = len(kd_i), len(ct_i)
    if n1 == 0 or n2 == 0:
        return (np.nan, np.nan, np.nan)
    d = np.empty(n_boot)
    for b in range(n_boot):
        a = rng.integers(0, n1, n1)
        c = rng.integers(0, n2, n2)
        p1 = psi_from_counts(kd_i[a].sum(), kd_s[a].sum(), il, sl)
        p2 = psi_from_counts(ct_i[c].sum(), ct_s[c].sum(), il, sl)
        d[b] = (p1 - p2) if (not np.isnan(p1) and not np.isnan(p2)) else np.nan
    d = d[~np.isnan(d)]
    if d.size < 100:
        return (np.nan, np.nan, np.nan)
    return (float(np.median(d)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)))


def leave_one_out_dpsi(row):
    """Her replikayi tek tek disarida birakip dPSI'nin nasil degistigini olcer."""
    il, sl = row["IncFormLen"], row["SkipFormLen"]
    kd_i = np.array(row["IJC_SAMPLE_1_list"], dtype=float)
    kd_s = np.array(row["SJC_SAMPLE_1_list"], dtype=float)
    ct_i = np.array(row["IJC_SAMPLE_2_list"], dtype=float)
    ct_s = np.array(row["SJC_SAMPLE_2_list"], dtype=float)
    res = []
    for i in range(len(kd_i)):
        m = np.ones(len(kd_i), bool); m[i] = False
        p1 = np.nanmean([psi_from_counts(a, b, il, sl) for a, b in zip(kd_i[m], kd_s[m])])
        p2 = np.nanmean([psi_from_counts(a, b, il, sl) for a, b in zip(ct_i, ct_s)])
        res.append((f"KD_rep{i+1}_cikarildi", p1 - p2))
    for j in range(len(ct_i)):
        m = np.ones(len(ct_i), bool); m[j] = False
        p1 = np.nanmean([psi_from_counts(a, b, il, sl) for a, b in zip(kd_i, kd_s)])
        p2 = np.nanmean([psi_from_counts(a, b, il, sl) for a, b in zip(ct_i[m], ct_s[m])])
        res.append((f"CTRL_rep{j+1}_cikarildi", p1 - p2))
    return res


def load_gene_groups():
    """Ca2+ gen gruplarini kurar (Grup 1-4, kumulatif)."""
    fp = ("/Users/elmas/Desktop/TEZ/07_ANALIZ_PAKETLERI/"
          "TDP43_SOCE_TRP_ANLAMLI_ANALIZLER_2026-04-25/"
          "90_IKINCIL_full_TDP43KD_DEG_significant_only/calcium_gene_lists_used.tsv")
    t = pd.read_csv(fp, sep="\t")
    t["gene"] = t["gene"].astype(str).str.upper().str.strip()
    by = {k: set(v["gene"]) for k, v in t.groupby("list_name")}
    cur = by.get("Curated_258", set())
    soce_cats = {"SOCE_core", "SOCE_modulators", "TRP_channels"}
    g1 = set(t[(t.list_name == "Curated_258") & (t.category.isin(soce_cats))]["gene"])
    g1 |= by.get("SOCE_core_manual", set())
    cats2 = soce_cats | {"ER_channels", "PLC_IP3", "PM_channels_other", "VGCC",
                         "Endolysosomal", "Ca_pumps"}
    g2 = g1 | set(t[(t.list_name == "Curated_258") & (t.category.isin(cats2))]["gene"])
    g3 = g2 | cur
    g4 = g3 | by.get("KEGG_calcium_cAMP", set()) | by.get("Calcium_genes_analysis", set())
    return {"Grup1": g1, "Grup2": g2, "Grup3": g3, "Grup4": g4}


if __name__ == "__main__":
    g = load_gene_groups()
    for k, v in g.items():
        print(f"{k}: {len(v)} gen")

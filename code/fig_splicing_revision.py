"""Figures for the 25 September 2026 revision.

Figure 4                 splicing evidence in SOCE-related genes (main text)
Supplementary Figure S3  CBARP locus: BAM-based Sashimi plot and rMATS events
Figure 6                 working model (arrows now separate measurements from hypotheses)

Inputs are package files only: tables/Table3_robust_SOCE_splicing_events.csv,
supplementary/S15_STIM2.1_exon_six_datasets.csv, supplementary/S3_rMATS_significant_events.csv.gz
and source_data/CBARP_locus/* (written by code/cbarp_bam_extract.py from the alignments).
"""
from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MPath
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "figures" / "main"
SUPP = ROOT / "figures" / "supplementary"
LOCUS = ROOT / "source_data" / "CBARP_locus"
BLUE, ORANGE, INK, MUTED, GRID = "#0072B2", "#E69F00", "#222222", "#6B6B6B", "#D9D9D9"
GROUP = {"Control": BLUE, "TDP-43 KD": ORANGE}

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.titlesize": 8.5, "axes.labelsize": 7.5,
    "axes.linewidth": 0.6, "xtick.labelsize": 6.8, "ytick.labelsize": 6.5,
    "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.facecolor": "white",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.bf": "Arial:bold",
})


def save(fig, folder, stem, pdf=False):
    for ext in ("png", "svg") + (("pdf",) if pdf else ()):
        path = folder / f"{stem}.{ext}"
        fig.savefig(path, dpi=400)
        if ext == "svg":
            path.write_text("\n".join(l.rstrip() for l in path.read_text().splitlines()) + "\n")
    plt.close(fig)


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", length=2.5, width=0.6)


def group_points(ax, x0, values, group, s=24):
    xs = x0 + np.linspace(-0.12, 0.12, len(values))
    ax.scatter(xs, values, s=s, zorder=3, lw=0.6,
               facecolor="white" if group == "Control" else ORANGE,
               edgecolor=BLUE if group == "Control" else INK)
    ax.plot([x0 - 0.25, x0 + 0.25], [np.mean(values)] * 2, color=GROUP[group], lw=1.5, zorder=2)


# ------------------------------------------------------------------ CBARP locus data
junc = pd.read_csv(LOCUS / "CBARP_junction_counts_per_library.tsv", sep="\t")
covj = json.loads((LOCUS / "CBARP_coverage_per_library.json").read_text())
C0 = covj["meta"]["start"]
samples = pd.DataFrame(covj["meta"]["samples"])
sh = junc[junc.dataset == "SH-SY5Y"].groupby(["group", "start", "end"])["count"].sum()
# Same reads as the manuscript's local-splicing-variation test (Methods 2.5, Section 3.5).
assert sh[("Control", 1235343, 1235500)] == 26 and sh[("TDP-43 KD", 1235343, 1235500)] == 33
assert junc[(junc.dataset == "SH-SY5Y") & (junc.end == 1235500)].groupby("group")["count"].sum().to_dict() \
    == {"Control": 148, "TDP-43 KD": 36}


def junction_b_share():
    rows = []
    for (d, g, s), sub in junc.groupby(["dataset", "group", "sample"]):
        total = sub.loc[sub.end == 1235500, "count"].sum()
        b = sub.loc[(sub.start == 1235343) & (sub.end == 1235500), "count"].sum()
        rows.append(dict(dataset=d, group=g, sample=s, junction_b=int(b), exon4_donor_total=int(total),
                         share_b=b / total))
    t = pd.DataFrame(rows)
    t.to_csv(LOCUS / "CBARP_junction_b_share_per_library.tsv", sep="\t", index=False, float_format="%.4f")
    return t


# ------------------------------------------------------------------ Figure 4
def figure4():
    t3 = pd.read_csv(ROOT / "tables" / "Table3_robust_SOCE_splicing_events.csv").set_index("gene")
    s15 = pd.read_csv(ROOT / "supplementary" / "S15_STIM2.1_exon_six_datasets.csv")
    share = junction_b_share()

    fig = plt.figure(figsize=(7.1, 5.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.08], hspace=0.55)
    top = gs[0].subgridspec(1, 4, wspace=0.55)
    bot = gs[1].subgridspec(1, 3, width_ratios=[0.30, 1.25, 1], wspace=0.45)

    for i, gene in enumerate(["STIMATE", "ORAI3", "STIM2", "STIM1"]):
        ax = fig.add_subplot(top[i])
        r = t3.loc[gene]
        c = np.array([float(x) for x in r.PSI_control.split(";")])
        k = np.array([float(x) for x in r.PSI_knockdown.split(";")])
        group_points(ax, 0, c, "Control")
        group_points(ax, 1, k, "TDP-43 KD")
        top_y = max(c.max(), k.max())
        ax.set_ylim(-0.03, top_y * 1.45 if top_y > 0.3 else 0.5)
        ax.set_xlim(-0.5, 1.5)
        ax.set_xticks([0, 1], ["Control", "KD"])
        ax.set_ylabel("PSI" if i == 0 else "")
        ax.set_title(("A  " if i == 0 else "") + gene, loc="left", weight="bold")
        label = f"ΔPSI {r.delta_PSI:+.3f}\n95% CI {r.CI95_low:+.3f} to {r.CI95_high:+.3f}".replace("-", "−")
        ax.text(0.5, 0.985, label,
                transform=ax.transAxes, ha="center", va="top", fontsize=6.2, color=INK)
        clean(ax)

    # B: STIM2.1 SOAR exon, fixed-effect meta-analysis
    ax = fig.add_subplot(bot[1])
    names = {"GSE296712_SHSY5Y": "SH-SY5Y", "GSE230647_iPSC_koloni": "iPSC colonies",
             "GSE77702_iPSC_MN": "iPSC motor neurons", "GSE27394_mouse_SE": "Mouse striatum",
             "Mouse_PE_C2C12": "C2C12", "Mouse_PE_NSC34": "NSC34"}
    d = s15.delta_PSI.to_numpy(); v = s15.variance.to_numpy(); w = 1 / v

    def pool(mask):
        est = np.sum(w[mask] * d[mask]) / np.sum(w[mask]); se = np.sqrt(1 / np.sum(w[mask]))
        return est, est - 1.96 * se, est + 1.96 * se

    all_ = pool(np.ones(len(d), bool)); human = pool((s15.species == "human").to_numpy())
    assert abs(all_[0] - 0.0013) < 5e-4 and abs(human[0] - 0.031) < 5e-4
    for i, (name, est, var, ww) in enumerate(zip(s15.dataset, d, v, w)):
        half = 1.96 * np.sqrt(var)
        ax.plot([est - half, est + half], [i, i], color=MUTED, lw=0.9, zorder=2)
        ax.scatter(est, i, marker="s", s=12 + 60 * ww / w.max(), color=INK, zorder=3)
    for j, (lab, (est, lo, hi)) in enumerate([("Pooled, six datasets", all_), ("Pooled, human only", human)]):
        y = len(d) + 0.6 + j
        ax.add_patch(Polygon([[lo, y], [est, y - 0.28], [hi, y], [est, y + 0.28]], closed=True,
                             color=INK if j == 0 else MUTED, zorder=3))
    ax.axvline(0, color=MUTED, lw=0.6)
    ax.set_xlim(-0.30, 0.40); ax.set_ylim(len(d) + 2.0, -0.7)
    ax.set_yticks(list(range(len(d))) + [len(d) + 0.6, len(d) + 1.6],
                  [names[n] for n in s15.dataset] + ["Pooled, six datasets", "Pooled, human only"])
    ax.get_yticklabels()[len(d)].set_fontweight("bold")
    ax.tick_params(axis="y", length=0, labelsize=6.6)
    ax.set_xlabel("ΔPSI of the 24-nt SOAR exon (KD − control)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="x", direction="out", length=2.5, width=0.6)
    ax.set_title("B  STIM2.1 exon across six datasets", loc="left", weight="bold")

    # C: CBARP junction b per library
    ax = fig.add_subplot(bot[2])
    order = [("SH-SY5Y", "Control", 0), ("SH-SY5Y", "TDP-43 KD", 1),
             ("iPSC colonies", "Control", 2.5), ("iPSC colonies", "TDP-43 KD", 3.5)]
    for dset, grp, x0 in order:
        vals = share[(share.dataset == dset) & (share.group == grp)].sort_values("share_b").share_b.to_numpy()
        group_points(ax, x0, vals, grp)
    ax.set_xticks([0, 1, 2.5, 3.5], ["Control", "KD", "Control", "KD"])
    for xc, lab in [(0.5, "SH-SY5Y"), (3.0, "iPSC colonies")]:
        ax.text(xc, -0.2, lab, ha="center", va="top", fontsize=6.8, transform=ax.get_xaxis_transform())
    ax.set_ylim(-0.03, 1.08); ax.set_xlim(-0.55, 4.05)
    ax.set_ylabel("Share of CBARP exon-4 reads\nusing the alternative 3′ site")
    ax.set_title("C  CBARP exon 4 junction", loc="left", weight="bold")
    clean(ax)

    fig.subplots_adjust(left=0.075, right=0.985, top=0.95, bottom=0.12)
    save(fig, MAIN, "Figure4_SOCE_splicing")


# ------------------------------------------------------------------ Supplementary Figure S3
X_LEFT, X_RIGHT = 1235600, 1234960
E4, E5, ALT_ACC = (1235501, 1235565), (1235001, 1235145), 1235342
KEY = {(1235146, 1235500): ("a", "above"), (1235343, 1235500): ("b", "above"), (1235146, 1235311): ("c", "below")}


def arc(ax, x1, x2, y1, y2, height, lw):
    top = max(y1, y2) + height if height > 0 else min(y1, y2) + height
    verts = [(x1, y1), (x1, top), (x2, top), (x2, y2)]
    ax.add_patch(PathPatch(MPath(verts, [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]),
                           fill=False, lw=lw, edgecolor=INK, capstyle="round", zorder=4))
    return (x1 + x2) / 2, y1 + 0.75 * (top - y1)


def sashimi_track(ax, dataset, group):
    ids = samples[(samples.dataset == dataset) & (samples.group == group)]["sample"]
    m = np.array([covj["coverage"][s]["all"] for s in ids], dtype=float).mean(axis=0)
    x = np.arange(C0, C0 + len(m))
    ymax = m[(x <= X_LEFT) & (x >= X_RIGHT)].max()
    ax.fill_between(x, 0, m, step="mid", color=GROUP[group], alpha=0.85, lw=0, zorder=2)
    counts = junc[(junc.dataset == dataset) & (junc.group == group)].groupby(["start", "end"])["count"].sum()
    donor_tot = counts[counts.index.get_level_values("end") == 1235500].sum()
    acc_tot = counts[counts.index.get_level_values("start") == 1235146].sum()
    ax.set_ylim(-0.62 * ymax, 1.62 * ymax)
    for (s, e), (lab, where) in KEY.items():
        c = int(counts.get((s, e), 0))
        lw = 0.5 + 3.2 * c / (acc_tot if lab == "c" else donor_tot)
        xa, xb = s - 1, e + 1
        if where == "above":
            tx, ty = arc(ax, xb, xa, np.interp(xb, x, m), np.interp(xa, x, m),
                         (0.55 if lab == "a" else 0.30) * ymax, lw)
        else:
            tx, ty = arc(ax, xb, xa, 0, 0, -0.40 * ymax, lw)
        ax.text(tx, ty, f"{lab} {c:,}", ha="center", va="center", fontsize=6.4, color=INK, zorder=6,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.9))
    ax.axhline(0, color=MUTED, lw=0.4, zorder=1)
    ax.set_yticks([0, round(ymax)]); ax.tick_params(axis="y", length=2, pad=1)
    ax.spines[["top", "right", "bottom"]].set_visible(False); ax.spines["left"].set_bounds(0, ymax)
    ax.set_xlim(X_LEFT, X_RIGHT); ax.set_xticks([])
    ax.text(0.005, 0.98, f"{dataset} · {group} (n = {len(ids)})", transform=ax.transAxes,
            ha="left", va="top", fontsize=7, weight="bold", color=INK)
    ax.axvspan(1235146 - 0.5, ALT_ACC + 0.5, color=GRID, alpha=0.35, lw=0, zorder=0)


def gene_models(ax):
    ax.set_xlim(X_LEFT, X_RIGHT); ax.set_ylim(-0.5, 4.3); ax.axis("off")
    rows = [("MANE Select (ENST00000650044)", [E4, E5], False),
            ("Retained intron (ENST00000589260)", [E4, (1235001, 1235342)], False),
            ("NMD-annotated (ENST00000648750)", [E4, (1235337, 1235390), E5], False),
            ("Structure implied by junction c", [E4, (1235312, 1235342), E5], True)]
    for i, (name, exons, dashed) in enumerate(rows):
        y = 3.55 - i * 1.1
        ax.plot([max(e[1] for e in exons), min(e[0] for e in exons)], [y, y], color=MUTED, lw=0.6,
                ls=(0, (2, 1.5)) if dashed else "-", zorder=1)
        for a, b in exons:
            ax.add_patch(Rectangle((a - 0.5, y - 0.22), b - a + 1, 0.44, fc="white" if dashed else "#4D4D4D",
                                   ec="#4D4D4D", lw=0.7, ls="--" if dashed else "-", zorder=2))
        ax.text(1235598, y + 0.30, name, ha="left", va="bottom", fontsize=6.2, color=MUTED)
    ax.axvspan(1235146 - 0.5, ALT_ACC + 0.5, color=GRID, alpha=0.35, lw=0, zorder=0)


def supplementary_s3():
    fig = plt.figure(figsize=(7.1, 7.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[5.0, 2.0], hspace=0.28)
    top = gs[0].subgridspec(6, 1, height_ratios=[0.42, 1, 1, 1, 1, 1.45], hspace=0.12)
    ax_h = fig.add_subplot(top[0]); ax_h.set_xlim(X_LEFT, X_RIGHT); ax_h.set_ylim(0, 1); ax_h.axis("off")
    ax_h.text(0.0, 1.0, "A  CBARP exon 4–exon 5 region: coverage and junction reads", transform=ax_h.transAxes,
              ha="left", va="top", weight="bold", fontsize=8.5)
    for (a, b), lab in [(E4, "Exon 4"), (E5, "Exon 5")]:
        ax_h.text((a + b) / 2, 0.02, lab, ha="center", va="bottom", fontsize=6.8, color=INK)
    ax_h.text((1235146 + ALT_ACC) / 2, 0.02, "alternative 3′ region", ha="center", va="bottom",
              fontsize=6.6, color=MUTED)
    axes = []
    for i, (d, g) in enumerate([("SH-SY5Y", "Control"), ("SH-SY5Y", "TDP-43 KD"),
                                ("iPSC colonies", "Control"), ("iPSC colonies", "TDP-43 KD")]):
        ax = fig.add_subplot(top[i + 1]); sashimi_track(ax, d, g); axes.append(ax)
    axes[0].set_ylabel("Mean depth", labelpad=2); axes[2].set_ylabel("Mean depth", labelpad=2)
    ax_m = fig.add_subplot(top[5]); gene_models(ax_m)
    ax_x = ax_m.secondary_xaxis(-0.02)
    ticks = [1235600, 1235500, 1235400, 1235300, 1235200, 1235100, 1235000]
    ax_x.set_xticks(ticks, [f"{t:,}" for t in ticks]); ax_x.tick_params(length=2, labelsize=6.3)
    ax_x.set_xlabel("chr19 position (GRCh38), drawn 5′→3′ for the minus-strand gene", fontsize=6.8, labelpad=2)

    ev = pd.read_csv(ROOT / "supplementary" / "S3_rMATS_significant_events.csv.gz")
    q = ev[(ev.model == "JC") & (ev.gene.str.upper() == "CBARP") & (ev.passes_coverage_filter == True)
           & (ev.FDR < 0.05) & (ev.delta_PSI.abs() >= 0.10)]
    assert len(q) == 32 and q.dataset.nunique() == 5
    ax = fig.add_subplot(gs[1])
    order = ["SH-SY5Y (GSE296712)", "iPSC colonies (GSE230647)", "Mouse striatum (GSE27394)",
             "C2C12 (GSE171714)", "NSC34 (GSE171714)"]
    labels = ["SH-SY5Y", "iPSC colonies", "Mouse striatum", "C2C12", "NSC34"]
    for yi, name in enumerate(order):
        vals = q.loc[q.dataset == name, "delta_PSI"].sort_values().to_numpy()
        offs = np.linspace(-0.2, 0.2, len(vals)) if len(vals) > 1 else np.array([0.])
        ax.scatter(vals, yi + offs, s=18, zorder=3, lw=0.4, edgecolor=INK, color="#8C8C8C")
        ax.text(0.93, yi, f"n = {len(vals)}", va="center", ha="right", fontsize=6.4, color=MUTED)
    ax.axvline(0, color=MUTED, lw=0.6)
    ax.set_xlim(-0.85, 0.95); ax.set_ylim(4.55, -0.55)
    ax.set_yticks(range(5), labels)
    ax.set_xlabel("rMATS event ΔPSI (KD − control)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", direction="out", length=2.5, width=0.6)
    ax.set_title("B  Coverage-qualified rMATS CBARP events", loc="left", weight="bold")
    fig.subplots_adjust(left=0.13, right=0.985, top=0.99, bottom=0.06)
    save(fig, SUPP, "Supplementary_Figure_S3_CBARP_locus", pdf=True)


# ------------------------------------------------------------------ Figure 6
def box(ax, x, y, w, h, title, body, edge):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.015",
                                edgecolor=edge, facecolor="white", linewidth=0.9))
    ax.text(x + w / 2, y + h * 0.68, title, ha="center", va="center", weight="bold", fontsize=8)
    ax.text(x + w / 2, y + h * 0.30, body, ha="center", va="center", fontsize=7.2)


def arrow(ax, p1, p2, dashed=False, color=INK):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=9, lw=0.9, color=color,
                                 linestyle=(0, (3, 2)) if dashed else "-", shrinkA=0, shrinkB=0))


def figure6():
    fig, ax = plt.subplots(figsize=(7.1, 3.9))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.01, 0.985, "Measured in the laboratory cultures", va="top", weight="bold", fontsize=8.5)
    box(ax, 0.36, 0.76, 0.28, 0.15, "shRNA knockdown", "TARDBP mRNA ↓ 94%", BLUE)
    box(ax, 0.02, 0.50, 0.30, 0.17, "RT-qPCR", "TRPC1 / STIM1 / ORAI1 / ATP2A3 ↑", BLUE)
    box(ax, 0.68, 0.50, 0.30, 0.17, "Fura-2", "Ca$^{2+}$ readdition ↓; ER release ↓ (n.s.)", BLUE)
    arrow(ax, (0.42, 0.76), (0.22, 0.67)); arrow(ax, (0.58, 0.76), (0.78, 0.67))
    ax.text(0.01, 0.42, "Candidate mechanisms from public RNA-seq (hypotheses)", va="top", weight="bold",
            fontsize=8.5)
    box(ax, 0.02, 0.06, 0.30, 0.19, "Parallel Ca$^{2+}$ entry route", "CBARP splicing → CaV (not measured)", ORANGE)
    box(ax, 0.35, 0.06, 0.30, 0.19, "ER refilling", "ATP2A2 ↓ / bioenergetics", ORANGE)
    box(ax, 0.68, 0.06, 0.30, 0.19, "Channel composition and feedback", "ORAI3 share ↑ / SARAF ↑", ORANGE)
    arrow(ax, (0.55, 0.25), (0.76, 0.50), dashed=True, color=MUTED)
    arrow(ax, (0.83, 0.25), (0.83, 0.50), dashed=True, color=MUTED)
    ax.plot([0.66, 0.71], [0.965, 0.965], color=INK, lw=0.9)
    ax.text(0.715, 0.965, "measured link", va="center", fontsize=6.6)
    ax.plot([0.83, 0.88], [0.965, 0.965], color=MUTED, lw=0.9, ls=(0, (3, 2)))
    ax.text(0.885, 0.965, "hypothesis to test", va="center", fontsize=6.6, color=MUTED)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.02)
    save(fig, MAIN, "Figure6_working_model")


if __name__ == "__main__":
    figure4()
    supplementary_s3()
    figure6()
    print("written: Figure4_SOCE_splicing, Supplementary_Figure_S3_CBARP_locus, Figure6_working_model")

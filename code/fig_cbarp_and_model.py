"""Source-backed CBARP main figure and a labelled working-model schematic.

Figure 5 reads Supplementary Table S3. No junction curve or trace is
reconstructed from aggregate counts. Figure 6 is conceptual and uses no
synthetic measurements.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "main"
OUT.mkdir(exist_ok=True)
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREY = "#585858"
LIGHT = "#E6E6E6"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "axes.linewidth": 0.7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "svg.fonttype": "none",
    "savefig.facecolor": "white",
})


def save(fig, stem):
    for ext in ("png", "svg"):
        path = OUT / f"{stem}.{ext}"
        fig.savefig(path, dpi=400, bbox_inches="tight", pad_inches=0.09)
        if ext == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)


def cb_splicing():
    all_events = pd.read_csv(ROOT / "supplementary" / "S3_rMATS_significant_events.csv.gz")
    q = all_events[
        (all_events.model == "JC")
        & (all_events.gene.str.upper() == "CBARP")
        & (all_events.passes_coverage_filter == True)
        & (all_events.FDR < 0.05)
        & (all_events.delta_PSI.abs() >= 0.10)
    ].copy()
    assert len(q) == 32 and q.dataset.nunique() == 5
    assert abs(q.delta_PSI.abs().median() - 0.37) < 0.001
    ev = q[(q.dataset == "iPSC colonies (GSE230647)") & (q.event_class == "SE") & (q.ID == 98689)]
    assert len(ev) == 1
    ev = ev.iloc[0]
    exon_length = int(ev.exonEnd - ev.exonStart_0base)
    genomic_label = f"{ev.chr}:{int(ev.exonStart_0base) + 1:,}–{int(ev.exonEnd):,}  ({ev.strand} strand)"
    ctrl = np.array([float(x) for x in ev.PSI_control_per_replicate.split(",")])
    kd = np.array([float(x) for x in ev.PSI_knockdown_per_replicate.split(",")])
    assert len(ctrl) == len(kd) == 4

    fig = plt.figure(figsize=(7.1, 4.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 2.7], width_ratios=[1, 1.38],
                          hspace=0.17, wspace=0.24)
    ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[1, 0]); ax_c = fig.add_subplot(gs[:, 1])

    ax_a.set_xlim(0, 10); ax_a.set_ylim(0, 2); ax_a.axis("off")
    ax_a.text(0, 2.0, "A  CBARP skipped-exon structure", va="top", weight="bold", fontsize=9)
    for x, w, label in [(0.6, 1.8, "Exon"), (4.15, 1.65, f"{exon_length} bp"), (7.8, 1.7, "Exon")]:
        ax_a.add_patch(Rectangle((x, 0.72), w, 0.43, facecolor="white", edgecolor="black", lw=0.8))
        ax_a.text(x + w / 2, 0.93, label, ha="center", va="center", fontsize=6.6)
    ax_a.plot([2.4, 4.15], [0.93, 0.93], color="black", lw=0.8)
    ax_a.plot([5.8, 7.8], [0.93, 0.93], color="black", lw=0.8)
    ax_a.text(5, 0.25, genomic_label, ha="center", fontsize=7)

    x_ctrl = np.array([0.86, 0.96, 1.06, 1.16]); x_kd = np.array([1.86, 1.96, 2.06, 2.16])
    ax_b.scatter(x_ctrl, ctrl, facecolor="white", edgecolor=BLUE, s=35, lw=1.15, zorder=3)
    ax_b.scatter(x_kd, kd, facecolor=ORANGE, edgecolor="black", s=35, lw=0.45, zorder=3)
    ax_b.plot([0.82, 1.20], [ctrl.mean()] * 2, color=BLUE, lw=1.4)
    ax_b.plot([1.82, 2.20], [kd.mean()] * 2, color=ORANGE, lw=1.4)
    ax_b.set_xlim(0.55, 2.5); ax_b.set_ylim(-0.04, 0.79)
    ax_b.set_xticks([1.01, 2.01], ["Control", "TDP-43 KD"])
    ax_b.set_ylabel("Exon inclusion (PSI)")
    ax_b.set_title("B  iPSC colony libraries", loc="left", weight="bold")
    ax_b.text(1.46, 0.73, f"ΔPSI {ev.delta_PSI:+.3f}", ha="center", fontsize=8)
    ax_b.spines[["top", "right"]].set_visible(False)
    ax_b.tick_params(direction="out", length=3, width=0.7)

    order = ["SH-SY5Y (GSE296712)", "iPSC colonies (GSE230647)",
             "Mouse striatum (GSE27394)", "C2C12 (GSE171714)", "NSC34 (GSE171714)"]
    labels = ["SH-SY5Y", "iPSC colonies", "Mouse striatum", "C2C12", "NSC34"]
    for yi, name in enumerate(order):
        vals = q.loc[q.dataset == name, "delta_PSI"].sort_values().to_numpy()
        offsets = np.linspace(-0.16, 0.16, len(vals)) if len(vals) > 1 else np.array([0.])
        for value, off in zip(vals, offsets):
            ax_c.scatter(value, yi + off, color=ORANGE if value > 0 else BLUE,
                         edgecolor="black", linewidth=0.3, s=29, zorder=3)
        ax_c.text(0.81, yi, f"{len(vals)} events", va="center", fontsize=7)
    ax_c.axvline(0, color=GREY, lw=0.8)
    ax_c.set_xlim(-0.86, 0.95); ax_c.set_ylim(-0.5, 4.5)
    ax_c.set_yticks(range(5), labels)
    ax_c.invert_yaxis()
    ax_c.set_xlabel("Event-level ΔPSI (TDP-43 KD − control)")
    ax_c.set_title("C  Coverage-qualified CBARP events", loc="left", weight="bold")
    ax_c.spines[["top", "right", "left"]].set_visible(False)
    ax_c.tick_params(axis="y", length=0)
    ax_c.tick_params(axis="x", direction="out", length=3, width=0.7)
    save(fig, "Figure5_CBARP_splicing")


def box(ax, x, y, w, h, title, body, edge=GREY, fill="white"):
    rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.018",
                          edgecolor=edge, facecolor=fill, linewidth=0.9)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h*0.70, title, ha="center", va="center", weight="bold", fontsize=8)
    ax.text(x + w/2, y + h*0.30, body, ha="center", va="center", fontsize=7.4)


def arrow(ax, x1, y1, x2, y2, dashed=False, color=GREY):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=10,
                                 linewidth=1, linestyle="--" if dashed else "-", color=color))


def model():
    fig, ax = plt.subplots(figsize=(7.1, 3.75))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.02, 0.98, "Laboratory observations", va="top", weight="bold", fontsize=9)
    box(ax, 0.37, 0.75, 0.26, 0.18, "Laboratory shRNA", "TARDBP mRNA ↓", edge=BLUE)
    box(ax, 0.02, 0.45, 0.27, 0.20, "Target-gene qPCR", "TRPC1 / STIM1 / ORAI1 ↑", edge=BLUE)
    box(ax, 0.71, 0.45, 0.27, 0.20, "Fura-2 readdition", "Ca2+ amplitude ↓", edge=BLUE)
    arrow(ax, 0.42, 0.75, 0.22, 0.65, color=BLUE)
    arrow(ax, 0.58, 0.75, 0.78, 0.65, color=BLUE)
    ax.text(0.02, 0.38, "Public RNA-seq context and testable mechanisms", va="top",
            weight="bold", fontsize=9)
    box(ax, 0.02, 0.10, 0.27, 0.20, "ER refilling", "ATP2A2 ↓ / bioenergetics", edge=ORANGE)
    box(ax, 0.37, 0.10, 0.26, 0.20, "Channel regulation", "ORAI3 share ↑ / SARAF ↑", edge=ORANGE)
    box(ax, 0.71, 0.10, 0.27, 0.20, "Parallel pathway", "CBARP splicing / CaV", edge=ORANGE)
    fig.subplots_adjust(left=0.015, right=0.985, top=0.985, bottom=0.06)
    save(fig, "Figure6_working_model")


if __name__ == "__main__":
    cb_splicing()
    model()

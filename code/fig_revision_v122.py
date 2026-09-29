"""Figures redrawn for release v1.0.8 from files in this repository.

Figure 2                 Fura-2: representative recordings on common axes, and the three wells per group
                         for ER release, Ca2+ readdition and their ratio (one culture plate)
Figure 5                 ALS tissue heat map with the number of cases and controls in each row label
Supplementary Figure S1  example detection-power simulation, with its assumptions printed on the figure
Supplementary Figure S7  TRPC1 across neurological cohorts, grouped by disease; MS re-analyses separated
Supplementary Figure S9  the two original Fura-2 recordings at their own axis ranges

Inputs: supplementary/S1_laboratory_source_data.xlsx, source_data/fura2_traces/*,
tables/Table5_cross_disease_comparison.csv, supplementary/S16b_multiple_sclerosis_donor_level.csv,
source_data/power_simulation_S1.csv.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "figures" / "main"
SUPP = ROOT / "figures" / "supplementary"
BLUE, ORANGE, INK, MUTED = "#0072B2", "#D55E00", "#222222", "#6B6B6B"
PALE_BLUE, PALE_ORANGE = "#96C4DE", "#F1BE96"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.titlesize": 8.5, "axes.titleweight": "bold",
    "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.7,
    "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
    "pdf.fonttype": 42, "savefig.facecolor": "white", "figure.facecolor": "white",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
})


def save(fig, folder, stem, exts=("png", "svg"), dpi=400):
    for ext in exts:
        path = folder / f"{stem}.{ext}"
        fig.savefig(path, dpi=dpi)
        if ext == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)


# ------------------------------------------------------------------ Figure 2
def figure2():
    src = ROOT / "source_data" / "fura2_traces"
    per_well = pd.read_excel(ROOT / "supplementary/S1_laboratory_source_data.xlsx", sheet_name="Fura2_per_culture")

    def trace(name):
        d = np.loadtxt(src / f"PRISM_RAW_{name}_130626.txt")
        t, r = d[:, 0], d[:, 1]
        assert np.allclose(np.median(np.diff(t)), 1.0)
        k = 15
        smooth = np.convolve(r, np.ones(k) / k, mode="same")
        onset = t[k + np.argmax(np.gradient(smooth, t)[k:-k])]
        return t - onset, r, onset

    ct, cr, con = trace("CONTROL")
    kt, kr, kon = trace("KD")
    assert 2600 < con < 2700 and 1850 < kon < 1950

    fig = plt.figure(figsize=(6.05, 5.55))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], hspace=0.55, wspace=0.55,
                          left=0.095, right=0.985, top=0.955, bottom=0.115)
    ax = fig.add_subplot(gs[0, :])
    ax.plot(ct, cr, color=BLUE, lw=0.55, label="Non-targeting shRNA control (one well)")
    ax.plot(kt, kr, color=ORANGE, lw=0.55, label="shTDP-43 (one well)")
    ax.axvline(0, color=MUTED, lw=0.6, ls=(0, (3, 2)))
    ax.text(18, 2.62, "Ca$^{2+}$ readdition", fontsize=6.8, color=MUTED, va="top")
    ax.set_xlim(-1300, 950)
    ax.set_ylim(0.6, 2.7)
    ax.set_xlabel("Time relative to the steepest point of the readdition rise (s)")
    ax.set_ylabel("F$_{340}$/F$_{380}$")
    ax.tick_params(direction="out", length=2.5, width=0.6)
    ax.legend(frameon=False, loc="upper left", fontsize=6.8, handlelength=1.6)
    ax.set_title("A  Representative recordings on common axes", loc="left")

    labels = ["Non-targeting\nshRNA", "shTDP-43"]
    groups = ["Non-targeting shRNA control", "shTDP-43"]
    specs = [("ER_release", "B  ER Ca$^{2+}$ release", "Δ(F$_{340}$/F$_{380}$)", (0, 0.5)),
             ("readdition", "C  Ca$^{2+}$ readdition", "Δ(F$_{340}$/F$_{380}$)", (0, 2.3)),
             ("readdition_to_release_ratio", "D  Readdition / release", "ratio per well", (0, 9.5))]
    for j, (col, title, ylab, ylim) in enumerate(specs):
        ax = fig.add_subplot(gs[1, j])
        for i, (g, colr) in enumerate(zip(groups, (BLUE, ORANGE))):
            v = per_well.loc[per_well.group == g, col].to_numpy(float)
            assert len(v) == 3
            ax.bar(i, v.mean(), 0.6, color=colr, alpha=0.9)
            ax.errorbar(i, v.mean(), yerr=v.std(ddof=1) / np.sqrt(3), color=INK, capsize=3, lw=0.9)
            ax.scatter(i + np.linspace(-0.13, 0.13, 3), v, s=17, facecolor="white", edgecolor=INK, lw=0.7, zorder=3)
        ax.set_xticks([0, 1], labels)
        ax.set_ylim(*ylim)
        ax.set_ylabel(ylab)
        ax.set_title(title, loc="left")
        ax.tick_params(direction="out", length=2.5, width=0.6)
    fig.text(0.5, 0.012, "B–D: mean ± well-to-well SEM; dots, the three wells per group of one culture plate; descriptive",
             ha="center", fontsize=6.8, color=MUTED)
    save(fig, MAIN, "Figure2_calcium_responses")


# ------------------------------------------------------------------ Supplementary Figure S9
def supp_s9():
    im = Image.open(ROOT / "source_data" / "fura2_traces" / "Figure2_panels_A-D.png").convert("RGB")
    top = im.crop((0, 40, im.width, 1150))
    top.save(SUPP / "Supplementary_Figure_S9_original_Fura2_recordings.png", dpi=(300, 300))
    top.save(SUPP / "Supplementary_Figure_S9_original_Fura2_recordings.pdf", resolution=300)
    return top.size


# ------------------------------------------------------------------ Supplementary Figure S1
def supp_s1():
    tsv = pd.read_csv("/Users/elmas/Desktop/MAKALE/03_TABLOLAR/guc_simulasyonu.tsv", sep="\t")
    dps = [0.05, 0.10, 0.15, 0.20, 0.30]
    long = pd.DataFrame([dict(reads_per_sample=int(r.derinlik), true_delta_PSI=d, power=r[f"dPSI_{d:.2f}"])
                         for _, r in tsv.iterrows() for d in dps])
    long["n_simulations"] = 2000
    long.to_csv(ROOT / "source_data" / "power_simulation_S1.csv", index=False)
    P = long.pivot(index="reads_per_sample", columns="true_delta_PSI", values="power")

    fig = plt.figure(figsize=(5.9, 4.95))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.6, 1.15], hspace=0.27, left=0.115, right=0.975, top=0.94, bottom=0.02)
    ax = fig.add_subplot(gs[0])
    for depth, colr in zip([10, 20, 50, 100], ["#000000", "#0072B2", "#009E73", "#D55E00"]):
        ax.plot(dps, [P.loc[depth, d] for d in dps], "o-", color=colr, lw=2, ms=5, label=f"{depth} reads/sample")
    ax.axhline(0.8, color="#333", ls="--", lw=1)
    ax.annotate("80% power", xy=(0.055, 0.82), fontsize=8)
    ax.axvline(0.10, color="black", ls=":", lw=1.2)
    ax.annotate("|ΔPSI| = 0.10", xy=(0.105, 0.41), fontsize=8)
    ax.set_xlabel("true ΔPSI (knockdown − control)")
    ax.set_ylabel("simulated detection power")
    ax.set_ylim(0, 1.02)
    ax.set_title("Example detection-power simulation, 3 + 3 design", loc="left")
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    tx = fig.add_subplot(gs[1])
    tx.axis("off")
    text = ("Simulation assumptions (an example, not the rMATS analysis)\n"
            "• three versus three replicates; replicate PSI drawn from a normal distribution around 0.5 ± ΔPSI/2,\n"
            "   between-replicate SD 0.05, limited to 0.01–0.99\n"
            "• binomial read sampling at a fixed depth of 10, 20, 50 or 100 informative reads per sample\n"
            "• two-sample t test (equal variances) on PSI; a hit is p < 0.05 without multiple-testing correction\n"
            "• 2,000 simulations per point; no FDR correction and no |ΔPSI| ≥ 0.10 threshold")
    tx.text(0.0, 1.0, text, va="top", ha="left", fontsize=7, linespacing=1.45, color=INK)
    save(fig, SUPP, "Supplementary_Figure_S1_detection_power", exts=("png", "svg", "pdf"))


# ------------------------------------------------------------------ Figure 5
def figure5():
    d = pd.read_csv(ROOT / "tables" / "Table5_cross_disease_comparison.csv")
    d = d[(d.group == "ALS") & d.gene.isin(["TRPC1", "SARAF", "CBARP"])].copy()
    regions = ["Hippocampus", "Cerebellum", "Motor cortex (lateral)", "Frontal cortex", "Motor cortex (medial)",
               "Temporal cortex", "Occipital cortex", "Spinal cord (lumbar)", "Spinal cord (cervical)",
               "Spinal cord (thoracic)"]
    genes = ["TRPC1", "SARAF", "CBARP"]
    eff = d.pivot(index="region", columns="gene", values="cliffs_delta").loc[regions, genes]
    q = d.pivot(index="region", columns="gene", values="q_value").loc[regions, genes]
    n = d.drop_duplicates("region").set_index("region").loc[regions, ["n_case", "n_control"]].astype(int)
    assert eff.notna().all().all() and q.notna().all().all()
    fig, ax = plt.subplots(figsize=(7.1, 5.5), layout="constrained")
    cmap = LinearSegmentedColormap.from_list("accessible_diverging", ["#95C4DF", "#FFFFFF", "#F1B481"])
    im = ax.imshow(eff.to_numpy(), vmin=-1, vmax=1, cmap=cmap, aspect="auto")
    for i, region in enumerate(regions):
        for j, g in enumerate(genes):
            mark = "*" if q.loc[region, g] < 0.05 else ""
            ax.text(j, i, f"{eff.loc[region, g]:+.2f}{mark}".replace("-", "−"), ha="center", va="center",
                    fontsize=7.8, weight="bold" if mark else "normal")
    ax.set_xticks(range(3), genes, fontsize=8, weight="bold")
    ax.set_yticks(range(len(regions)), [f"{r}  (n = {n.loc[r, 'n_case']}/{n.loc[r, 'n_control']})" for r in regions])
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, length=0)
    ax.axhline(6.5, color="#555555", lw=0.65)
    ax.set_title("ALS tissue versus non-neurological controls, unadjusted regional effects", loc="left", pad=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.042, pad=0.035, ticks=[-1, -0.5, 0, 0.5, 1])
    cb.set_label("Cliff's δ (case − control)")
    fig.supxlabel("n = ALS cases/controls in that region. * Benjamini–Hochberg q < 0.05 within the region "
                  "(correction across the genes tested there; Table 5).", fontsize=6.8, color=MUTED, x=0.02, ha="left")
    save(fig, MAIN, "Figure5_ALS_expression", dpi=300)


# ------------------------------------------------------------------ Supplementary Figure S7
def supp_s7():
    t5 = pd.read_csv(ROOT / "tables" / "Table5_cross_disease_comparison.csv")
    t5 = t5[t5.gene == "TRPC1"]
    s16b = pd.read_csv(ROOT / "supplementary" / "S16b_multiple_sclerosis_donor_level.csv")

    def row(label, delta, n, sig=False, dagger=False):
        return dict(label=label, delta=float(delta), n=n, sig=bool(sig), dagger=bool(dagger))

    als_regions = ["Hippocampus", "Cerebellum", "Motor cortex (lateral)", "Frontal cortex", "Motor cortex (medial)",
                   "Temporal cortex", "Occipital cortex", "Spinal cord (lumbar)", "Spinal cord (cervical)",
                   "Spinal cord (thoracic)"]
    blocks = []
    a_ = t5[t5.group == "ALS"].set_index("region")
    blocks.append(("ALS (NYGC cohort)", [row(r, a_.loc[r, "cliffs_delta"], f"{int(a_.loc[r,'n_case'])}/{int(a_.loc[r,'n_control'])}",
                                             a_.loc[r, "q_value"] < 0.05) for r in als_regions]))
    o_ = t5[t5.group == "ONd"].set_index("region")
    blocks.append(("Other neurological disorders (NYGC cohort; diagnoses not public)",
                   [row(r, o_.loc[r, "cliffs_delta"], f"{int(o_.loc[r,'n_case'])}/{int(o_.loc[r,'n_control'])}",
                        o_.loc[r, "q_value"] < 0.05) for r in ["Cerebellum", "Frontal cortex", "Temporal cortex"]]))
    ad = t5[t5.group == "AD"].iloc[0]
    pdz = t5[t5.group == "PD"].iloc[0]
    blocks.append(("Alzheimer's disease (GSE125583)",
                   [row("Fusiform gyrus", ad.cliffs_delta, f"{int(ad.n_case)}/{int(ad.n_control)}", ad.q_value < 0.05)]))
    blocks.append(("Parkinson's disease (GSE68719)",
                   [row("Prefrontal cortex (BA9)", pdz.cliffs_delta, f"{int(pdz.n_case)}/{int(pdz.n_control)}", pdz.q_value < 0.05)]))
    m = t5[t5.group == "MS"].set_index("region")
    don = m.loc["Donor level (10 MS vs 5 control donors)"]
    pooled = m.loc["All five regions pooled (centred within region)"]
    blocks.append(("Multiple sclerosis (two cohorts)",
                   [row("GSE138614 white matter, donor level", don.cliffs_delta, "10/5 donors", don.q_value < 0.05),
                    row("GSE123496, five regions pooled", pooled.cliffs_delta, "25/25 samples", pooled.q_value < 0.05)]))
    lesion = m.loc["MS lesions vs control white matter"]
    nawm = m.loc["Normal-appearing white matter vs control white matter"]
    g = lambda comp, unit: s16b[(s16b.comparison == comp) & (s16b.gene == "TRPC1") & (s16b.unit == unit)].iloc[0]
    nawm_d = g("NAWM vs control WM", "donor")
    adj_d = g("MS vs control, myelin+glia-adjusted TRPC1", "donor")
    sens = [("GSE138614: re-analyses of the samples and donors of the row above (not independent replicates)",
             [row("Lesions vs control WM, sample level", lesion.cliffs_delta, "52/25", lesion.q_value < 0.05),
              row("NAWM vs control WM, sample level", nawm.cliffs_delta, "21/25", nawm.q_value < 0.05),
              row("NAWM vs control WM, donor level", nawm_d.cliffs_delta, "7/5 donors", False, nawm_d.p < 0.05),
              row("All MS, MBP/PLP1/GFAP-adjusted, donor level", adj_d.cliffs_delta, "10/5 donors", False, adj_d.p < 0.05)])]
    assert abs(don.cliffs_delta + 0.84) < 1e-6 and abs(nawm_d.cliffs_delta + 0.771) < 1e-6

    left, right = 0.375, 0.875
    hA = sum(len(r) + 1.15 for _, r in blocks)
    hB = sum(len(r) + 1.15 for _, r in sens)
    fig = plt.figure(figsize=(7.1, 9.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[hA, hB + 0.8], hspace=0.24, left=left, right=right, top=0.945, bottom=0.115)
    x_hdr = (0.012 - left) / (right - left)

    def draw(ax, blks):
        y = 0
        ticks, labs = [], []
        for title, rows in blks:
            ax.text(x_hdr, y, title, ha="left", va="center", fontsize=7.3, weight="bold", transform=ax.get_yaxis_transform())
            y += 0.85
            for r in rows:
                v = r["delta"]
                pos = v > 0
                colr = (ORANGE if pos else BLUE) if r["sig"] else (PALE_ORANGE if pos else PALE_BLUE)
                ax.barh(y, v, 0.68, color=colr, edgecolor=(ORANGE if pos else BLUE) if r["dagger"] else "none", lw=0.9)
                mark = "*" if r["sig"] else ("†" if r["dagger"] else "")
                ax.text(v + (0.02 if pos else -0.02), y, f"{v:+.3f}{mark}".replace("-", "−"), ha="left" if pos else "right",
                        va="center", fontsize=6.8, weight="bold" if r["sig"] else "normal")
                ax.text(1.03, y, r["n"], transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6.4, color=MUTED)
                ticks.append(y)
                labs.append(r["label"])
                y += 1
            y += 0.3
        ax.set_yticks(ticks, labs, fontsize=6.9)
        ax.set_ylim(y - 0.2, -0.6)
        ax.axvline(0, color="#333", lw=0.9)
        ax.set_xlim(-1.2, 0.95)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.text(1.03, -0.6, "n, case/control", transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=6.4, color=MUTED)

    ax1 = fig.add_subplot(gs[0])
    draw(ax1, blocks)
    ax2 = fig.add_subplot(gs[1])
    draw(ax2, sens)
    ax2.set_xlabel("TRPC1 Cliff's δ (case − control)")
    ax1.set_xlabel("TRPC1 Cliff's δ (case − control)")
    ax1.text(0.0, 1.012, "← decreased", transform=ax1.transAxes, ha="left", va="bottom", fontsize=7, color=MUTED)
    ax1.text(1.0, 1.012, "increased →", transform=ax1.transAxes, ha="right", va="bottom", fontsize=7, color=MUTED)
    fig.text(0.012, ax1.get_position().y1 + 0.028, "A  Case–control comparisons by disease", fontsize=8.5, weight="bold")
    fig.text(0.012, ax2.get_position().y1 + 0.028, "B  Multiple-sclerosis sensitivity analyses", fontsize=8.5, weight="bold")
    fig.text(0.012, 0.012,
             "Dark bars: Benjamini–Hochberg q < 0.05 (* within each cohort's comparison family); pale bars: not significant; "
             "† uncorrected p < 0.05 only.\n"
             "NAWM, normal-appearing white matter; WM, white matter; BA9, Brodmann area 9. Cohorts differ in tissue compartment, sample size and\n"
             "processing; the display is descriptive and the values are not pooled.",
             fontsize=6.4, color=MUTED, va="bottom", ha="left", linespacing=1.4)
    save(fig, SUPP, "Supplementary_Figure_S7_cross_disease_TRPC1", exts=("png", "svg", "pdf"), dpi=300)


if __name__ == "__main__":
    figure2()
    print("S9 size:", supp_s9())
    supp_s1()
    figure5()
    supp_s7()
    print("written: Figure 2, Figure 5, Supplementary S1, S7, S9")

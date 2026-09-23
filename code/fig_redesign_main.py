#!/usr/bin/env python3
"""v3 — main figures for the submission package.

Figure 1 (detection power), Figure 2 (TRPC1 robustness), Figure 3 (robust SOCE
candidates) and Figure 7 (transcript-family abundance, composition-adjusted TPM).
Numbering follows MANUSCRIPT_v4_SUBMISSION; figure numbers are not burned into the images."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from importlib.machinery import SourceFileLoader

core = SourceFileLoader("core", "/Users/elmas/Desktop/MAKALE/04_KOD/01_rmats_core.py").load_module()
TAB = os.path.join(core.OUT, "03_TABLOLAR")
FIG = "/Users/elmas/Desktop/MAKALE/11_NEUROCHEMISTRY_INTERNATIONAL_FIGURE_REVISION/figures/redesign_work"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.titlesize": 9, "axes.titleweight": "bold",
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "savefig.facecolor": "white", "figure.facecolor": "white",
})
KD_C, CT_C, ACC, WARN = "#D55E00", "#0072B2", "#009E73", "#E69F00"

# ================================================================ FIG 2: TRPC1
fig, axes = plt.subplots(3, 1, figsize=(7.1, 8.8), layout="constrained",
                         gridspec_kw={"height_ratios": [1.0, 1.0, 1.15]})

# --- 1A: replika duzeyi PSI
ax = axes[0]
kd = [1.000, 0.778, 1.000]; ct = [1.000, 1.000, 0.455]
ct_x = [-0.08, 0.05, 0.18]          # ust uste binen noktalari ayirmak icin jitter
kd_x = [0.82, 0.95, 1.08]
ax.scatter(ct_x, ct, s=70, color=CT_C, zorder=3, label="Control")
ax.scatter(kd_x, kd, s=70, color=KD_C, zorder=3, label="TDP-43 KD")
ax.hlines(np.mean(ct), -0.15, 0.25, color=CT_C, lw=2)
ax.hlines(np.mean(kd), 0.75, 1.15, color=KD_C, lw=2)
for xs, v, n in [(ct_x, ct, ["2/0", "10/0", "10/6"]),
                 (kd_x, kd, ["25/0", "7/1", "28/0"])]:
    for xi, yi, ni in zip(xs, v, n):
        ax.annotate(ni, (xi, yi), textcoords="offset points", xytext=(0, 11),
                    fontsize=7, color="#555", ha="center")
ax.set_xlim(-0.4, 1.45); ax.set_ylim(0.3, 1.12)
ax.set_xticks([0.05, 0.95]); ax.set_xticklabels(["Control", "TDP-43 KD"])
ax.set_ylabel("PSI (inclusion level)")
ax.set_title("A  Replicate-level PSI", loc="left")
ax.text(0.02, 0.05, "labels: IJC/SJC", transform=ax.transAxes, fontsize=7, color="#777")

# --- 1B: bootstrap dagilimi
ax = axes[1]
df = core.load_event("GSE296712_SHSY5Y", "SE", "JC")
df["gene_up"] = df["geneSymbol"].str.upper()
r = df[(df.gene_up == "TRPC1") & (df.ID == 33987)].iloc[0]
rng = np.random.default_rng(42)
il, sl = r["IncFormLen"], r["SkipFormLen"]
kdi = np.array(r["IJC_SAMPLE_1_list"]); kds = np.array(r["SJC_SAMPLE_1_list"])
cti = np.array(r["IJC_SAMPLE_2_list"]); cts = np.array(r["SJC_SAMPLE_2_list"])
boot = []
for _ in range(20000):
    a = rng.integers(0, 3, 3); b = rng.integers(0, 3, 3)
    p1 = core.psi_from_counts(kdi[a].sum(), kds[a].sum(), il, sl)
    p2 = core.psi_from_counts(cti[b].sum(), cts[b].sum(), il, sl)
    if not (np.isnan(p1) or np.isnan(p2)):
        boot.append(p1 - p2)
boot = np.array(boot)
lo, hi = np.percentile(boot, [2.5, 97.5])
ax.hist(boot, bins=60, color="#b8c9c7", edgecolor="none")
ax.axvline(0, color="#333", lw=1.2, ls="--")
ax.axvline(0.108, color=KD_C, lw=2, label="rMATS ΔPSI = +0.108")
ax.axvspan(lo, hi, color=ACC, alpha=0.14)
ax.set_xlabel("bootstrap ΔPSI"); ax.set_ylabel("frequency")
ax.set_title("B  Bootstrap confidence interval", loc="left")
ax.legend(fontsize=7, frameon=False, loc="upper left")
ax.annotate(f"95% CI: [{lo:+.2f}, {hi:+.2f}]", xy=(0.98, 0.93), xycoords="axes fraction",
            fontsize=8, color="black", ha="right", va="top")

# --- 1C: birini disarida birak
ax = axes[2]
loo = core.leave_one_out_dpsi(r)
labs = [k.replace("_cikarildi", "").replace("_", " ") for k, _ in loo]
vals = [v for _, v in loo]
cols = [KD_C if k.startswith("KD") else CT_C for k, _ in loo]
y = np.arange(len(vals))
ax.barh(y, vals, color=cols, height=0.62)
ax.axvline(0, color="#333", lw=1)
ax.axvline(0.108, color="#999", lw=1, ls=":")
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=8)
ax.invert_yaxis(); ax.set_xlabel("ΔPSI with that replicate removed")
ax.set_title("C  Leave-one-out sensitivity", loc="left")


fig.savefig(os.path.join(FIG, "Figure2_TRPC1_robustness.png"))
fig.savefig(os.path.join(FIG, "Figure2_TRPC1_robustness.pdf"))
fig.savefig(os.path.join(FIG, "Figure2_TRPC1_robustness.svg"))
plt.close(fig)
print("Figure 2 (TRPC1 robustness) written")

# ======================================= FIG 7: transcript-family abundance (TPM)
# TPM adjusted for library composition (build_family_abundance.py): unadjusted TPM would
# show every gene about a quarter lower in knockdown, because a few abundant transcripts
# take a larger share of the fixed TPM total there.
S = pd.read_csv(os.path.join(core.OUT, "03_TABLOLAR", "v3", "Table4_family_TPM.csv"))
S = S[(S.Gene != "FAMILY TOTAL") & S.TPM_control_adjusted.notna()]
fams = ["ORAI (CRAC channel)", "SERCA (Ca2+ re-uptake into ER)", "STIM (ER Ca2+ sensor)",
        "SOCE regulators", "Mitochondrial Ca2+ uptake"]
TITLES = {"ORAI (CRAC channel)": "ORAI (CRAC channel)",
          "SERCA (Ca2+ re-uptake into ER)": "SERCA (Ca$^{2+}$ re-uptake into ER)",
          "STIM (ER Ca2+ sensor)": "STIM (ER Ca$^{2+}$ sensor)",
          "SOCE regulators": "SOCE regulators",
          "Mitochondrial Ca2+ uptake": "Mitochondrial Ca$^{2+}$ uptake"}
fig, _axs = plt.subplots(3, 2, figsize=(7.1, 8.2), layout="constrained")
axes = [_axs[0,0], _axs[0,1], _axs[1,0], _axs[1,1], _axs[2,0]]
_axs[2,1].axis("off")
RTQ = {"TRPC1", "STIM1", "ORAI1", "ATP2A3"}
for ax, fam in zip(axes, fams):
    sub = S[S.Family == fam].sort_values("TPM_control_adjusted", ascending=False)
    x = np.arange(len(sub)); w = 0.38
    ax.bar(x - w/2, sub.TPM_control_adjusted, w, color=CT_C, label="Control")
    ax.bar(x + w/2, sub.TPM_KD_adjusted, w, color=KD_C, label="TDP-43 KD")
    ax.set_ylim(0, 1.2 * max(sub.TPM_control_adjusted.max(), sub.TPM_KD_adjusted.max()))
    labs = [f"$\\bf{{{g}}}$*" if g in RTQ else g for g in sub.Gene]
    ax.set_xticks(x); ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=8)

    if ax is axes[0]:
        ax.set_ylabel("TPM, adjusted for library composition")
    net = sub.TPM_KD_adjusted.sum() - sub.TPM_control_adjusted.sum()
    pct = 100 * net / sub.TPM_control_adjusted.sum()
    ax.set_title(f"{TITLES[fam]}  |  net {pct:+.1f}%", loc="left", fontsize=8.5)
_axs[2,1].legend(handles=[Patch(facecolor=CT_C, label="Control"),
                    Patch(facecolor=KD_C, label="TDP-43 KD")],
           loc="center", frameon=False, fontsize=8)
fig.savefig(os.path.join(FIG, "Figure7_transcript_family_abundance.png"))
fig.savefig(os.path.join(FIG, "Figure7_transcript_family_abundance.pdf"))
plt.close(fig)
print("Figure 7 written")

# ================================================================ FIG 1: detection power
P = pd.read_csv(os.path.join(TAB, "guc_simulasyonu.tsv"), sep="\t")
fig, ax = plt.subplots(figsize=(5.8, 4.0), layout="constrained")
dps = [0.05, 0.10, 0.15, 0.20, 0.30]
cmap = ["#000000", "#0072B2", "#009E73", "#D55E00"]
for (i, row), c in zip(P.iterrows(), cmap):
    ys = [row[f"dPSI_{d:.2f}"] for d in dps]
    ax.plot(dps, ys, "o-", color=c, lw=2, ms=5, label=f"{int(row.derinlik)} reads/sample")
ax.axhline(0.8, color="#333", ls="--", lw=1)
ax.annotate("80% power", xy=(0.055, 0.82), fontsize=8)
ax.axvline(0.10, color="black", ls=":", lw=1.2)
ax.annotate("|ΔPSI| = 0.10", xy=(0.105, 0.41), fontsize=8, color="black")
ax.set_xlabel("true ΔPSI"); ax.set_ylabel("detection power")
ax.set_ylim(0, 1.02)
ax.set_title("Detection power in a 3 + 3 design", loc="left")
# lower right: the upper left holds the "80% power" label
ax.legend(fontsize=7.5, frameon=False, loc="lower right")
fig.savefig(os.path.join(FIG, "Figure1_detection_power.png"))
fig.savefig(os.path.join(FIG, "Figure1_detection_power.pdf"))
fig.savefig(os.path.join(FIG, "Figure1_detection_power.svg"))
plt.close(fig)
print("Figure 1 (detection power) written")

# ================================================================ FIG 3: saglam adaylar
B = pd.read_csv(os.path.join(TAB, "SOCE_izoform_bootstrap_GA.tsv"), sep="\t")
B = B[B.veri_seti == "GSE296712_SHSY5Y"].copy()
B["etiket"] = B.gen + " (" + B.olay + ")"
B = B.sort_values("dPSI")
fig, ax = plt.subplots(figsize=(7.1, 4.1), layout="constrained")
y = np.arange(len(B))
for i, (_, r) in enumerate(B.iterrows()):
    good = r.sifir_GA_icinde == "hayir"
    ax.plot([r.GA_alt, r.GA_ust], [i, i], color=(ACC if good else "#bbb"), lw=3,
            solid_capstyle="round")
    ax.plot(r.dPSI, i, "o", ms=8, color=(ACC if good else "#999"),
            markeredgecolor="white", markeredgewidth=1.2, zorder=3)
# TRPC1'i referans olarak ekle
ax.plot([-0.093, 0.522], [len(B), len(B)], color=KD_C, lw=3, alpha=1.0,
        solid_capstyle="round")
ax.plot(0.108, len(B), "o", ms=8, color=KD_C, markeredgecolor="white", markeredgewidth=1.2)
labels = list(B.etiket) + ["TRPC1 (SE), the event of Figure 2"]
ax.axvline(0, color="#333", lw=1, ls="--")
ax.set_yticks(list(y) + [len(B)]); ax.set_yticklabels(labels, fontsize=8.5)
ax.get_yticklabels()[-1].set_color("black")
ax.get_yticklabels()[-1].set_weight("bold")
ax.set_xlabel("ΔPSI (bootstrap 95% CI)")
ax.set_title("SOCE splicing candidates in SH-SY5Y", loc="left")
fig.savefig(os.path.join(FIG, "Figure3_robust_candidates.png"))
fig.savefig(os.path.join(FIG, "Figure3_robust_candidates.pdf"))
plt.close(fig)
print("Figure 3 written")
print("\nFigures ->", FIG)

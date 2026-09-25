"""Figure 2: keep panels A-D unchanged and add panel E, both recordings on one scale.

Panels A-D are the existing figure (source_data/fura2_traces/Figure2_panels_A-D.png).
Panel E replots the original ratio/time exports of the two representative recordings
(PRISM_RAW_CONTROL_130626.txt, PRISM_RAW_KD_130626.txt; 1-s sampling) on identical y and
time axes. Time is expressed relative to the steepest point of the Ca2+-readdition rise,
located on a 15-s moving average; no values are smoothed in the plotted traces.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data" / "fura2_traces"
BLUE, ORANGE, INK, MUTED = "#0072B2", "#D55E00", "#222222", "#6B6B6B"
plt.rcParams.update({"font.family": "Arial", "font.size": 7.5, "axes.linewidth": 0.6,
                     "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "savefig.facecolor": "white",
                     "svg.fonttype": "none", "mathtext.fontset": "custom", "mathtext.rm": "Arial"})


def load(name):
    d = np.loadtxt(SRC / f"PRISM_RAW_{name}_130626.txt")
    t, r = d[:, 0], d[:, 1]
    assert np.allclose(np.median(np.diff(t)), 1.0)
    k = 15
    smooth = np.convolve(r, np.ones(k) / k, mode="same")
    slope = np.gradient(smooth, t)
    onset = t[k + np.argmax(slope[k:-k])]
    return t - onset, r, onset


ctrl_t, ctrl_r, ctrl_on = load("CONTROL")
kd_t, kd_r, kd_on = load("KD")
assert 2600 < ctrl_on < 2700 and 1850 < kd_on < 1950

img = mpimg.imread(SRC / "Figure2_panels_A-D.png")
h, w = img.shape[:2]
width = 6.05
top_h = width * h / w
fig = plt.figure(figsize=(width, top_h + 1.95))
ax_img = fig.add_axes([0, 1.95 / (top_h + 1.95), 1, top_h / (top_h + 1.95)])
ax_img.imshow(img, interpolation="lanczos"); ax_img.axis("off")

ax = fig.add_axes([0.105, 0.045, 0.86, 1.45 / (top_h + 1.95)])
ax.plot(ctrl_t, ctrl_r, color=BLUE, lw=0.55, label="Non-targeting shRNA")
ax.plot(kd_t, kd_r, color=ORANGE, lw=0.55, label="shTDP-43")
ax.axvline(0, color=MUTED, lw=0.6, ls=(0, (3, 2)))
ax.text(18, 2.62, "Ca$^{2+}$ readdition", fontsize=6.8, color=MUTED, va="top")
ax.set_xlim(-1300, 950); ax.set_ylim(0.6, 2.7)
ax.set_xlabel("Time relative to the readdition rise (s)")
ax.set_ylabel("F$_{340}$/F$_{380}$")
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(direction="out", length=2.5, width=0.6)
ax.legend(frameon=False, loc="upper left", fontsize=6.8, handlelength=1.6)
fig.text(0.012, 1.95 / (top_h + 1.95) - 0.004, "E", weight="bold", fontsize=11, va="top")
for ext in ("png", "svg"):
    fig.savefig(ROOT / "figures" / "main" / f"Figure2_calcium_responses.{ext}", dpi=400)
print(f"readdition rise: control {ctrl_on:.1f} s, shTDP-43 {kd_on:.1f} s")

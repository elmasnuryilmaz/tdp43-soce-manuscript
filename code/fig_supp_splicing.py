#!/usr/bin/env python3
"""Supplementary Figure S3: per-library PSI for four SOCE splicing events."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
PKG=Path("/Users/elmas/Desktop/MAKALE/11_NEUROCHEMISTRY_INTERNATIONAL_FIGURE_REVISION")
OUT=PKG/"figures"/"supplementary"; OUT.mkdir(parents=True,exist_ok=True)
BLUE,ORANGE="#0072B2","#D55E00"
plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.labelsize":8,
 "axes.titlesize":9,"axes.titleweight":"bold","xtick.labelsize":7.5,
 "ytick.labelsize":7.5,"axes.spines.top":False,"axes.spines.right":False,
 "pdf.fonttype":42,"ps.fonttype":42,"savefig.dpi":300,
 "mathtext.fontset":"custom","mathtext.rm":"Arial","mathtext.it":"Arial:italic",
 "mathtext.bf":"Arial:bold","figure.facecolor":"white","savefig.facecolor":"white"})
t=pd.read_csv(PKG/"tables"/"Table2_robust_SOCE_splicing_events.csv").set_index("gene")
order=["STIMATE","ORAI3","STIM2","STIM1"]
fig,axes=plt.subplots(2,2,figsize=(7.1,6.0),layout="constrained")
for label,ax,g in zip("ABCD",axes.flat,order):
 r=t.loc[g]
 c=np.array([float(x) for x in r.PSI_control.split(";")])
 k=np.array([float(x) for x in r.PSI_knockdown.split(";")])
 for x,v,color in [(0,c,BLUE),(1,k,ORANGE)]:
  ax.scatter(x+np.linspace(-.10,.10,len(v)),v,s=35,color=color,zorder=3)
  ax.hlines(v.mean(),x-.20,x+.20,color=color,lw=1.5,zorder=2)
 ax.set_ylim(-.04,1.06);ax.set_xlim(-.45,1.45)
 ax.set_xticks([0,1],["Control","TDP-43 depletion"])
 ax.set_ylabel("PSI")
 ax.set_title(f"{label}  {g} skipped exon",loc="left")
 ax.text(.98,.96,f"rMATS ΔPSI {r.delta_PSI:+.3f}\n95% CI {r.CI95_low:+.3f} to {r.CI95_high:+.3f}",
         transform=ax.transAxes,ha="right",va="top",fontsize=7.5,color="black")
fig.savefig(OUT/"Supplementary_Figure_S3_splicing_replicates.pdf")
fig.savefig(OUT/"Supplementary_Figure_S3_splicing_replicates.png")
print("S3 written")

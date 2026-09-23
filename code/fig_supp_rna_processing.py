#!/usr/bin/env python3
"""Supplementary Figures S4, S5, S6 and S8 from the analysis source tables."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

PKG = Path("/Users/elmas/Desktop/MAKALE/11_NEUROCHEMISTRY_INTERNATIONAL_FIGURE_REVISION")
OUT = PKG / "figures" / "supplementary"
OUT.mkdir(parents=True, exist_ok=True)
D = Path("/Users/elmas/Desktop/MAKALE/07_DISK_ANALIZLERI")
RESULTS = D / "sonuclar"
PANEL = Path("/Users/elmas/Desktop/TEZ/output/calcium_panel_overlap_analysis_2026-04-25/01_gene_panels")
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#999999"
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.labelsize": 8,
    "axes.titlesize": 9, "axes.titleweight": "bold",
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .7, "xtick.major.width": .7, "ytick.major.width": .7,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "savefig.dpi": 300, "figure.facecolor": "white", "savefig.facecolor": "white",
})
def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png")
    plt.close(fig)

# S4. Working cryptic-junction positive controls, cumulative Ca-gene panels,
# and the RNA-binding-protein specificity comparison.
kr = pd.read_csv(RESULTS / "YUKSEK_GUVEN_OWN_SH_SY5Y.tsv", sep="\t")
pos = ["STMN2","ACTL6B","UNC13A","PFKP","HDGFL2","AGRN","ATG4B",
       "ELAVL3","SETD5","ARHGAP32","GPSM2","RSF1"]
d = kr[kr.gene.isin(pos)].sort_values("dPSI").drop_duplicates("gene", keep="last")
fig, axes = plt.subplots(3,1,figsize=(7.1,8.8),layout="constrained",
                        gridspec_kw={"height_ratios":[1.30,.80,1.0]})
ax=axes[0]
ax.barh(d.gene,d.dPSI,color=BLUE,height=.7)
for i,(v,k,c) in enumerate(zip(d.dPSI,d.okuma_KD,d.okuma_CTRL)):
    ax.text(v+.014,i,f"{int(k)}/{int(c)}",va="center",fontsize=7,color="black")
ax.set_xlim(0,1.32);ax.set_xlabel("ΔPSI (knockdown − control)")
ax.set_title("A  Literature positive controls recovered de novo",loc="left")
names=["Tier1_SOCE_TRP_51","Tier2_Channel_Release_Transport_117",
       "Tier3_Curated_Calcium_Handling_258","Tier4_Expanded_Calcium_Associated_732"]
sizes=[51,117,258,732]
counts=[]
for name in names:
    members=set(pd.read_csv(PANEL/f"{name}.csv").gene_upper)
    counts.append(kr[kr.gene.isin(members)].gene.nunique())
ax=axes[1]
x=np.arange(4)
ax.bar(x,counts,color=BLUE,width=.6)
for i,n in enumerate(counts):
    ax.text(i,n+.08,str(n),ha="center",fontsize=8,color="black")
ax.set_xticks(x,[f"Tier {i+1}\n({sizes[i]} genes)" for i in range(4)])
ax.set_ylabel("Genes with a cryptic event")
ax.set_ylim(0,max(counts)+1.2)
ax.set_title("B  Cumulative Ca$^{2+}$ gene panels",loc="left")
posgenes={"STMN2","UNC13A","HDGFL2","ACTL6B","AGRN","KALRN","ARHGAP32","PFKP",
          "ATG4B","SETD5","ELAVL3","POLDIP3","CAMK2B","RSF1","GPSM2","SYNJ2"}
order=["SH_SY5Y","SH_SY5Y_DOZ25","iPSC_koloni","K562_mRNA","iPSC_MN",
       "K562_totalRNA","iPSC_MN_FUS","iPSC_MN_TAF15"]
label={"SH_SY5Y":"SH-SY5Y 75 ng/mL","SH_SY5Y_DOZ25":"SH-SY5Y 25 ng/mL",
       "iPSC_koloni":"iPSC colonies","K562_mRNA":"K562 poly(A)+",
       "iPSC_MN":"iPSC-MN · TDP-43","K562_totalRNA":"K562 total RNA",
       "iPSC_MN_FUS":"iPSC-MN · FUS","iPSC_MN_TAF15":"iPSC-MN · TAF15"}
values=[];labels=[];cols=[]
for key in order:
    f=RESULTS/f"RT_KRIPTIK_{key}.tsv"
    if not f.exists(): continue
    found=set(pd.read_csv(f,sep="\t").gene.astype(str).str.upper())
    values.append(len(found & posgenes));labels.append(label[key])
    cols.append(GREY if key.endswith(("FUS","TAF15")) else
                (ORANGE if key.startswith("SH_SY5Y") else BLUE))
ax=axes[2]
yy=np.arange(len(values))[::-1]
ax.barh(yy,values,color=cols,height=.7)
for ypos,v in zip(yy,values):
    ax.text(v+.2,ypos,str(v),va="center",fontsize=7.5,color="black")
ax.set_yticks(yy,labels);ax.set_xlim(0,18)
ax.set_xlabel("Positive-control genes recovered (of 16)")
ax.set_title("C  Human comparison datasets (permissive definition)",loc="left")
save(fig,"Supplementary_Figure_S4_cryptic_controls")

# S5. Dependence-aware exploratory NMD interaction.
nmd=pd.read_csv(RESULTS/"NMD_etkilesim_paylasimli_kontrol_t4.tsv",sep="\t")
old=pd.read_csv(RESULTS/"NMD_etkilesim_tum_genler.tsv",sep="\t",usecols=["ens","sembol"])
nmd=nmd.merge(old,on="ens",how="left")
nmd["symbol"]=nmd.sembol.fillna(nmd.ens)
cryptic={g for g in kr.gene.astype(str) if not g.startswith("ENSG") and g!="."}
background=nmd.loc[~nmd.symbol.isin(cryptic),"interaction_log2"]
on=nmd.loc[nmd.symbol.isin(cryptic),"interaction_log2"]
cb=float(nmd.loc[nmd.symbol=="CBARP","interaction_log2"].iloc[0])
fig,ax=plt.subplots(figsize=(7.1,4.0),layout="constrained")
bins=np.linspace(-2,2,60)
ax.hist(background,bins=bins,density=True,color="#C9CDCF",
        label=f"Other genes (n={len(background):,})")
ax.hist(on,bins=bins,density=True,histtype="step",lw=1.8,color=ORANGE,
        label=f"Cryptic-junction genes (n={len(on)})")
ax.axvline(cb,color=BLUE,lw=1.4)
ax.text(cb-.04,1.16,f"CBARP {cb:+.2f}",ha="right",va="top",color="black",fontsize=8)
ax.set_ylim(0,1.40)
ax.set_xlabel("TDP-43-specific NMD interaction (log$_2$)")
ax.set_ylabel("Density")
ax.set_title("Exploratory NMD interaction",loc="left")
ax.legend(frameon=False,fontsize=7.5,loc="upper left")
save(fig,"Supplementary_Figure_S5_NMD_interaction")

# S6. Coverage-index candidates with their existing bootstrap intervals,
# plus the STMN2 coverage positive control.
a=pd.read_csv(RESULTS/"APA_corrected_full_core_summary.tsv",sep="\t")
a=a[a.delta.abs()>=.05].sort_values("delta").reset_index(drop=True)
fig,(ax,ax2)=plt.subplots(2,1,figsize=(7.1,9.0),layout="constrained",
                         gridspec_kw={"height_ratios":[1.8,1.0]})
for i,r in a.iterrows():
    color=GREEN if r.measure=="IPA_index" else BLUE
    ax.plot([r.boot_low,r.boot_high],[i,i],color=color,lw=1.3)
    ax.plot(r.delta,i,"o",color=color,ms=4.8)
ax.axvline(0,color="#444444",lw=.8)
lab=[f"{r.gene} · {r.unit.replace('termexon','terminal exon')}"
     +(" †" if (r.gene,r.unit)==("STIM2","intron13") else "")
     for _,r in a.iterrows()]
ax.set_yticks(range(len(a)),lab,fontsize=7)
ax.set_xlabel("Δ coverage index (knockdown − control)")
ax.set_title("A  Depth-qualified coverage-index changes",loc="left")
ax.legend(handles=[Patch(color=GREEN,label="Intronic index"),
                   Patch(color=BLUE,label="Distal 3′UTR index")],
          frameon=False,fontsize=7.5,loc="upper left")
man=pd.read_csv(D/"kod"/"ornekler.tsv",sep="\t")
sub=man[man.dataset=="SH_SY5Y"]
kd=list(sub[sub.group=="KD"]["sample"]);ct=list(sub[sub.group=="CTRL"]["sample"])
b=pd.read_csv(RESULTS/"apa_bedcov_SH_SY5Y_corrected_full.tsv",sep="\t",header=None,
              names=["chrom","start","end","isim"]+kd+ct)
st=b[b.isim.str.startswith("STMN2|intron2|")]
assert len(st)==2
i5=st[st.isim.str.contains("I5")].iloc[0]
i3=st[st.isim.str.contains("I3")].iloc[0]
for j in range(3):
    x=np.array([0,1])
    c=[float(i5[ct[j]])/(i5["end"]-i5["start"]),
       float(i3[ct[j]])/(i3["end"]-i3["start"])]
    k=[float(i5[kd[j]])/(i5["end"]-i5["start"]),
       float(i3[kd[j]])/(i3["end"]-i3["start"])]
    ax2.plot(x,c,"-o",color=BLUE,ms=4.5,label="Control" if j==0 else None)
    ax2.plot(x+.07,k,"-o",color=ORANGE,ms=4.5,label="TDP-43 depletion" if j==0 else None)
    if k[1]==0:
        ax2.plot(1.07,0,"o",mfc="white",mec=ORANGE,mew=1.4,ms=8,zorder=5)
        ax2.annotate("Measured zero",(1.07,0),xytext=(8,2),textcoords="offset points",
                     color="black",fontsize=7)
ax2.set_xlim(-.12,1.42)
ax2.set_xticks([.035,1.035],["Intron 5′ end","Intron 3′ end"])
ax2.set_ylabel("Mean coverage depth")
ax2.set_title("B  STMN2 intron 2 positive control",loc="left")
ax2.legend(frameon=False,fontsize=7.5)
save(fig,"Supplementary_Figure_S6_APA_coverage")

# S8. Cryptic STMN2 frequency and region-level correlations of two measures.
s=pd.read_csv(RESULTS/"NYGC_kriptik_PSI_ornek_duzeyi.tsv",sep="\t",low_memory=False)
s=s[s.STMN2_toplam>=20].copy()
s["control"]=s.grup.astype(str).str.contains("Non-Neurological")
s["als"]=s.grup.astype(str).str.strip().eq("ALS Spectrum MND")
longregions=["Spinal Cord Lumbar","Spinal Cord Cervical","Spinal Cord Thoracic",
 "Cortex Motor Unspecified","Cortex Motor Lateral","Cortex Motor Medial",
 "Cortex Temporal","Hippocampus","Cortex Frontal","Cortex Occipital","Cerebellum"]
short=["Lumbar cord","Cervical cord","Thoracic cord","Motor cortex, other",
       "Motor cortex, lateral","Motor cortex, medial","Temporal cortex",
       "Hippocampus","Frontal cortex","Occipital cortex","Cerebellum"]
vals=[];labels=[]
for raw,lab in zip(longregions,short):
    als=s[s.als & s.doku.eq(raw)]
    control=s[s.control & s.doku.eq(raw)]
    if len(als)<20:continue
    vals.append((100*(als.STMN2_kriptik>0).mean(),
                 100*(control.STMN2_kriptik>0).mean() if len(control) else 0))
    labels.append(f"{lab} (ALS {len(als)}; control {len(control)})")
regions=[r for r in short if any(x.startswith(r+" (") for x in labels)]
corr=pd.read_csv(PKG/"supplementary"/"S9_cryptic_PSI_correlations_within_ALS.csv")
genes=["SNAP25","TRPC1","SARAF","ATP2A2","CBARP"]
proxies=["gene-level STMN2","cryptic STMN2 PSI"]
fig=plt.figure(figsize=(7.1,8.8),layout="constrained")
gs=fig.add_gridspec(2,2,height_ratios=[1.15,1.0])
ax=fig.add_subplot(gs[0,:])
y=np.arange(len(vals));w=.36
ax.barh(y-w/2,[v[0] for v in vals],w,color=ORANGE,label="ALS")
ax.barh(y+w/2,[v[1] for v in vals],w,color=BLUE,label="Non-neurological control")
ax.set_yticks(y,labels,fontsize=7)
ax.invert_yaxis();ax.set_xlim(0,80)
ax.set_xlabel("Samples with cryptic STMN2 junction (%)")
ax.set_title("A  Cryptic STMN2 by region",loc="left")
ax.legend(frameon=False,fontsize=7.5,loc="lower right")
cmap=LinearSegmentedColormap.from_list("rho",
                                      ["#95C4DF","#FFFFFF","#F1B481"])
for j,proxy in enumerate(proxies):
    ax=fig.add_subplot(gs[1,j])
    mat=corr[(corr.proxy==proxy)&corr.target_gene.isin(genes)].pivot(
        index="region",columns="target_gene",values="spearman_rho")
    qmat=corr[(corr.proxy==proxy)&corr.target_gene.isin(genes)].pivot(
        index="region",columns="target_gene",values="q_value")
    # Full, shared eleven-region order; the cohort table uses the short names.
    reg_order=[r for r in longregions if r in mat.index]
    mat=mat.reindex(reg_order)[genes];qmat=qmat.reindex(reg_order)[genes]
    im=ax.imshow(mat.to_numpy(),vmin=-1,vmax=1,cmap=cmap,aspect="auto")
    for iy in range(len(reg_order)):
        for ix in range(len(genes)):
            if pd.notna(qmat.iat[iy,ix]) and qmat.iat[iy,ix]<.05:
                ax.plot(ix,iy,"o",ms=2.8,color="black")
    ax.set_xticks(range(len(genes)),genes,rotation=45,ha="right",fontsize=7)
    ax.set_yticks(range(len(reg_order)),[short[longregions.index(r)] for r in reg_order] if j==0 else [],fontsize=6.7)
    ax.set_title(f"{'B' if j==0 else 'C'}  {proxy}",loc="left")
cb=fig.colorbar(im,ax=fig.axes[1:3],fraction=.025,pad=.025,ticks=[-1,0,1])
cb.set_label("Spearman ρ")
save(fig,"Supplementary_Figure_S8_STMN2_proxies")
print("S4, S5, S6, S8 written")

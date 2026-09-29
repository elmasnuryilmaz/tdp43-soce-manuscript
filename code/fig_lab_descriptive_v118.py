"""Current Figure 1: biological RT-qPCR replicates and within-experiment WST-1 wells."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'supplementary/S1_laboratory_source_data.xlsx'
td=pd.read_excel(source,sheet_name='TARDBP_qPCR')
qp=pd.read_excel(source,sheet_name='Target_qPCR_rel')
ws=pd.read_excel(source,sheet_name='WST1')
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':10,
    'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,
    'svg.fonttype':'none','savefig.dpi':300,'savefig.facecolor':'white'})
blue,orange,grey='#0072B2','#D55E00','#999999'
fig,axes=plt.subplots(1,3,figsize=(9.0,3.8),gridspec_kw={'width_ratios':[1.1,1.65,1]},layout='constrained')
def draw(ax,x,v,color,width=.55):
    v=np.asarray(v,float)
    ax.bar(x,v.mean(),width,color=color,alpha=.9)
    ax.errorbar(x,v.mean(),yerr=v.std(ddof=1)/np.sqrt(len(v)),color='#333333',capsize=3)
    ax.scatter(x+np.linspace(-.09,.09,len(v)),v,facecolors='white',edgecolors='#222222',s=23,zorder=3)

ax=axes[0]
for i,(g,c) in enumerate(zip(['Untransduced control','Non-targeting shRNA control','shTDP-43'],[grey,blue,orange])):
    draw(ax,i,td.loc[td.group==g,'rel_expression'],c)
ax.set_xticks([0,1,2],['No virus','Control\nshRNA','shTDP-43'])
ax.set_ylabel('TARDBP mRNA (2$^{-\\Delta\\Delta Ct}$)')
ax.set_ylim(0,1.3); ax.set_title('A TARDBP mRNA',loc='left')
ax.text(.5,.97,'4 biological replicates/group',transform=ax.transAxes,ha='center',va='top',fontsize=7.5)
ax=axes[1]
genes=['TRPC1','STIM1','ORAI1','ATP2A3']
for i,g in enumerate(genes):
    for off,grp,c in [(-.18,'Non-targeting shRNA control',blue),(.18,'shTDP-43',orange)]:
        draw(ax,i+off,qp.loc[(qp.gene==g)&(qp.group==grp),'rel_expression'],c,.35)
ax.set_xticks(range(4),genes); ax.set_ylim(0,4.8)
ax.set_ylabel('Relative mRNA (2$^{-\\Delta\\Delta Ct}$)')
ax.set_title('B Calcium-related mRNAs',loc='left')
ax.legend([plt.Rectangle((0,0),1,1,color=blue),plt.Rectangle((0,0),1,1,color=orange)],['Non-targeting shRNA','shTDP-43'],frameon=False,fontsize=7.5,loc='upper left')
ax.text(.5,.77,'4 biological replicates/group',transform=ax.transAxes,ha='center',fontsize=7.5)
ax=axes[2]
for i,(grp,c) in enumerate([('Non-targeting shRNA control',blue),('shTDP-43',orange)]):
    draw(ax,i,ws.loc[ws.group==grp,'signal_pct_of_control'],c)
ax.axhline(100,color=grey,ls='--',lw=.9); ax.set_ylim(0,125)
ax.set_xticks([0,1],['Non-targeting\nshRNA','shTDP-43'])
ax.set_ylabel('WST-1 signal (% of control)'); ax.set_title('C WST-1 at 48 h',loc='left')
ax.text(.5,.96,'4 wells; one experiment',transform=ax.transAxes,ha='center',va='top',fontsize=7.5)
fig.supxlabel('Mean ± SEM; A–B: biological replicates; C: wells from one experiment',fontsize=8)
for ext in ('png','svg'): fig.savefig(ROOT/f'figures/main/Figure1_functional_consequences.{ext}')
plt.close(fig)

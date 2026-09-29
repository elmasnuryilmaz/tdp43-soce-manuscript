"""Current Supplementary Figure S5 from the descriptive NMD table."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
nmd=pd.read_csv(ROOT/'source_data/nmd_descriptive_all_genes.csv.gz')
kr=pd.read_csv(ROOT/'supplementary/S5_high_confidence_cryptic_events.csv')
cryptic=set(kr.loc[kr.comparison=='SH-SY5Y 75 ng/mL (MAPQ-filtered set)','gene'])
if not cryptic:
    raw=pd.read_csv(ROOT.parent/'07_DISK_ANALIZLERI/sonuclar/YUKSEK_GUVEN_OWN_SH_SY5Y.tsv',sep='\t')
    cryptic=set(raw.gene)
cryptic={g for g in cryptic if not str(g).startswith('ENSG') and g!='.'}
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':10,'axes.titleweight':'bold',
    'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','savefig.dpi':300})
fig,ax=plt.subplots(figsize=(7.1,4.0),layout='constrained')
background=nmd.loc[~nmd.gene.isin(cryptic),'interaction_log2']; on=nmd.loc[nmd.gene.isin(cryptic),'interaction_log2']
bins=np.linspace(-2,2,60)
ax.hist(background,bins=bins,density=True,color='#C9CDCF',label=f'Other genes (n={len(background):,})')
ax.hist(on,bins=bins,density=True,histtype='step',lw=1.8,color='#D55E00',label=f'Cryptic-junction genes (n={len(on)})')
cb=float(nmd.loc[nmd.gene=='CBARP','interaction_log2'].iloc[0]); ax.axvline(cb,color='#0072B2',lw=1.4)
ax.text(cb-.04,1.16,f'CBARP {cb:+.2f}',ha='right',va='top')
ax.set_ylim(0,1.4); ax.set_xlabel('Mean TDP-43 × NMD-inhibition interaction (log$_2$)'); ax.set_ylabel('Density')
ax.set_title('Descriptive NMD interaction',loc='left'); ax.legend(frameon=False,loc='upper left')
fig.supxlabel('Four contrasts share baseline libraries; no inferential significance assigned',fontsize=8)
for ext in ('png','svg'): fig.savefig(ROOT/f'figures/supplementary/Supplementary_Figure_S5_NMD_interaction.{ext}')
plt.close(fig)

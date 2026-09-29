"""Rebuild corrected S13 and descriptive NMD tables, and expose STIM1 eligibility.

NMD point estimates retain the original normalization, expression filter and four
contrasts. No test treats the shared-baseline contrasts as independent replicates.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DISK = ROOT.parent / '07_DISK_ANALIZLERI'
SUP = ROOT / 'supplementary'
SRC = ROOT / 'source_data'

def s13():
    names = {'SH_SY5Y':'SH-SY5Y 75 ng/mL','SH_SY5Y_DOZ25':'SH-SY5Y 25 ng/mL',
             'iPSC_koloni':'iPSC colonies','iPSC_MN':'iPSC-MN, TDP-43 KD',
             'iPSC_MN_FUS':'iPSC-MN, FUS KD','iPSC_MN_TAF15':'iPSC-MN, TAF15 KD',
             'K562_totalRNA':'K562 total RNA','K562_mRNA':'K562 poly(A)+ mRNA',
             'Fare_striatum':'Mouse striatum'}
    cols = {'sinif':'junction_class','lsv_tipi':'lsv_type','psi_KD':'PSI_knockdown',
            'psi_CTRL':'PSI_control','dPSI':'delta_PSI','GA_alt':'CI95_low','GA_ust':'CI95_high',
            'q':'q_value','okuma_KD':'reads_knockdown','okuma_CTRL':'reads_control',
            'toplam_KD':'LSV_total_knockdown','toplam_CTRL':'LSV_total_control'}
    out=[]
    for f in sorted((DISK/'sonuclar').glob('RT_LSV_*.tsv.gz')):
        ds=f.name[7:-7]
        mouse=ds in ('C2C12','NSC34','Fare_striatum')
        chrom,start,end,gene=('5',54110115,54110139,'Stim2') if mouse else ('4',27007982,27008006,'STIM2')
        d=pd.read_csv(f,sep='\t',low_memory=False)
        d['chrom']=d.chrom.astype(str).str.removeprefix('chr')
        d=d[(d.gene==gene)&(d.chrom==chrom)]
        d=d[(d.end==start)|(d.start==end+1)].copy()
        d['flank']=np.where(d.end==start,'upstream','downstream')
        d['comparison']=names.get(ds,ds)
        d=d.rename(columns=cols)
        d['lsv_type']=d.lsv_type.replace({'verici':'donor','alici':'acceptor'})
        d['junction_class']=d.junction_class.replace({'anotasyonlu':'annotated','yeni_kombinasyon':'novel combination','yeni_bolge':'novel splice site','tamamen_yeni':'fully novel'})
        out.append(d[['comparison','flank','chrom','start','end','junction_class','lsv_type',*list(cols.values())[2:]]])
    result=pd.concat(out,ignore_index=True)
    result.to_csv(SUP/'S13_STIM2_SOAR_exon_junction_level.csv',index=False)
    print('S13:',len(result),'rows;',result.flank.value_counts().to_dict())

def nmd():
    base=DISK/'nmd_GSE307054'
    counts=pd.read_csv(base/'GSE307054_counts.csv.gz',index_col=0)
    meta=pd.read_csv(base/'GSE307054_meta.csv.gz',index_col=0)
    cols=counts.columns.intersection(meta.index)
    norm=counts[cols].div(meta.loc[cols,'sizeFactor'],axis=1)
    lg=np.log2(norm+1)
    conds=['X','XS','XU','US']
    baseline=lg[['CON_1','CON_2']].mean(axis=1)-lg[['TDPKD_1','TDPKD_2']].mean(axis=1)
    values={c:lg[[f'T{c}_1',f'T{c}_2']].mean(axis=1)-lg[[f'C{c}_1',f'C{c}_2']].mean(axis=1)+baseline for c in conds}
    d=pd.DataFrame(values)
    keep=(norm[[f'T{c}_{i}' for c in conds for i in (1,2)]].mean(axis=1)>=10)&(norm[['TDPKD_1','TDPKD_2']].mean(axis=1)>=5)
    d=d[keep].copy()
    old=pd.read_csv(DISK/'sonuclar/NMD_etkilesim_tum_genler.tsv',sep='\t',usecols=['ens','sembol']).drop_duplicates('ens').set_index('ens')
    d.insert(0,'gene',old.sembol.reindex(d.index).fillna(pd.Series(d.index,index=d.index)))
    d.insert(1,'interaction_log2',d[conds].mean(axis=1))
    d.insert(2,'n_conditions_positive',(d[conds]>0).sum(axis=1))
    d.insert(3,'min_condition_interaction',d[conds].min(axis=1))
    d.insert(4,'max_condition_interaction',d[conds].max(axis=1))
    d=d.rename(columns={'X':'interaction_XRN1','XS':'interaction_XRN1_SMG6','XU':'interaction_XRN1_UPF1','US':'interaction_UPF1_SMG6'})
    d.index.name='ensembl_gene_id'
    d=d.reset_index().sort_values('interaction_log2',ascending=False)
    d.to_csv(SRC/'nmd_descriptive_all_genes.csv.gz',index=False)
    previous=pd.read_csv(SUP/'S10_NMD_interaction_SOCE_panel.csv')
    selected=d[d.gene.isin(previous.gene)]
    selected.to_csv(SUP/'S10_NMD_interaction_SOCE_panel.csv',index=False)
    panels=pd.read_csv(SUP/'S2_calcium_gene_panels.csv')
    # Read the published long-form panel definitions, not precomputed p values.
    print('panel columns:',list(panels.columns))
    pcol=next(c for c in panels if 'panel' in c.lower() or 'tier' in c.lower())
    gcol=next(c for c in panels if c.lower() in ('gene','symbol','gene_symbol'))
    rows=[]
    for panel,g in panels.groupby(pcol,sort=False):
        z=d[d.gene.isin(g[gcol])].interaction_log2
        rows.append(dict(panel=panel,n_genes=len(z),median_interaction_log2=z.median(),interpretation='descriptive; not an NMD-sensitivity test'))
    pos={'STMN2','UNC13A','HDGFL2','ACTL6B','AGRN','KALRN','ARHGAP32','PFKP','ATG4B','SETD5','CAMK2B','ELAVL3','POLDIP3','RSF1','GPSM2','SYNJ2'}
    z=d[d.gene.isin(pos)].interaction_log2
    rows.append(dict(panel='Cryptic_splicing_reference_genes_16',n_genes=len(z),median_interaction_log2=z.median(),interpretation='cryptic-splicing references; not a validated NMD-positive panel'))
    summary=pd.DataFrame(rows)
    summary.to_csv(SUP/'S10b_NMD_panel_descriptive_summary.csv',index=False)
    for oldfile in [SUP/'S10b_NMD_panel_level_tests.csv',SRC/'nmd_panel_t4.csv']:
        if oldfile.exists(): oldfile.unlink()  # obsolete derived tests; raw measurements retained
    print('NMD:',len(d),'genes; no inferential p or q values')
    print(selected[selected.gene.isin(['CBARP','STIM1'])].to_string(index=False))

def stim1():
    base=Path('/Users/elmas/Desktop/TEZ/output/reanalysis_corrected_full_2026-07-22/isoform_0_vs_75/drimseq_stager')
    f=pd.read_csv(base/'drimseq_feature_results.csv')
    s=pd.read_csv(base/'stageR_adjusted_pvalues.csv')
    print('STIM1 input columns',list(f.columns),list(s.columns))
    f=f[f.gene_id=='STIM1'].copy()
    s=s[s.geneID=='STIM1']
    f=f.merge(s[['txID','transcript']],left_on='feature_id',right_on='txID',how='left')
    f=f.rename(columns={'transcript':'stageR_confirmation_adjusted_p'})
    f['stageR_eligible']=f.pvalue.notna()
    f['stageR_status']=np.where(f.stageR_eligible,'evaluated; not confirmed','not evaluated; missing DRIMSeq p value')
    f.to_csv(SRC/'STIM1_transcript_test_eligibility.csv',index=False)
    print(f[['feature_id','pvalue','stageR_status']].to_string(index=False))

if __name__=='__main__':
    s13(); nmd(); stim1()

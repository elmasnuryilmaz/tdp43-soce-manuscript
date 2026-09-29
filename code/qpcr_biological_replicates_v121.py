"""RT-qPCR inference after author confirmation of four biological replicates/group.
Two-sided Welch tests on Delta Ct; Holm correction within the four targets and,
separately, the two TARDBP knockdown-versus-control comparisons.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
ROOT=Path(__file__).resolve().parents[1]
def holm(p):
    p=np.asarray(p,float); order=np.argsort(p); result=np.empty(len(p))
    result[order]=np.minimum(1,np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
    return result

def calculate():
    source=ROOT/'supplementary/S1_laboratory_source_data.xlsx'
    td=pd.read_excel(source,sheet_name='TARDBP_qPCR')
    tg=pd.read_excel(source,sheet_name='Target_qPCR_Ct')
    rows=[]
    for family,comparisons in [('TARDBP',[(td,'TARDBP',g) for g in ['Untransduced control','Non-targeting shRNA control']]),('four targets',[(tg[tg.gene==g],g,'Non-targeting shRNA control') for g in ['TRPC1','STIM1','ORAI1','ATP2A3']])]:
        start=len(rows)
        for z,g,control in comparisons:
            c=z.loc[z.group==control,'dCt'].to_numpy(float); k=z.loc[z.group=='shTDP-43','dCt'].to_numpy(float)
            assert len(c)==len(k)==4
            test=ttest_ind(k,c,equal_var=False)
            rows.append(dict(gene=g,control=control,n_control=len(c),n_knockdown=len(k),experimental_unit='biological replicate',test='two-sided Welch t-test on Delta Ct',correction_family=family,mean_dCt_control=c.mean(),mean_dCt_knockdown=k.mean(),t_statistic=test.statistic,df=test.df,p_value=test.pvalue))
        for row,p in zip(rows[start:],holm([r['p_value'] for r in rows[start:]])):row['holm_adjusted_p']=p
    return pd.DataFrame(rows)
if __name__=='__main__':
    d=calculate(); d.to_csv(ROOT/'source_data/qpcr_biological_replicate_tests.csv',index=False)
    print(d[['gene','control','p_value','holm_adjusted_p']].to_string(index=False))

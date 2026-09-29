"""Donor-level sensitivity of Table 5 and S18 expression comparisons.

Keep the original per-comparison normalization; average normalized expression
within donor/region before fitting marker adjustments or comparing groups.
BH families are the regions within each group/gene/adjustment combination.
"""
import pandas as pd
from scipy import stats
import nygc_composition_and_cryptic_s18 as s

rows=[]
for group,regions in [(s.ALS,s.BRAIN+s.CORD),(s.OND,s.OND_REGIONS)]:
    for region in regions:
        result=s.comparison(region,group)
        if result is None: continue
        lg,case,ctrl=result
        ids=s.meta.loc[lg.columns,'denek']
        assert ids.notna().all()
        labels=s.meta.loc[lg.columns,['denek','grup']].drop_duplicates()
        assert not labels.denek.duplicated().any(), 'donor in both groups within region'
        labels=labels.set_index('denek').grup
        donor=lg.T.groupby(ids).mean().T
        ca=[x for x in donor.columns if labels[x]==group]
        co=[x for x in donor.columns if labels[x]==s.CONTROL]
        for gene in ('TRPC1','SARAF','CBARP'):
            models=s.MODELS if gene=='TRPC1' and region in s.BRAIN else {'none':[]}
            for model,markers in models.items():
                v=pd.Series(s.residual(donor.loc[gene].values,[donor.loc[m].values for m in markers]),index=donor.columns)
                rows.append(dict(group=s.LABEL[group],region=region,gene=gene,adjustment=model,
                    n_case_samples=len(case),n_control_samples=len(ctrl),n_case_donors=len(ca),n_control_donors=len(co),
                    cliffs_delta=s.cliff(v[ca],v[co]),p_value=stats.mannwhitneyu(v[ca],v[co],alternative='two-sided').pvalue))
d=pd.DataFrame(rows)
for _,g in d.groupby(['group','gene','adjustment']): d.loc[g.index,'q_value']=s.bh(g.p_value.values)
d.to_csv(s.ROOT/'supplementary/S18d_NYGC_donor_level_sensitivity.csv',index=False)
print('Donor level TRPC1:')
print(d[(d.gene=='TRPC1')&(d.adjustment=='none')].to_string(index=False))
print('Repeated donor comparisons:')
print(d[(d.n_case_samples!=d.n_case_donors)|(d.n_control_samples!=d.n_control_donors)].to_string(index=False))

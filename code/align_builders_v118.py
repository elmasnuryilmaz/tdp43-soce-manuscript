"""One-time maintenance of current entrypoints after the scientific audit."""
from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=R/'code/build_tables_v3.py'; s=p.read_text()
a=s.index('# S10 NMD'); b=s.index('# S11 APA',a)
s=s[:a]+'''# S10 NMD: current descriptive contrasts, no independent-contrast tests.
from audit_tables_v118 import nmd as build_nmd_descriptive, s13 as build_s13_corrected
build_nmd_descriptive()
'''+s[b:]
a=s.index('# S13 SOAR exon junction level'); b=s.index('# S14',a)
s=s[:a]+'''# S13 includes both exon flanks with chromosome-prefix normalization.
build_s13_corrected()
'''+s[b:]
p.write_text(s)
p=R/'code/recalc_nmd_shared_control.R'
p.write_text('''#!/usr/bin/env Rscript
# Current NMD entrypoint: descriptive contrasts only. The previous four-contrast
# t/sign tests omitted shared-baseline covariance and must not be regenerated.
args <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args[grepl("^--file=", args)][1])
script <- file.path(dirname(normalizePath(self)), "audit_tables_v118.py")
python <- Sys.getenv("PYTHON", "python3")
status <- system2(python, shQuote(script))
quit(status=status)
''')
p=R/'code/stager_confirmation.py'
p.write_text('''#!/usr/bin/env python3
"""Show nine STIM1 features, seven evaluable confirmations and two missing tests."""
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]
d=pd.read_csv(root/'source_data/STIM1_transcript_test_eligibility.csv')
print(d[['feature_id','pvalue','stageR_confirmation_adjusted_p','stageR_status']].to_string(index=False))
print("Total:",len(d),"Evaluable:",d.stageR_eligible.sum(),"Not evaluated:",(~d.stageR_eligible).sum())
''')
p=R/'code/build_source_data.py'; s=p.read_text()
a=s.index('F, pa = stats.f_oneway'); b=s.index('for g in genes:',s.index('running = 0.0',a))
s=s[:a]+'''rows.append(dict(panel="1A", measurement="TARDBP silencing", group="shTDP-43 vs controls",
                 n=4, mean=round(100 * (1 - kd.mean() / nt.mean()), 1), SEM="",
                 test="descriptive only; four technical RT-qPCR repeats per group", p_value=""))
genes = ["TRPC1", "STIM1", "ORAI1", "ATP2A3"]
'''+s[b:]
s=s.replace('test="two-tailed Student\'s t-test; Holm-adjusted across four targets",\n                     p_value=f"{holm[g]:.6g}"','test="descriptive only; four technical RT-qPCR repeats per group", p_value=""')
s=s.replace('RT-qPCR of TARDBP in three groups (n = 4 each); raw Ct for ',
            'RT-qPCR of TARDBP in three groups (four technical repeats each); raw Ct for ')
s=s.replace('Experiment dates','Recorded assay dates')
s=s.replace('Summary_stats     group means and SEM for Figures 1 and 2; Figure 1 tests',
            'Summary_stats     group means and technical SEM; all laboratory data descriptive')
p.write_text(s)

# Current figure entrypoints invoke the corrected builders; original source versions
# remain recoverable through Git, without regenerating superseded claims.
(R/'code/fig_main_lab.py').write_text('"""Current Figure 1, descriptive technical repeats."""\nfrom fig_lab_descriptive_v118 import *\n')
p=R/'code/fig_supp_rna_processing.py'; s=p.read_text()
s=s.replace('PKG = Path("/Users/elmas/Desktop/MAKALE/11_NEUROCHEMISTRY_INTERNATIONAL_FIGURE_REVISION")','PKG = Path(__file__).resolve().parents[1]')
a=s.index('# S5.'); b=s.index('# S6.',a)
s=s[:a]+'''# S5 uses the corrected descriptive table, without independence claims.
import fig_nmd_descriptive_v118

'''+s[b:]; p.write_text(s)
print('Current table, NMD, stageR and figure entrypoints aligned')

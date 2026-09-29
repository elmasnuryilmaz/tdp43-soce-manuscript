#!/usr/bin/env python3
"""Show nine STIM1 features, seven evaluable confirmations and two missing tests."""
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]
d=pd.read_csv(root/'source_data/STIM1_transcript_test_eligibility.csv')
print(d[['feature_id','pvalue','stageR_confirmation_adjusted_p','stageR_status']].to_string(index=False))
print("Total:",len(d),"Evaluable:",d.stageR_eligible.sum(),"Not evaluated:",(~d.stageR_eligible).sum())

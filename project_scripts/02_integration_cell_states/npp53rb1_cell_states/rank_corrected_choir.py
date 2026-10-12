"""Rank fixed Leiden candidates against an internally audited corrected reference."""
import hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score,adjusted_rand_score

project=Path(__file__).resolve().parents[1]
source=project/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42'
out=source/'choir_label_fix_seed42'
assert (out/'CHOIR_LABEL_FIX_VERIFIED.json').exists()
assert json.loads((out/'CHOIR_LABEL_AUDIT.json').read_text())['valid']
assert not (out/'RANKED.json').exists()
plan=json.loads((out/'CORRECTION_PLAN.json').read_text())
for name,digest in plan['input_sha256'].items():
    assert hashlib.sha256((out/name).read_bytes()).hexdigest()==digest
    assert hashlib.sha256((source/name).read_bytes()).hexdigest()==digest
assign=pd.read_csv(source/'resolution_assignments_PRIVATE.csv',index_col=0,dtype=str)
labels=pd.read_csv(out/'choir_labels_PRIVATE.csv',index_col=0,dtype=str)
assert assign.index.equals(labels.index) and assign.shape==(42806,36)
assert not labels.isna().any().any()
rows=[]
for r in np.round(np.arange(.1,1.21,.1),1):
    groups=[assign[f'leiden_r{r:.1f}_seed{s}'] for s in [42,19,888]]
    amis=[adjusted_mutual_info_score(labels.CHOIR_cluster,x,average_method='arithmetic') for x in groups]
    rows.append(dict(resolution=float(r),clusters=groups[0].nunique(),AMI=amis[0],AMI_seed19=amis[1],AMI_seed888=amis[2]))
rank=pd.DataFrame(rows).sort_values(['AMI','resolution'],ascending=[False,True])
rank['rank']=range(1,13)
rank.to_csv(out/'AMI_ranking.csv',index=False)
rank.head(5).to_csv(out/'scSHC_candidates.csv',index=False)
previous=pd.read_csv(source/'choir_labels_PRIVATE.csv',index_col=0,dtype=str)
assert previous.index.equals(labels.index)
pd.crosstab(previous.CHOIR_cluster,labels.CHOIR_cluster).to_csv(out/'warned_vs_corrected_overlap.csv')
manifest=dict(cells=42806,seed=42,candidates=rank.head(5).resolution.tolist(),
    ranking='Arithmetic AMI; descending score, ascending resolution tie break; corrected reference only',
    label_sha256=hashlib.sha256((out/'choir_labels_PRIVATE.csv').read_bytes()).hexdigest(),
    assignment_sha256=hashlib.sha256((source/'resolution_assignments_PRIVATE.csv').read_bytes()).hexdigest(),
    warned_vs_corrected_ARI=adjusted_rand_score(previous.CHOIR_cluster,labels.CHOIR_cluster),
    previous_reference_rejected=True)
(out/'RANKED.json').write_text(json.dumps(manifest,indent=2))
print(rank.to_string(index=False))
print(json.dumps(manifest,indent=2))

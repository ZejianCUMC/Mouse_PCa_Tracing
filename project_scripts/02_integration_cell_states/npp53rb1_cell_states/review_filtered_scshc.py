"""Post-test composition and cross-candidate agreement; no automatic state labels."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

project=Path(__file__).resolve().parents[1]
out=project/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42'
assert (out/'COMPUTATION_VERIFIED.json').exists()
dest=out/'scshc_review'
dest.mkdir(exist_ok=True)
assert not (dest/'COMPLETE.json').exists()
meta=pd.read_csv(out/'metadata_PRIVATE.csv',index_col=0)
scores=pd.read_csv(out/'expression_review/signature_scores_PRIVATE.csv',index_col=0)
screen=pd.read_csv(out/'screening_labels_PRIVATE.csv',index_col=0,dtype=str)
assert scores.index.equals(meta.index) and screen.index.equals(meta.index)
rank=pd.read_csv(out/'top5_evidence.csv')
partitions={};rows=[];means=[]
for item in rank.itertuples():
    tag=f'r{item.resolution:.1f}'
    labels=pd.read_csv(out/f'scSHC_{tag}_labels_PRIVATE.csv',index_col=0,dtype=str)
    assert labels.index.equals(meta.index)
    group=labels.scSHC_cluster
    partitions[tag]=group
    for cl in sorted(group.unique()):
        mask=group.eq(cl);sub=meta.loc[mask]
        rows.append(dict(rank=int(item.rank),resolution=item.resolution,scSHC_cluster=cl,cells=len(sub),
            input_clusters=labels.loc[mask,'leiden_cluster'].nunique(),samples=sub.sample_id.nunique(),
            runs=sub.run_id.nunique(),dominant_sample=sub.sample_id.value_counts().index[0],
            dominant_sample_fraction=sub.sample_id.value_counts().iloc[0]/len(sub),
            dominant_run_fraction=sub.run_id.value_counts().iloc[0]/len(sub),
            median_genes=sub.n_genes_by_counts.median(),median_UMI=sub.total_counts.median(),
            median_mt=sub.pct_counts_mt.median(),
            reporter_positive=int(((sub.raw_GFP>0)|(sub.raw_mCherry>0)).sum()),
            historical_cl2_fraction=screen.loc[mask,'cluster'].eq('2').mean()))
        means.append(dict(resolution=item.resolution,scSHC_cluster=cl,**scores.loc[mask].mean().to_dict()))
    for field in ['sample_id','run_id','treatment']:
        pd.crosstab(group,meta[field]).to_csv(dest/f'{tag}_{field}.csv')
quality=pd.DataFrame(rows)
assert quality.groupby('resolution').cells.sum().eq(42806).all()
quality.to_csv(dest/'group_QC.csv',index=False)
pd.DataFrame(means).to_csv(dest/'group_signature_means.csv',index=False)
agreement=pd.DataFrame({name:{other:adjusted_rand_score(v,w) for other,w in partitions.items()}
                        for name,v in partitions.items()})
assert np.allclose(agreement,agreement.T) and np.allclose(np.diag(agreement),1)
agreement.to_csv(dest/'cross_candidate_ARI.csv')
(dest/'COMPLETE.json').write_text(json.dumps(dict(cells=42806,candidates=5,
    role='Post-selection composition, QC and agreement review; no labels assigned',
    caveat='Agreement is not biological accuracy. Signature scores are descriptive; sample imbalance and retained low-complexity cells require interpretation.'),indent=2))
print(agreement.to_string())
print(quality.to_string(index=False))

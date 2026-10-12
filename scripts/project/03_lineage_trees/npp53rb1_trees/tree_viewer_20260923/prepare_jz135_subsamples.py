"""User-authorized nested, state-stratified JZ135 sampling; no allele resolution."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy import sparse
from prepare_tree_viewer_inputs import sha
P=Path(__file__).resolve().parents[1]
R=P/'lineage_tracing/tree_viewer_20260923_seed42'
C=P/'lineage_tracing/cassiopeia_eight_20260923_seed42'
S=R/'subsample_cassiopeia'
assert not (S/'PREPARED.json').exists()
meta=pd.read_csv(C/'inputs/JZ135/cell_metadata.csv',index_col=0)
a=pd.read_csv(C/'inputs/JZ135/allele_table.csv.gz',dtype={'cellBC':str,'intBC':str})
counts=meta.putative_state.value_counts().sort_index()
rng=np.random.default_rng(42)
orders={s:rng.permutation(meta.index[meta.putative_state==s]) for s in counts.index}
x=sparse.load_npz(R/'inputs/JZ135/counts.npz')
cells=pd.read_csv(R/'inputs/JZ135/cells.csv',index_col=0).index
emb=pd.read_csv(R/'inputs/JZ135/umap.csv',index_col=0)
manifest=json.loads((R/'PREPARED.json').read_text())
manifest['samples'].remove('JZ135');manifest['inputs']=[i for i in manifest['inputs'] if i['sample']!='JZ135']
records=[];previous=set()
for n in [1000,2000,3000]:
    sample=f'JZ135_n{n}'
    quota=(n-len(counts))*counts/counts.sum()
    allocation=np.floor(quota).astype(int)+1
    for state in (quota-np.floor(quota)).sort_values(ascending=False,kind='stable').index[:n-allocation.sum()]:allocation[state]+=1
    chosen=set(c for s in counts.index for c in orders[s][:allocation[s]])
    assert len(chosen)==n and previous<=chosen;previous=chosen
    sub=meta.loc[meta.index.isin(chosen)].copy();suba=a[a.cellBC.isin(chosen)].copy()
    out=S/'inputs'/sample;out.mkdir(parents=True,exist_ok=False)
    sub.to_csv(out/'cell_metadata.csv',index_label='cell')
    suba.to_csv(out/'allele_table.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    qc=dict(sample_id='JZ135',analysis_id=sample,matched_cells=n,allele_rows=len(suba),
        allele_table_sha256=sha(out/'allele_table.csv.gz'),metadata_sha256=sha(out/'cell_metadata.csv'),
        ambiguity_policy='Native tuples, unchanged; no dominant-allele resolution',
        sampling='Proportional by current state; at least one per observed state; nested selections; seed42',
        seed=42,state_counts=allocation.to_dict(),full_state_counts=counts.to_dict())
    (out/'INPUT_QC.json').write_text(json.dumps(qc,indent=2))
    dest=R/'inputs'/sample;dest.mkdir(parents=True,exist_ok=False)
    ix=cells.get_indexer(sub.index);assert (ix>=0).all()
    sparse.save_npz(dest/'counts.npz',x[ix])
    pd.DataFrame(index=sub.index).to_csv(dest/'cells.csv',index_label='cell')
    (dest/'genes.csv').write_bytes((R/'inputs/JZ135/genes.csv').read_bytes())
    emb.loc[sub.index].to_csv(dest/'umap.csv')
    manifest['samples'].append(sample)
    manifest['inputs'].append(dict(sample=sample,cells=n,genes=x.shape[1],files={f.name:sha(f) for f in dest.iterdir()},no_state_labels=True))
    records.append(qc)
    print(sample,n,allocation.to_dict(),flush=True)
(S/'PREPARED.json').write_text(json.dumps(dict(samples=[r['analysis_id'] for r in records],seed=42,subsamples=records),indent=2))
manifest['full_JZ135_hold']='User requested stopping full-cell run; only n1000/n2000/n3000 fit. Full input retained for audit only.'
(R/'PREPARED.json').write_text(json.dumps(manifest,indent=2))
(R/'JZ135_SUBSAMPLING.json').write_text(json.dumps(records,indent=2))

"""Freeze expression/coordinates, deliberately excluding state labels from fitting inputs."""
import json,hashlib
from pathlib import Path
import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse

P=Path(__file__).resolve().parents[3]   # scripts -> run dir -> lineage_tracing -> project
RUN=P/'lineage_tracing/tree_viewer_mjz026_028_20260925_seed42'
CAS=P/'lineage_tracing/cassiopeia_mjz026_028_20260925_seed42'
FINAL=P/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/exclude_current_c17_c19_seed42'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

if __name__=='__main__':
    assert not (RUN/'PREPARED.json').exists()
    RUN.mkdir(parents=True,exist_ok=True)
    source=P/'analysis/min400_20260921/merged_counts.h5ad'
    a=ad.read_h5ad(source,backed='r')
    umap=pd.read_csv(FINAL/'umap_PRIVATE.csv',index_col=0)
    samples=json.loads((CAS/'PREPARED.json').read_text())['samples']
    records=[]
    for sample in samples:
        meta=pd.read_csv(CAS/'inputs'/sample/'cell_metadata.csv',index_col=0)
        ix=a.obs_names.get_indexer(meta.index);assert (ix>=0).all()
        # Backed AnnData layers can be materialized by the reader; sparse throughout.
        x=sparse.csr_matrix(a.layers['counts'][ix,:])
        assert x.shape==(len(meta),33702) and np.isfinite(x.data).all() and (x.data>=0).all()
        out=RUN/'inputs'/sample;out.mkdir(parents=True,exist_ok=False)
        sparse.save_npz(out/'counts.npz',x)
        pd.DataFrame(index=meta.index).to_csv(out/'cells.csv',index_label='cell')
        pd.DataFrame(index=a.var_names).to_csv(out/'genes.csv',index_label='gene_id')
        umap.loc[meta.index].to_csv(out/'umap.csv')
        rec=dict(sample=sample,cells=len(meta),genes=x.shape[1],min_detected_genes=int((x>0).sum(axis=1).min()),
                 files={f.name:sha(f) for f in out.iterdir()},no_state_labels=True)
        records.append(rec);print(rec['sample'],rec['cells'],rec['min_detected_genes'],flush=True)
    a.file.close()
    (RUN/'PREPARED.json').write_text(json.dumps(dict(samples=samples,inputs=records,seed=42,
        source_sha256=sha(source),umap_sha256=sha(FINAL/'umap_PRIVATE.csv'),
        scope='14 MJZ026-028 Np53Rb1 mice (hold lifted 2026-09-23); states only for display'),indent=2))

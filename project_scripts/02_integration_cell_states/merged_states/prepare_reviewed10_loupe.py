"""Reuse verified, identical-cohort counts/projections; replace annotations only."""
from pathlib import Path
import json, hashlib, shutil, sys
import anndata as ad
import numpy as np
import pandas as pd
p=Path(__file__).resolve().parents[1]
m=p/'analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42/approved_merges'
d=m/'reviewed10'; old=m/'loupe_export'; out=d/'loupe_export'
if '--oxidative' in sys.argv:out=d/'loupe_export_oxidative'
assert not out.exists() or not any(out.iterdir())
a=ad.read_h5ad(d/'named_graph_PRIVATE.h5ad');b=ad.read_h5ad(m/'named_graph_PRIVATE.h5ad')
assert a.obs_names.equals(b.obs_names)
for key in ['X_scVI','X_umap']:np.testing.assert_array_equal(a.obsm[key],b.obsm[key])
assert (a.obsp['connectivities']!=b.obsp['connectivities']).nnz==0
assert json.loads((old/'FINAL_EXPORT_VERIFIED.json').read_text())['status']=='PASS'
meta=pd.read_csv(old/'metadata_PRIVATE.csv',index_col=0)
assert meta.index.equals(a.obs_names)
for col in ['putative_state','final_cluster']:meta[col]=a.obs[col].astype(str)
assert meta.loc[meta.final_cluster.eq('c3'),'putative_state'].unique().tolist()==['c3_Lum-basal-hybrid']
out.mkdir(exist_ok=True)
for name in ['all_counts_PRIVATE.h5','features.csv','source_feature_mapping.csv','COUNTS_COMPLETE.json','umap_PRIVATE.csv','tsne_PRIVATE.csv','diffmap_PRIVATE.csv','diffmap_all_components_PRIVATE.csv']:
    shutil.copy2(old/name,out/name)
    assert hashlib.file_digest((old/name).open('rb'),'sha256').hexdigest()==hashlib.file_digest((out/name).open('rb'),'sha256').hexdigest()
meta.to_csv(out/'metadata_PRIVATE.csv')
r=json.loads((old/'PROJECTIONS_COMPLETE.json').read_text())
r.update(source_sha256=hashlib.file_digest((d/'named_graph_PRIVATE.h5ad').open('rb'),'sha256').hexdigest(),reused_identical_cohort_latent_graph_projections=True)
(out/'PROJECTIONS_COMPLETE.json').write_text(json.dumps(r,indent=2))
print('PASS: same cells/counts/latent/graph/projections; current 10-state labels prepared')

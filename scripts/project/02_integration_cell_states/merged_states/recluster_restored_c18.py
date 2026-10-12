"""Recluster all restored-cohort cells on its verified graph, without geometry changes."""
from pathlib import Path
import json,time,hashlib
import numpy as np,pandas as pd,scanpy as sc
from sklearn.metrics import adjusted_mutual_info_score,adjusted_rand_score
p=Path(__file__).resolve().parents[1];o=p/'analysis/merged_20260923_seed42';src=o/'restored_old_c18_seed42';d=src/'recluster_r1_seed42';d.mkdir(exist_ok=True)
assert not (d/'CLUSTERING_VERIFIED.json').exists()
assert json.loads((src/'VERIFIED.json').read_text())['status']=='PASS'
a=sc.read_h5ad(src/'restored_graph_PRIVATE.h5ad');assert a.n_obs==55133
start=time.time();a.obs['previous_group']=a.obs.final_cluster.astype(str);a.obs['previous_label']=a.obs.putative_state.astype(str)
for seed in [42,19,888]:
 sc.tl.leiden(a,resolution=1.0,random_state=seed,key_added=f'restored_leiden_seed{seed}')
 print('seed',seed,'groups',a.obs[f'restored_leiden_seed{seed}'].nunique(),flush=True)
col='restored_leiden_seed42';a.obs['raw_cluster']=a.obs[col].astype(str)
ct=pd.crosstab(a.obs.raw_cluster,a.obs.previous_group);ct.to_csv(d/'raw_previous_crosswalk.csv');a.obs.to_csv(d/'raw_metadata_PRIVATE.csv');a.write_h5ad(d/'clustered_graph_PRIVATE.h5ad',compression='gzip')
metrics={str(s):{'groups':int(a.obs[f'restored_leiden_seed{s}'].nunique()),'ARI_to_42':adjusted_rand_score(a.obs[col],a.obs[f'restored_leiden_seed{s}']),'AMI_to_42':adjusted_mutual_info_score(a.obs[col],a.obs[f'restored_leiden_seed{s}'])} for s in [42,19,888]}
b=sc.read_h5ad(src/'restored_graph_PRIVATE.h5ad');assert a.obs_names.equals(b.obs_names)
for k in b.obsm:np.testing.assert_array_equal(a.obsm[k],b.obsm[k])
for k in b.obsp:assert (a.obsp[k]!=b.obsp[k]).nnz==0
r=dict(status='PASS',cells=a.n_obs,resolution=1.0,seeds=metrics,latent_graph_umap_unchanged=True,new_leiden_assignments=True,scSHC_rerun=False,CHOIR_rerun=False,selection='Reuse user-selected resolution1.0; original-cohort CHOIR/AMI/scSHC evidence is historical, not validation of restored-cohort clustering',elapsed_seconds=time.time()-start)
(d/'CLUSTERING_VERIFIED.json').write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True);print(ct.to_string(),flush=True)

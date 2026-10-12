"""Fresh tumour-only HVGs; keep the user-defined cohort fixed. Stop downstream at UMAP."""
import os, json
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.io import mmwrite

run=Path(os.environ['RB1_RUN']); out=run/'tumor_umap_20260921';out.mkdir(exist_ok=True)
assert not (out/'tumor_hvg_counts.h5ad').exists(), 'Refuse overwriting prepared input'
flags=pd.read_csv(run/'harmony_major_types/cluster_tumor_cell_classification.csv',index_col=0)
selected=flags.index[flags.tumor_by_cluster.astype(str).str.lower()=='true']
assert len(selected)==43057 and not flags.loc[selected,'harmony_cluster'].isin([4,5,6,7]).any()
source=sc.read_h5ad(run/'major_types_full.h5ad',backed='r')
b=source[selected,:].to_memory();source.file.close()
obs=b.obs.copy();obs['run_id']=obs.library.astype(str)
obs['parent_harmony_cluster']=flags.loc[b.obs_names,'harmony_cluster'].astype(str)
for col in list(obs):
    if col.startswith('clusters_') or col in ['primary_cluster','major_type','putative_tumor']:
        del obs[col]
obs['treatment']=obs.treatment.astype(str).replace({'Castration+Enzalutamide':'Castration+ENZA'})
assert obs.sample_id.nunique()==22 and obs.run_id.nunique()==6
assert obs.groupby('sample_id',observed=True).run_id.nunique().max()==1
a=ad.AnnData(X=b.layers['counts'].copy(),obs=obs,var=b.var.copy())
a.obs['total_counts']=np.asarray(a.X.sum(axis=1)).ravel()
for col in list(a.var):
    if col.startswith('highly_variable') or col in ['means','dispersions','dispersions_norm','variances','variances_norm']:
        del a.var[col]
del b
assert np.isfinite(a.X.data).all() and (a.X.data>=0).all()
assert np.allclose(a.X.data,np.round(a.X.data)), 'scVI requires count inputs'
assert int(((a.obs.raw_GFP>0)|(a.obs.raw_mCherry>0)).sum())==39683
sc.pp.filter_genes(a,min_cells=3)
eligible=~a.var_names.str.lower().isin(['gfp','mcherry','egfp','tdtomato','cre','luciferase'])
h=a[:,eligible].copy()
sc.pp.highly_variable_genes(h,n_top_genes=3000,flavor='seurat_v3',batch_key='run_id')
genes=h.var_names[h.var.highly_variable].tolist()
assert len(genes)==3000 and not {'GFP','mCherry'}&set(genes)
a.var['highly_variable']=a.var_names.isin(genes)
a.write_h5ad(out/'tumor_full_counts.h5ad')
h=a[:,genes].copy();h.layers['counts']=h.X.copy();h.write_h5ad(out/'tumor_hvg_counts.h5ad')
mmwrite(out/'tumor_hvg_counts.mtx',h.X.T)
(out/'hvg.txt').write_text('\n'.join(genes)+'\n')
a.obs.to_csv(out/'metadata.csv')
pd.crosstab(a.obs.sample_id,a.obs.run_id).to_csv(out/'sample_run_design.csv')
pd.crosstab(a.obs.run_id,a.obs.treatment).to_csv(out/'run_treatment_design.csv')
# Raw-expression programs are diagnostics only, never correction covariates.
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
programs={'Basal':['Krt5','Krt14','Trp63','Krt17'],
          'Luminal':['Krt8','Krt18','Krt19','Epcam'],
          'Mesenchymal':['Vim','Zeb1','Col1a1','Col3a1'],
          'Cycling':['Mki67','Top2a']}
scores=pd.DataFrame(index=a.obs_names)
for name,markers in programs.items():
    scores[name]=np.asarray(a[:,[g for g in markers if g in a.var_names]].X.mean(axis=1)).ravel()
scores.to_csv(out/'expression_programs.csv')
manifest={'cells':a.n_obs,'samples':22,'runs':6,'hvgs':len(genes),'HVG_method':'seurat_v3, batch_key=run_id',
 'definition':'All cells outside full-cohort Harmony theta1 Leiden0.3 clusters4/5/6/7; reporter-independent; cluster9 retained',
 'tumor_reporter_positive':39683,'tumor_reporter_negative':3374,'stop':'UMAP only; no new clustering/annotation/DE',
 'source':str(run/'major_types_full.h5ad'),'scanpy_version':sc.__version__,'seed':19}
(out/'prepared.json').write_text(json.dumps(manifest,indent=2));print(manifest,flush=True)

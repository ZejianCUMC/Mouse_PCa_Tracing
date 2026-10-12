"""Preserve current latent/graph/UMAP; compute additional Loupe projections."""
import os,time,json,hashlib,sys
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/npp_merged_numba')
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
from scipy.sparse.csgraph import connected_components
from sklearn.manifold import TSNE
p=Path(__file__).resolve().parents[1]
c=p/'analysis/merged_20260923_seed42/filtered_no_c17_c18_seed42/resolution_1p0/final_cleanup'
merged='--merged' in sys.argv
if merged:c=p/'analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42/approved_merges'
d=c/'loupe_export';d.mkdir(exist_ok=True)
assert not (d/'PROJECTIONS_COMPLETE.json').exists()
source=c/('named_graph_PRIVATE.h5ad' if merged else 'final_named_graph_PRIVATE.h5ad');digest=hashlib.file_digest(source.open('rb'),'sha256').hexdigest()
start=time.time();a=sc.read_h5ad(source);n=55133 if merged else 50241;assert a.shape==(n,0)
z=a.obsm['X_scVI'];assert z.shape==(n,30)
ncomp,_=connected_components(a.obsp['connectivities']);assert ncomp==1,'Disconnected graph requires review'
pd.DataFrame(a.obsm['X_umap'],index=a.obs_names,columns=['UMAP_1','UMAP_2']).to_csv(d/'umap_PRIVATE.csv')
print('Computing diffusion map on unchanged saved scVI graph',flush=True)
sc.tl.diffmap(a,n_comps=15,random_state=42)
ev=a.uns['diffmap_evals'];dc=a.obsm['X_diffmap']
assert np.isclose(ev[0],1,atol=1e-5) and ev[1]<.999999 and np.isfinite(dc).all()
pd.DataFrame(dc[:,1:3],index=a.obs_names,columns=['DC_1','DC_2']).to_csv(d/'diffmap_PRIVATE.csv')
pd.DataFrame(dc,index=a.obs_names,columns=['DC_stationary']+[f'DC_{i}' for i in range(1,15)]).to_csv(d/'diffmap_all_components_PRIVATE.csv')
print('Computing tSNE from all 30 current scVI dimensions',flush=True)
model=TSNE(n_components=2,perplexity=30,init='pca',learning_rate='auto',max_iter=1000,random_state=42,n_jobs=4,verbose=2)
ts=model.fit_transform(z);assert ts.shape==(n,2) and np.isfinite(ts).all()
pd.DataFrame(ts,index=a.obs_names,columns=['tSNE_1','tSNE_2']).to_csv(d/'tsne_PRIVATE.csv')
assert digest==hashlib.file_digest(source.open('rb'),'sha256').hexdigest()
record=dict(status='PASS',cells=n,seed=42,latent_dimensions=30,UMAP='exact saved current coordinates',diffmap_graph='unchanged saved current 15-neighbor scVI graph',connected_components=int(ncomp),diffmap_eigenvalues=ev.tolist(),stationary_component_excluded=True,tsne_perplexity=30,tsne_iterations=1000,tsne_init='pca',tsne_learning_rate='auto',tsne_kl_divergence=float(model.kl_divergence_),source_sha256=digest,seconds=time.time()-start)
(d/'PROJECTIONS_COMPLETE.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2),flush=True)

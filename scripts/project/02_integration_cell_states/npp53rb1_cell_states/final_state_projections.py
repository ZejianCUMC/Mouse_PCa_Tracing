"""Fresh tSNE and diffusion map for the final 40,141-cell subset."""
import os,tempfile
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.environ.setdefault('NUMBA_CACHE_DIR',tempfile.mkdtemp(prefix='rb1-final-numba-'))
from pathlib import Path
import json,time,hashlib
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from sklearn.manifold import TSNE
p=Path(__file__).resolve().parents[1]
out=p/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/exclude_current_c17_c19_seed42'
start=time.time()
assert not (out/'FINAL_PROJECTIONS_COMPLETE.json').exists()
latent=pd.read_csv(out/'scvi_latent_PRIVATE.csv',index_col=0)
assert latent.shape==(40141,30) and latent.index.is_unique
protected=['scvi_latent_PRIVATE.csv','state_labels_PRIVATE.csv','umap_PRIVATE.csv','connectivities_PRIVATE.npz','distances_PRIVATE.npz']
hashes={f:hashlib.sha256((out/f).read_bytes()).hexdigest() for f in protected}
print('Diffusion map from existing final-cohort graph',flush=True)
a=ad.AnnData(np.zeros((len(latent),1),dtype=np.float32),obs=pd.DataFrame(index=latent.index))
a.obsm['X_scVI']=latent.to_numpy(dtype=np.float32)
a.obsp['connectivities']=sparse.load_npz(out/'connectivities_PRIVATE.npz')
a.obsp['distances']=sparse.load_npz(out/'distances_PRIVATE.npz')
a.uns['neighbors']={'connectivities_key':'connectivities','distances_key':'distances','params':{'n_neighbors':15,'method':'umap','random_state':42,'metric':'euclidean','use_rep':'X_scVI'}}
ncomp,_=connected_components(a.obsp['connectivities'])
assert ncomp==1, 'Disconnected graph needs explicit diffusion-component review'
sc.tl.diffmap(a,n_comps=15,random_state=42)
ev=a.uns['diffmap_evals'];dc=a.obsm['X_diffmap']
assert np.isclose(ev[0],1,atol=1e-5) and ev[1]<.999999 and np.isfinite(dc).all()
# Skip stationary first eigenvector; display two first nonstationary components.
pd.DataFrame(dc[:,1:3],index=latent.index,columns=['DC_1','DC_2']).to_csv(out/'diffmap_PRIVATE.csv')
pd.DataFrame(dc,index=latent.index,columns=['DC_stationary']+[f'DC_{i}' for i in range(1,15)]).to_csv(out/'diffmap_all_components_PRIVATE.csv')
print('Diffusion map complete; computing tSNE from all30 retained scVI dimensions',flush=True)
ts=TSNE(n_components=2,perplexity=30,init='pca',learning_rate='auto',max_iter=1000,random_state=42,n_jobs=4,verbose=1).fit_transform(latent.to_numpy())
assert ts.shape==(40141,2) and np.isfinite(ts).all()
pd.DataFrame(ts,index=latent.index,columns=['tSNE_1','tSNE_2']).to_csv(out/'tsne_PRIVATE.csv')
assert hashes=={f:hashlib.sha256((out/f).read_bytes()).hexdigest() for f in protected}
(out/'FINAL_PROJECTIONS_COMPLETE.json').write_text(json.dumps(dict(cells=40141,seed=42,tsne_perplexity=30,tsne_iterations=1000,tsne_init='pca',tsne_learning_rate='auto',diffmap_graph='saved final40141cell scVI15-neighbor graph',diffmap_components=15,diffmap_display_columns=[1,2],stationary_component_excluded=True,graph_connected_components=int(ncomp),diffmap_eigenvalues=ev.tolist(),source_sha256=hashes,elapsed_seconds=time.time()-start),indent=2))
print('PASS: new final-cohort tSNE/diffmap; fixed UMAP/labels/latent protected.',flush=True)

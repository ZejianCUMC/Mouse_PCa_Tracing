"""Restore exactly old-c18; rebuild graph/UMAP, preserving latent and labels."""
from pathlib import Path
import json,time,os,hashlib,textwrap
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
p=Path(__file__).resolve().parents[1];o=p/'analysis/merged_20260923_seed42';c=o/'filtered_no_c17_c18_seed42/resolution_1p0/final_cleanup';d=o/'restored_old_c18_seed42';d.mkdir(exist_ok=True)
assert not (d/'VERIFIED.json').exists()
start=time.time();current=sc.read_h5ad(c/'final_named_graph_PRIVATE.h5ad');old=pd.read_csv(o/'states_seed42/final/final_metadata_PRIVATE.csv',index_col=0);restore=old[old.putative_state.str.startswith('c18 ')].copy()
assert len(current)==50241 and len(restore)==4892 and not set(current.obs_names)&set(restore.index)
restore['final_cluster']='old-c18';restore['putative_state']='old-c18 Restored mesenchymal (provisional)';restore['restoration_status']='Restored old-c18';meta=current.obs.copy();meta['restoration_status']='Previously retained';meta=pd.concat([meta,restore],axis=0)
assert len(meta)==55133 and meta.index.is_unique
fit=json.loads((o/'INTEGRATION_SELECTION.json').read_text())['selected_fit'];latent=pd.read_csv(o/fit/'latent_PRIVATE.csv',index_col=0);z=latent.loc[meta.index].to_numpy()
np.testing.assert_array_equal(z[:50241],current.obsm['X_scVI'])
assert not set(meta.index)&set(old.index[old.putative_state.str.startswith('c17 ')])
excluded=pd.read_csv(c/'excluded_cells_PRIVATE.csv',index_col=0);assert not set(meta.index)&set(excluded.index)
b=ad.AnnData(obs=meta);b.obsm['X_scVI']=z
(d/'STARTED.json').write_text(json.dumps(dict(pid=os.getpid(),seed=42,neighbors=15,min_dist=.3,cells=len(meta),restored_cells=4892,estimated_minutes=15,estimated_ram_gb=8),indent=2))
print('Building graph',b.n_obs,flush=True);sc.pp.neighbors(b,use_rep='X_scVI',n_neighbors=15,random_state=42,metric='euclidean');print('Rebuilding UMAP',flush=True);sc.tl.umap(b,random_state=42,min_dist=.3)
assert not np.allclose(b.obsm['X_umap'][:50241],current.obsm['X_umap'])
np.testing.assert_array_equal(b.obsm['X_scVI'],z);assert np.isfinite(b.obsm['X_umap']).all()
names=list(pd.read_csv(c/'state_order.csv').state)+['old-c18 Restored mesenchymal (provisional)'];pal={s:matplotlib.colors.to_hex(plt.get_cmap('tab20')(i)) for i,s in enumerate(names)};pal[names[-1]]='#E7298A'
tr={'Intact':'#4c78a8','Castration':'#f58518','Castration+ENZA':'#54a24b'};rp={k:matplotlib.colors.to_hex(plt.get_cmap('tab20')(i)) for i,k in enumerate(sorted(meta.run_id.unique()))};xy=b.obsm['X_umap'];draw=np.random.default_rng(42).permutation(len(meta))
fig,axes=plt.subplots(1,3,figsize=(23,10))
for ax,col,palette,title in zip(axes,['putative_state','treatment','run_id'],[pal,tr,rp],['Inherited states + restored old-c18','Treatment','Library / run']):
 ax.scatter(xy[draw,0],xy[draw,1],c=[palette[v] for v in meta[col].iloc[draw]],s=3,alpha=.85,linewidths=0,rasterized=True);ax.set_title(title,fontsize=16);ax.set_xlabel('UMAP_1',fontsize=13);ax.set_ylabel('UMAP_2',fontsize=13);ax.set_xticks([]);ax.set_yticks([])
 if col=='putative_state':
  for k in palette:ax.text(*np.median(xy[meta[col].eq(k)],axis=0),k.split()[0],fontsize=9,ha='center',weight='bold',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
 ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color=v,markersize=5,label=textwrap.fill(k,36)) for k,v in palette.items()],loc='upper center',bbox_to_anchor=(.5,-.12),ncol=2 if col!='treatment' else 1,frameon=False,fontsize=8 if col=='putative_state' else 10)
fig.subplots_adjust(left=.04,right=.99,top=.93,bottom=.40,wspace=.22);fig.savefig(d/'restored_joint_umap.png',dpi=170,bbox_inches='tight');plt.close(fig)
b.write_h5ad(d/'restored_graph_PRIVATE.h5ad',compression='gzip');meta.to_csv(d/'metadata_PRIVATE.csv');pd.DataFrame(xy,index=meta.index,columns=['UMAP_1','UMAP_2']).to_csv(d/'umap_PRIVATE.csv')
meta.groupby(['final_cluster','putative_state','model','sample_id','run_id','treatment'],observed=True).size().rename('n_cells').reset_index().to_csv(d/'state_sample_counts.csv',index=False)
meta.groupby(['final_cluster','putative_state'],observed=True).size().rename('cells').reset_index().to_csv(d/'state_counts.csv',index=False)
# Independent output read-back checks.
v=sc.read_h5ad(d/'restored_graph_PRIVATE.h5ad');assert v.obs_names.equals(meta.index);np.testing.assert_array_equal(v.obsm['X_scVI'],z);pd.testing.assert_series_equal(v.obs.final_cluster.astype(str),meta.final_cluster.astype(str));assert (v.obsp['connectivities']-b.obsp['connectivities']).nnz==0
pd.testing.assert_series_equal(v.obs.loc[current.obs_names,'putative_state'].astype(str),current.obs.putative_state.astype(str));assert set(v.obs_names)-set(current.obs_names)==set(restore.index)
record=dict(status='PASS',cells=55133,retained_prior=50241,restored=4892,groups=13,seed=42,neighbors=15,min_dist=.3,latent_unchanged=True,graph_umap_rebuilt=True,reclustered=False,old_c17_still_excluded=True,other_227_exclusions_preserved=True,elapsed_seconds=time.time()-start,files={k:hashlib.sha256((d/k).read_bytes()).hexdigest() for k in ['restored_graph_PRIVATE.h5ad','restored_joint_umap.png','metadata_PRIVATE.csv','state_counts.csv']})
(d/'VERIFIED.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2),flush=True)

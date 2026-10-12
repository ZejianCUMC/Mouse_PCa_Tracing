"""Descriptive integration diagnostics; no automatic winner or biological labels."""
import sys,json,time
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
from scipy.stats import spearmanr
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]);variant=sys.argv[2] if len(sys.argv)>2 else '';assert variant in ['', 'endogenous']
suffix='_endogenous' if variant else ''
modes=[sys.argv[3]] if len(sys.argv)>3 else ['run','run_sample','sample']
dest=out/f'integration_review{suffix}';dest.mkdir(exist_ok=True)
assert not (dest/'COMPLETE.json').exists()
programs={
 'Luminal':['Krt8','Krt18','Epcam','Krt19','Foxa1','Nkx3-1','Ar','Pbsn'],
 'Basal_squamous':['Krt5','Krt14','Krt17','Trp63','Krt6a','Krt6b'],
 'Mesenchymal':['Vim','Col1a1','Col1a2','Col3a1','Fn1','Pdgfra','Zeb1','Sparc'],
 'Neuroendocrine':['Ascl1','Insm1','Chga','Chgb','Syp','Neurod1','Dll3'],
 'Cycling':['Mki67','Top2a','Cdk1','Ccnb1','Ube2c','Birc5'],
 'Myeloid':['Ptprc','Lyz2','Csf1r','Tyrobp','Adgre1','Lst1'],
 'Lymphoid':['Ptprc','Cd3d','Cd3e','Trbc2','Cd79a','Ms4a1'],
 'Endothelial':['Pecam1','Kdr','Cdh5','Emcn','Vwf'],
 'Ionocyte':['Foxi1','Cftr','Ascl3','Atp6v1b1']}
full=sc.read_h5ad(out/'shared_counts.h5ad');obs=full.obs.copy()
symbols=pd.Index(full.var.gene_symbol);allmarkers=list(dict.fromkeys(g for gs in programs.values() for g in gs if g in symbols))
idx=symbols.get_indexer(allmarkers)
counts=full.X[:,idx].toarray();norm=np.log1p(counts/obs.total_counts.values[:,None]*10000)
expr=pd.DataFrame(norm,index=obs.index,columns=allmarkers)
scores=pd.DataFrame({k:expr[[g for g in gs if g in expr]].mean(axis=1) for k,gs in programs.items()})
# Transparent marker expression summaries, not reference-ground-truth labels.
lineage=['Luminal','Basal_squamous','Mesenchymal','Neuroendocrine','Myeloid','Lymphoid','Endothelial','Ionocyte']
obs['marker_anchor']=scores[lineage].idxmax(axis=1)
scores.to_csv(dest/'marker_scores_PRIVATE.csv')
pd.DataFrame({'gene':allmarkers,'pct_positive':(counts>0).mean(0)*100,'mean_log1p':norm.mean(0)}).to_csv(dest/'canonical_markers_global.csv',index=False)
del full,counts,norm
records=[]
for mode in modes:
 path=out/f'scvi_{mode}_seed42{suffix}';assert (path/'COMPLETE.json').exists()
 z=pd.read_csv(path/'latent_PRIVATE.csv',index_col=0);assert z.index.equals(obs.index)
 a=ad.AnnData(X=np.zeros((len(obs),1),dtype=np.float32),obs=obs.copy());a.obsm['X_scVI']=z.to_numpy(dtype=np.float32)
 print('Neighbors and UMAP',mode,flush=True)
 sc.pp.neighbors(a,use_rep='X_scVI',n_neighbors=15,random_state=42)
 sc.tl.umap(a,random_state=42,min_dist=.3)
 pd.DataFrame(a.obsm['X_umap'],index=obs.index,columns=['UMAP_1','UMAP_2']).to_csv(path/'umap_PRIVATE.csv')
 a.write_h5ad(path/'graph_umap.h5ad',compression='gzip')
 graph=a.obsp['distances'].tocsr();rows=[]
 for i in range(len(obs)):
  nbr=graph.indices[graph.indptr[i]:graph.indptr[i+1]]
  same=obs.marker_anchor.values[nbr]==obs.marker_anchor.values[i]
  comparable=nbr[same]
  if len(comparable)>=5:
   freq=pd.Series(obs.run_id.values[comparable]).value_counts(normalize=True).values
   freq=freq[freq>0]
   entropy=float(-(freq*np.log(freq)).sum())
   crossrun=float(np.mean(obs.run_id.values[comparable]!=obs.run_id.values[i]))
  else:entropy=np.nan;crossrun=np.nan
  rows.append((float(same.mean()),entropy,crossrun))
 d=pd.DataFrame(rows,index=obs.index,columns=['marker_anchor_neighbor_agreement','within_anchor_run_entropy','within_anchor_other_run_fraction'])
 d['sample_id']=obs.sample_id;d['model']=obs.model;d['marker_anchor']=obs.marker_anchor
 d.groupby(['model','sample_id']).mean(numeric_only=True).to_csv(path/'sample_diagnostics.csv')
 samplemean=d.groupby('sample_id').mean(numeric_only=True).mean().to_dict()
 for key in scores:
  neighbor_mean=np.array([scores[key].values[graph.indices[graph.indptr[i]:graph.indptr[i+1]]].mean() for i in range(len(obs))])
  samplemean['program_neighbor_spearman_'+key]=float(spearmanr(scores[key],neighbor_mean).statistic)
 records.append({'method':mode,**samplemean})
 for columns,name in [(['model','treatment','run_id'],'design'),(['sample_id','original_state','marker_anchor'],'identity')]:
  fig,axes=plt.subplots(1,3,figsize=(24,6))
  for ax,col in zip(axes,columns):
   cats=sorted(obs[col].astype(str).unique());palette=plt.get_cmap('tab20' if len(cats)<=20 else 'gist_ncar',len(cats))
   for j,cat in enumerate(cats):
    mask=obs[col].astype(str).values==cat;xy=a.obsm['X_umap'][mask]
    ax.scatter(xy[:,0],xy[:,1],s=2,alpha=.7,c=[palette(j)],linewidths=0,rasterized=True,label=cat)
   ax.set_xlabel('UMAP_1',fontsize=12);ax.set_ylabel('UMAP_2',fontsize=12);ax.set_xticks([]);ax.set_yticks([]);ax.set_title(col,fontsize=14)
   ax.legend(loc='upper left',bbox_to_anchor=(1,1),fontsize=7,markerscale=3,frameon=False,ncol=1 if len(cats)<=24 else 2)
  fig.suptitle('scVI NB '+mode+' — seed 42',fontsize=16);fig.tight_layout()
  fig.savefig(dest/f'{mode}_{name}.png',dpi=160,bbox_inches='tight');plt.close(fig)
 print('Diagnostics complete',mode,flush=True)
pd.DataFrame(records).to_csv(dest/'integration_diagnostics.csv',index=False)
(dest/'COMPLETE.json').write_text(json.dumps({'cells':len(obs),'methods':[r['method'] for r in records],'seed':42,'note':'Descriptive diagnostics only. Anchors are marker-score argmax, not ground truth; unbalanced/confounded runs constrain interpretability. No automatic winner.'},indent=2))

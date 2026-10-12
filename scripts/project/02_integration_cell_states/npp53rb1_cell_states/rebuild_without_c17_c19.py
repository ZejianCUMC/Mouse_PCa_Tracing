"""Exclude current c17-c19; new graph/UMAP, unchanged scVI model and c1-c16 labels."""
import os,tempfile
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.environ.setdefault('NUMBA_CACHE_DIR',tempfile.mkdtemp(prefix='rb1-retained-numba-'))
os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='rb1-retained-mpl-'))
from pathlib import Path
import json,time,hashlib,shutil,sys,textwrap
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

project=Path(__file__).resolve().parents[1]
source=project/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42'
out=source/'exclude_current_c17_c19_seed42'
out.mkdir(exist_ok=True)
start=time.time()
protected=['scvi_latent_PRIVATE.csv','umap_PRIVATE.csv','resolution_assignments_PRIVATE.csv','cluster_state_inference/cell_state_labels_PRIVATE.csv']
hashes={f:hashlib.sha256((source/f).read_bytes()).hexdigest() for f in protected}
meta=pd.read_csv(source/'metadata_PRIVATE.csv',index_col=0)
lab=pd.read_csv(source/'cluster_state_inference/cell_state_labels_PRIVATE.csv',index_col=0)
latent=pd.read_csv(source/'scvi_latent_PRIVATE.csv',index_col=0)
assert meta.index.equals(lab.index) and latent.index.equals(lab.index)
excluded=['c17_Cycling-Basal','c18_Cycling-Stress','c19_Mes-lowComplexity']
keep=~lab.inferred_state.isin(excluded)
assert len(lab)==42806 and keep.sum()==40141
assert lab.loc[~keep,'original_cluster'].value_counts().to_dict()=={1:1978,6:651,16:36}
snapshot=out/'previous_report.html'
if not snapshot.exists():shutil.copy2(project/'Npp53Rb1_tracing_report_part1_preprocessing_tumor_states.html',snapshot)
audit=lab.copy();audit['retained']=keep;audit.to_csv(out/'retention_audit_PRIVATE.csv')
meta=meta.loc[keep];lab=lab.loc[keep];latent=latent.loc[keep]
assert lab.inferred_state.nunique()==16 and latent.shape==(40141,30)
manifest=dict(cells=40141,source_cells=42806,excluded_cells=2665,excluded_states=excluded,
    excluded_original_r1_1_clusters=[1,16,6],retained_states=16,primary_seed=42,samples=int(meta.sample_id.nunique()),
    runs=int(meta.run_id.nunique()),new_scvi_fit=False,new_clustering=False,new_neighbors=True,new_umap=True,
    neighbors=15,metric='euclidean',min_dist=.5,spread=1.0,threads=4,
    treatment_counts=meta.treatment.value_counts().to_dict(),run_counts=meta.run_id.value_counts().to_dict(),
    source_sha256=hashes,command=' '.join(sys.argv),pid=os.getpid())
(out/'STARTED.json').write_text(json.dumps(manifest,indent=2))
print('Retained 40,141 cells; building new graph and UMAP.',flush=True)
if not (out/'EMBEDDINGS_COMPLETE.json').exists():
 a=ad.AnnData(np.zeros((len(meta),1),dtype=np.float32),obs=meta.copy())
 a.obsm['X_scVI']=latent.to_numpy(dtype=np.float32)
 sc.pp.neighbors(a,use_rep='X_scVI',n_neighbors=15,metric='euclidean',random_state=42)
 sc.tl.umap(a,min_dist=.5,spread=1.,random_state=42)
 xy=a.obsm['X_umap'];assert xy.shape==(40141,2) and np.isfinite(xy).all()
 pd.DataFrame(xy,index=meta.index,columns=['UMAP_1','UMAP_2']).to_csv(out/'umap_PRIVATE.csv')
 sparse.save_npz(out/'connectivities_PRIVATE.npz',a.obsp['connectivities'])
 sparse.save_npz(out/'distances_PRIVATE.npz',a.obsp['distances'])
 meta.to_csv(out/'metadata_PRIVATE.csv');lab.to_csv(out/'state_labels_PRIVATE.csv')
 latent.to_csv(out/'scvi_latent_PRIVATE.csv')
 (out/'EMBEDDINGS_COMPLETE.json').write_text(json.dumps({**manifest,'elapsed_seconds':time.time()-start,'scanpy_version':sc.__version__},indent=2))
else:
 m=json.loads((out/'EMBEDDINGS_COMPLETE.json').read_text());assert m['source_sha256']==hashes
 xy=pd.read_csv(out/'umap_PRIVATE.csv',index_col=0).loc[meta.index].to_numpy()
print('Fresh retained-cell embedding complete; drawing three panels.',flush=True)
plt.rcParams.update({'font.size':15,'axes.titlesize':20,'axes.labelsize':16,'svg.fonttype':'none'})
fig=plt.figure(figsize=(23,12))
grid=fig.add_gridspec(2,3,height_ratios=[7,4.2],left=.045,right=.99,bottom=.025,top=.95,wspace=.22,hspace=.24)
axes=[fig.add_subplot(grid[0,i]) for i in range(3)]
legaxes=[fig.add_subplot(grid[1,i]) for i in range(3)]
for ax in legaxes:ax.axis('off')
state_names=sorted(lab.inferred_state.unique(),key=lambda x:int(x.split('_')[0][1:]))
state_colors={s:plt.get_cmap('tab20').colors[int(s.split('_')[0][1:])-1] for s in state_names}
treat_colors={'Intact':'#278c74','Castration':'#e99b27','Castration+ENZA':'#9569b7'}
run_names=sorted(meta.run_id.unique());run_colors=dict(zip(run_names,plt.get_cmap('tab10').colors))
order=np.random.default_rng(42).permutation(len(meta))
for i,(series,colors,title) in enumerate([(lab.inferred_state,state_colors,'Putative cell state'),(meta.treatment,treat_colors,'Treatment'),(meta.run_id,run_colors,'Library / run_id')]):
 ax=axes[i];ax.scatter(xy[order,0],xy[order,1],c=series.map(colors).to_numpy()[order],s=6,alpha=.9,linewidths=0,rasterized=True)
 ax.set(title=title,xlabel='UMAP_1',ylabel='UMAP_2',xticks=[],yticks=[])
 ax.spines[['top','right']].set_visible(False)
 handles=[]
 for key,col in colors.items():
  if i==0:
   rank,name=key.split('_',1)
   label=rank+' · '+name
  else:label=f'{key}  (n={int(series.eq(key).sum()):,})'
  handles.append(Line2D([],[],marker='o',linestyle='',color=col,markersize=8,label=label))
 legaxes[i].legend(handles=handles,loc='upper left',bbox_to_anchor=(0,1),frameon=False,
     fontsize=11.5 if i==0 else 15,ncol=2 if i==0 else 1,labelspacing=.8,handletextpad=.5,columnspacing=.8,borderaxespad=0)
# Short state IDs only, avoiding long overlapping names inside the UMAP.
texts=[]
for name in state_names:
 center=np.median(xy[lab.inferred_state.eq(name)],axis=0)
 texts.append(axes[0].text(*center,name.split('_')[0],ha='center',va='center',fontsize=13,weight='bold',
                         path_effects=[pe.withStroke(linewidth=3,foreground='white')]))
# No overall title/subtitle: scientific context belongs in the HTML caption.
fig.canvas.draw();renderer=fig.canvas.get_renderer()
# Check text and legend geometry in rendered display coordinates.
overlaps=[(i,j) for i in range(len(texts)) for j in range(i) if texts[i].get_window_extent(renderer).overlaps(texts[j].get_window_extent(renderer))]
for step in range(100):
 overlaps=[(i,j) for i in range(len(texts)) for j in range(i) if texts[i].get_window_extent(renderer).expanded(1.12,1.15).overlaps(texts[j].get_window_extent(renderer))]
 if not overlaps:break
 for i,j in overlaps:
  for idx,sign in [(i,1),(j,-1)]:
   display=axes[0].transData.transform(texts[idx].get_position())
   texts[idx].set_position(axes[0].transData.inverted().transform(display+np.array([0,sign*6])))
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
assert not [(i,j) for i in range(len(texts)) for j in range(i) if texts[i].get_window_extent(renderer).overlaps(texts[j].get_window_extent(renderer))], 'State text overlaps'
for plot_ax,legend_ax in zip(axes,legaxes):
 b=legend_ax.get_legend().get_window_extent(renderer)
 assert not b.overlaps(plot_ax.get_window_extent(renderer))
 assert fig.bbox.contains(b.x0,b.y0) and fig.bbox.contains(b.x1,b.y1), 'Legend clipping'
fig.savefig(out/'retained_state_treatment_run_umap.png',dpi=180)
fig.savefig(out/'retained_state_treatment_run_umap.pdf')
plt.close(fig)
assert hashes=={f:hashlib.sha256((source/f).read_bytes()).hexdigest() for f in protected}
oldxy=pd.read_csv(source/'umap_PRIVATE.csv',index_col=0).loc[meta.index].to_numpy()
assert not np.allclose(oldxy,xy), 'New UMAP must not be a cropped old embedding'
(out/'COMPLETE.json').write_text(json.dumps({**manifest,'elapsed_seconds':time.time()-start,'plot_panels':['state','treatment','run_id'],'legend_outside':True,'label_overlap_check':True},indent=2))
print('PASS: 40,141 retained cells; fresh UMAP; three panels and external legends; scientific source inputs unchanged.',flush=True)

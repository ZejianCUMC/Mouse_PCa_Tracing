"""Transfer reviewed state families by overlap and fresh markers, not raw IDs."""
from pathlib import Path
import json,gzip,sys
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmwrite
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
p=Path(__file__).resolve().parents[1];o=p/'analysis/merged_20260923_seed42';d=o/'restored_old_c18_seed42/recluster_r1_seed42'
# Split suffixes preserve distinct Leiden memberships; no forced mergers.
spec=[(7,'Lum-general','c1','Transferred'),(5,'Lum-IFN','c2','Transferred; pathway name inherited, not retested'),
 (12,'Lum-basal-hybrid-A','c3','Split; basal-enriched'),(15,'Lum-basal-hybrid-B','c3','Provisional split; luminal-enriched, only 327 cells'),
 (2,'Lum-basal-hybrid-Inflam','c4','Transferred; pathway name inherited, not retested'),
 (3,'Lum-Proximal','c6','Mixed: prior c6/c3/c2; provisional'),
 (0,'Basal-general','c5','Transferred'),(13,'Basal-general/Squamous','c5','Mixed prior c5/c7; no single previous identity'),
 (6,'Basal-Squamous','c7','Transferred'),(9,'Partial-EMT-A','c8','Split'),(14,'Partial-EMT-B','c8','Split; only 3 eligible samples'),
 (1,'Mes-Collagen','c9','Transferred; includes restored cells'),(8,'Mes-ECMremodeling','c10','Transferred'),
 (10,'Mes-GlycogenMetabolism','c11','Transferred; pathway name inherited, only 2 eligible samples'),
 (11,'Mes-Chemotactic','c12','Proposed ligand-based name; only 2 eligible samples, broad chemokine pathway not enriched'),
 (4,'Mes-Actomyosin','old-c18','Proposed process name; actomyosin/ROCK evidence, cycling retained')]
reviewed_merge='--merged' in sys.argv or '--reviewed10' in sys.argv
if reviewed_merge:
 d=d/'approved_merges';spec=json.loads((d/'NAME_SPEC.json').read_text())
 if '--reviewed10' in sys.argv:
  d=d/'reviewed10';spec=json.loads((d/'NAME_SPEC.json').read_text())
a=ad.read_h5ad(d/'clustered_graph_PRIVATE.h5ad');raw=a.obs.raw_cluster.astype(str)
summary=pd.read_csv(d/'profiles/cluster_summary.csv').set_index('cluster');rows=[];names={};ids={}
for i,(r,name,prev,note) in enumerate(spec,1):
 k=f'c{i}';ids[str(r)]=k;names[str(r)]=f'{k}_{name}';sub=a.obs[raw.eq(str(r))];q=summary.loc[r]
 markers=pd.read_csv(d/f'profiles/cluster_{r}_allgenes.csv');markers=markers[~markers.gene.str.match(r'^(Rpl|Rps|mt-|Gm\d)') & ~markers.gene.isin(['GFP','mCherry','Cas9','Cre'])]
 markers=markers[(markers.pct_in>=10)&(markers.log2fc_mean>=.5)].head(6)
 overlap=sub.previous_group.astype(str).value_counts()
 rows.append(dict(cluster=k,putative_name=name,cells=len(sub),Npp53=int(sub.model.eq('Npp53').sum()),Npp53Rb1=int(sub.model.eq('Npp53Rb1').sum()),raw_cluster=r,previous_family=prev,previous_overlap='; '.join(f'{x}: {n} ({100*n/len(sub):.1f}%)' for x,n in overlap.head(3).items()),representative_genes=', '.join(markers.gene),eligible_samples=int(q.eligible_samples),naming_status=note))
a.obs['final_cluster']=raw.map(ids);a.obs['putative_state']=raw.map(names)
assert a.n_obs==55133 and a.obs.final_cluster.nunique()==len(spec) and not a.obs.putative_state.isna().any()
a.obs.to_csv(d/'final_metadata_PRIVATE.csv');a.write_h5ad(d/'named_graph_PRIVATE.h5ad',compression='gzip')
pd.DataFrame(rows).to_csv(d/'state_summary.csv',index=False)
if reviewed_merge:
 cross=pd.read_csv(d/'merge_crosswalk.csv')
 cross['putative_name']=cross.new_cluster.map({r['cluster']:r['putative_name'] for r in rows})
 cross.to_csv(d/'merge_crosswalk.csv',index=False)
pd.crosstab(a.obs.final_cluster,a.obs.previous_group).reindex([r['cluster'] for r in rows]).to_csv(d/'previous_state_crosswalk.csv')
a.obs.groupby(['final_cluster','putative_state','model','sample_id','run_id','treatment'],observed=True).size().rename('n_cells').reset_index().to_csv(d/'state_sample_counts.csv',index=False)
(d/'PUTATIVE_NAMES.json').write_text(json.dumps(rows,indent=2))
pal={v:matplotlib.colors.to_hex(plt.get_cmap('tab20')(i)) for i,v in enumerate(names.values())}
tr={'Intact':'#4c78a8','Castration':'#f58518','Castration+ENZA':'#54a24b'};rp={k:matplotlib.colors.to_hex(plt.get_cmap('tab20')(i)) for i,k in enumerate(sorted(a.obs.run_id.unique()))}
xy=a.obsm['X_umap'];draw=np.random.default_rng(42).permutation(len(a));fig,axes=plt.subplots(1,3,figsize=(23,11))
for ax,col,palette in zip(axes,['putative_state','treatment','run_id'],[pal,tr,rp]):
 ax.scatter(xy[draw,0],xy[draw,1],c=[palette[v] for v in a.obs[col].iloc[draw]],s=3,alpha=.85,linewidths=0,rasterized=True)
 ax.set_xlabel('UMAP_1',fontsize=13);ax.set_ylabel('UMAP_2',fontsize=13);ax.set_xticks([]);ax.set_yticks([])
 if col=='putative_state':
  for k in palette:ax.text(*np.median(xy[a.obs[col].eq(k)],axis=0),k.split('_')[0],fontsize=10,ha='center',weight='bold',bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=1))
 ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color=v,markersize=5,label=k) for k,v in palette.items()],loc='upper center',bbox_to_anchor=(.5,-.12),ncol=1 if col in ['putative_state','treatment'] else 2,frameon=False,fontsize=9)
fig.subplots_adjust(left=.04,right=.99,top=.98,bottom=.47,wspace=.22);fig.savefig(d/'reclustered_joint_umap.png',dpi=170,bbox_inches='tight');plt.close(fig)
if reviewed_merge:
 from plot_reviewed_umap import render
 render(a,d/'reclustered_joint_umap.png')
c=o/'filtered_no_c17_c18_seed42/resolution_1p0/final_cleanup';genes=pd.read_csv(c/'dotplot_genes.csv').gene.tolist()
if reviewed_merge:
 # Keep the reviewed epithelial order; add discriminating programs, not only shared ECM genes.
 additions={'Cyp2f2':['Rsad2','Isg15','Irf7'], 'Aspn':['Sparcl1'], 'Fmod':['Cdh13'],
            'Crabp1':['Cxcl12','Ccl2','Ccl7','Il6','Cxcl1','Acta2','Tagln','Myl9','Tpm2']}
 genes=[g for old in genes for g in [old]+additions.get(old,[])]
 genes+=['Sox2','Syp','Chga']
 if '--reviewed10' in sys.argv:
  genes=[g for g in genes if g not in ['Fst','Cst6','Arg1']]
  ix=genes.index('Sprr2d');genes[ix:ix]=['Fst','Cst6','Arg1','Krt6a','Krt16']
  mes=['Aspn','Sparcl1','Col3a1','Col1a1','Col5a2','Mfap2',
       'Dcn','Igfbp2','Fmod','Cdh13','Mmp3','Prelp','Has2',
       'Tmeff2','Slc13a5','Pygl','Gbe1',
       'Acta2','Tagln','Myl9','Tpm2','Top2a','Pclaf','Birc5','Ube2c','Mki67']
  genes=genes[:genes.index('Aspn')]+mes+['Sox2','Syp','Chga']
 assert len(genes)==len(set(genes))
counts=ad.read_h5ad(o/'shared_counts.h5ad');ix=pd.Index(counts.var.gene_symbol).get_indexer(genes);assert (ix>=0).all()
X=counts[a.obs_names,ix].X;keep=np.asarray((X>0).sum(0)).ravel()>0;genes=np.array(genes)[keep].tolist();X=X[:,keep]
with gzip.open(d/'dotplot_counts_PRIVATE.mtx.gz','wb') as f:mmwrite(f,X.T,field='integer')
pd.DataFrame({'gene':genes}).to_csv(d/'dotplot_genes.csv',index=False);pd.DataFrame({'gene':genes,'panel':'Markers'}).to_csv(d/'dotplot_panels.csv',index=False)
if reviewed_merge:
 boundaries={'Epcam':'Luminal','Rsad2':'ISG','Wfdc18':'Proximal-like','Calml3':'Basal / epithelial',
             'Fst':'Differentiation','Vim':'EMT','Aspn':'Matrix / ECM','Has2':'Metabolism / state',
             'Cxcl12':'Inflammatory','Acta2':'Contractile','Top2a':'Cycling','Sox2':'Plasticity / NE'}
 if '--reviewed10' in sys.argv:
  boundaries.pop('Fst');boundaries['Krt6a']='Squamous'
  boundaries.pop('Has2');boundaries.pop('Cxcl12')
  boundaries.update({'Aspn':'Matrix / CAF-like','Dcn':'ECM remodeling',
                     'Tmeff2':'c9 markers / metabolism',
                     'Acta2':'Contractile'})
 groups=[];group=''
 for gene in genes:
  group=boundaries.get(gene,group);groups.append(group)
 pd.DataFrame({'gene':genes,'panel':groups}).to_csv(d/'dotplot_panels.csv',index=False)
pd.DataFrame({'state':list(names.values())}).to_csv(d/'state_order.csv',index=False)
print('Prepared',len(spec),'named clusters and',len(genes),'genes',flush=True)

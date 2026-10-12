"""Post-rerun marker-first review. Requires freshly reviewed cluster mapping JSON."""
import os,json
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import jensenshannon
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(os.environ['RB1_RUN']);r=out/'review';r.mkdir(exist_ok=True)
config=json.loads((out/'reviewed_cluster_mapping.json').read_text())
a=sc.read_h5ad(out/'major_types_full.h5ad');o=a.obs.copy();cl=o.primary_cluster.astype(str)
assert a.n_obs==56686 and o.sample_id.nunique()==22
assert not o.sample_id.isin(['JZ257','JZ275','JZ278','JZ281']).any()
assert o.groupby('sample_id',observed=True).size().min()>=400
selection=json.loads((out/'integration_selection.json').read_text())
assert config['method']==selection['candidate']
mapping=config['clusters'];assert set(cl.unique())==set(mapping)
o['marker_first_major_type']=cl.map({k:v['major_type'] for k,v in mapping.items()})
o['QC_review_flag']=cl.map({k:v['qc_flag'] for k,v in mapping.items()})
o['cycling_modifier']=cl.map({k:v.get('cycling',False) for k,v in mapping.items()})
o['all_reporter_positive_interest']=o.reporter_positive
for g in ['GFP','mCherry']:
    values=np.asarray(a.layers['counts'][:,a.var_names.get_loc(g)].toarray()).ravel()
    assert np.array_equal(values,o['raw_'+g].to_numpy())
assert np.array_equal(o.reporter_positive,((o.raw_GFP>0)|(o.raw_mCherry>0)))
assert np.array_equal(o.reporter_ge2,((o.raw_GFP>=2)|(o.raw_mCherry>=2)))
assert int(o.reporter_positive.sum())==39812 and int(o.reporter_ge2.sum())==36530
o['GFP_only']=(o.raw_GFP>0)&(o.raw_mCherry==0)
o['mCherry_only']=(o.raw_mCherry>0)&(o.raw_GFP==0)
o['double_positive']=(o.raw_GFP>0)&(o.raw_mCherry>0)
o.to_csv(r/'marker_first_all_cell_annotations.csv')
o.loc[o.reporter_positive].to_csv(r/'ALL_reporter_positive_cells_for_review.csv')
o.loc[o.reporter_ge2].to_csv(r/'ALL_reporter_ge2_cells_for_review.csv')
def summarize(keys,name):
    summary=o.groupby(keys,observed=True).agg(cells=('reporter_positive','size'),reporter_positive=('reporter_positive','sum'),reporter_ge2=('reporter_ge2','sum'),GFP_only=('GFP_only','sum'),mCherry_only=('mCherry_only','sum'),double_positive=('double_positive','sum'))
    summary['reporter_percent']=100*summary.reporter_positive/summary.cells
    assert summary.cells.sum()==a.n_obs
    summary.to_csv(r/name)
summarize(['marker_first_major_type'],'major_types_then_reporters.csv')
summarize(['primary_cluster','marker_first_major_type'],'clusters_then_reporters.csv')
summarize(['sample_id','marker_first_major_type'],'sample_type_reporters.csv')
summarize(['sample_id'],'sample_reporters.csv')
summary=o.groupby('primary_cluster',observed=True).agg(cells=('reporter_positive','size'),reporter_positive=('reporter_positive','sum'),median_genes=('n_genes_by_counts','median'),median_mt=('pct_counts_mt','median'))
summary['major_type']=summary.index.astype(str).map({k:v['major_type'] for k,v in mapping.items()})
summary['QC_flag']=summary.index.astype(str).map({k:v['qc_flag'] for k,v in mapping.items()})
summary['marker_reason']=summary.index.astype(str).map({k:v['reason'] for k,v in mapping.items()})
summary['dominant_mouse']=o.groupby('primary_cluster',observed=True).sample_id.agg(lambda x:x.value_counts().index[0])
summary['dominant_mouse_fraction']=o.groupby('primary_cluster',observed=True).sample_id.agg(lambda x:x.value_counts(normalize=True).iloc[0])
summary.to_csv(r/'reviewed_cluster_evidence.csv')
methods=selection['methods'];balance=[];support=[]
for (label,treatment),group in o[o.expression_anchor!='Unresolved'].groupby(['expression_anchor','treatment'],observed=True):
    bc=group.batch.value_counts();valid=bc[bc>=20].index;group=group[group.batch.isin(valid)]
    if len(valid)<2 or len(group)<100:continue
    support.append(dict(label=label,treatment=treatment,batches=';'.join(sorted(valid)),cells=len(group),min_batch_cells=int(bc.loc[valid].min())))
    for seed in [19,29,39]:
        rng=np.random.default_rng(seed);n=min(150,int(bc.loc[valid].min()))
        names=np.concatenate([rng.choice(group.index[group.batch==b].to_numpy(),n,replace=False) for b in sorted(valid)])
        ix=a.obs_names.get_indexer(names);b=o.loc[names,'batch'].to_numpy();expected=np.repeat(1/len(valid),len(valid))
        for method in methods:
            latent=a.obsm['X_'+method][ix]
            nn=NearestNeighbors(n_neighbors=min(16,len(ix))).fit(latent).kneighbors(latent,return_distance=False)[:,1:]
            local=np.column_stack([(b[nn]==v).mean(axis=1) for v in sorted(valid)])
            score=1-np.mean([jensenshannon(row,expected,base=2)**2 for row in local])
            balance.append(dict(method=method,label=label,treatment=treatment,seed=seed,cells=len(ix),one_minus_JS=float(score)))
pd.DataFrame(support).to_csv(r/'balanced_mixing_eligible_strata.csv',index=False)
balance=pd.DataFrame(balance);balance.to_csv(r/'balanced_mixing_details.csv',index=False)
balance.groupby(['method','seed']).one_minus_JS.mean().unstack().to_csv(r/'balanced_mixing_by_seed.csv')
a.obs=o;a.obsm['X_umap']=a.obsm['X_umap_'+config['method']].copy()
for g in ['GFP','mCherry']:a.obs['log1p_'+g]=np.log1p(o['raw_'+g])
fig,axes=plt.subplots(1,2,figsize=(22,8))
for field,ax in zip(['marker_first_major_type','primary_cluster'],axes):sc.pl.umap(a,color=field,ax=ax,show=False,size=3)
fig.tight_layout();fig.savefig(r/'major_types.png',dpi=140);plt.close(fig)
fig,axes=plt.subplots(2,3,figsize=(23,13))
for field,ax in zip(['log1p_GFP','log1p_mCherry','all_reporter_positive_interest','QC_review_flag','n_genes_by_counts','pct_counts_mt'],axes.flat):
    sc.pl.umap(a,color=field,ax=ax,show=False,size=3,**({'vmax':'p99'} if field.startswith('log1p') else {}))
fig.tight_layout();fig.savefig(r/'reporters_QC.png',dpi=140);plt.close(fig)
panels={'Epithelial':['Krt5','Krt14','Trp63','Krt8','Krt18','Krt19','Epcam','Cdh1'],
 'Mesenchymal':['Vim','Zeb1','Snai2','Col1a1','Col3a1','Dcn','Lum','Pdgfra','Acta2','Tagln'],
 'Immune / endothelial':['Ptprc','Lyz2','Tyrobp','Cd3e','Ms4a1','Pecam1','Cdh5'],
 'Cycling / reporters':['Mki67','Top2a','GFP','mCherry']}
panels={k:[g for g in v if g in a.var_names] for k,v in panels.items()}
sc.pl.dotplot(a,panels,groupby='primary_cluster',standard_scale='var',show=False)
plt.savefig(r/'marker_dotplot.png',dpi=150,bbox_inches='tight');plt.close('all')
genes=list(dict.fromkeys(sum(panels.values(),[])));expr=a[:,genes].to_df()
expr.groupby(cl,observed=True).mean().to_csv(r/'marker_means.csv')
(expr>0).groupby(cl,observed=True).mean().to_csv(r/'marker_fractions.csv')
coexp=((expr.Krt8>0)|(expr.Krt18>0))&(expr.Vim>0)
coexp.groupby(cl,observed=True).mean().to_csv(r/'keratin_vimentin_coexpression.csv')
pd.crosstab(o.sample_id,o.marker_first_major_type,normalize='index').plot.bar(stacked=True,figsize=(17,7),width=.9)
plt.ylabel('Fraction per mouse');plt.legend(bbox_to_anchor=(1,1),fontsize=8);plt.tight_layout();plt.savefig(r/'composition_by_mouse.png',dpi=140);plt.close('all')
a.write(r/'reviewed_all_cells.h5ad');a[o.reporter_positive].write(r/'ALL_reporter_positive_cells.h5ad')
result={'all_cells':a.n_obs,'mice':int(o.sample_id.nunique()),'all_reporter_positive':int(o.reporter_positive.sum()),'all_reporter_ge2':int(o.reporter_ge2.sum()),'method':config['method'],'excluded_samples':['JZ257','JZ275','JZ278','JZ281'],'marker_review':config,'status':'review computations complete; inspect balanced diagnostics and figures before publication'}
(r/'review_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));print('REVIEW COMPLETE')

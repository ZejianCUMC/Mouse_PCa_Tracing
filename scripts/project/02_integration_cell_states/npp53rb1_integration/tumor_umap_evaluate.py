"""UMAP-only comparison. No Leiden, differential expression or new state annotation."""
import os,json
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

out=Path(os.environ['RB1_RUN'])/'tumor_umap_20260921'
obs=pd.read_csv(out/'metadata.csv',index_col=0)
programs=pd.read_csv(out/'expression_programs.csv',index_col=0).loc[obs.index]
methods=['uncorrected','Harmony_run','Harmony_sample','Harmony_run_sample','scVI_run','scVI_sample','scVI_run_sample']
names={'uncorrected':'Uncorrected PCA','Harmony_run':'Harmony: run (theta 1)',
       'Harmony_sample':'Harmony: sample (theta 1)','Harmony_run_sample':'Harmony: run 1 + sample 0.25',
       'scVI_run':'scVI NB: run','scVI_sample':'scVI NB: sample','scVI_run_sample':'scVI NB: run + sample covariate'}
a=ad.AnnData(X=np.zeros((len(obs),1),dtype=np.float32),obs=obs.copy())
assert a.n_obs==43057 and not obs.parent_harmony_cluster.astype(str).isin(['4','5','6','7']).any()
# Independent-of-integration expression anchors inherited from the count-based panel.
anchor_names=['Basal_epithelial','Luminal_epithelial','Mesenchymal_fibroblast']
conf=obs.expression_anchor.isin(anchor_names)
rng=np.random.default_rng(19)
eval_ix=np.concatenate([rng.choice(np.flatnonzero(obs.expression_anchor.to_numpy()==v),
    min(1000,int((obs.expression_anchor==v).sum())),replace=False) for v in anchor_names])
anchors=obs.expression_anchor.to_numpy()
strata=[];support=[]
for batchkey in ['run_id','sample_id']:
    for (label,treat),g in obs.loc[conf].groupby(['expression_anchor','treatment'],observed=True):
        counts=g[batchkey].value_counts();valid=sorted(counts[counts>=20].index)
        if len(valid)<2:continue
        n=min(100,int(counts.loc[valid].min()))
        support.append(dict(batch_key=batchkey,anchor=label,treatment=treat,batches=len(valid),per_batch=n))
        for seed in [19,29,39]:
            rg=np.random.default_rng(seed)
            ids=np.concatenate([rg.choice(g.index[g[batchkey]==b].to_numpy(),n,replace=False) for b in valid])
            ix=obs.index.get_indexer(ids)
            strata.append((batchkey,label,treat,seed,ix,obs.iloc[ix][batchkey].to_numpy(),valid))
pd.DataFrame(support).to_csv(out/'mixing_support.csv',index=False)
metrics=[];mixrows=[];coord={}
for method in methods:
    print('Evaluate',method,flush=True)
    frame=pd.read_csv(out/f'latent_{method}.csv',index_col=0)
    assert frame.index.is_unique and set(frame.index)==set(obs.index)
    latent=frame.loc[obs.index].to_numpy();assert latent.shape==(43057,30) and np.isfinite(latent).all()
    a.obsm['X_'+method]=latent
    sc.pp.neighbors(a,use_rep='X_'+method,n_neighbors=15,random_state=19)
    sc.tl.umap(a,random_state=19,min_dist=.3)
    coord[method]=a.obsm['X_umap'].copy();a.obsm['X_umap_'+method]=coord[method]
    pd.DataFrame(coord[method],index=obs.index,columns=['UMAP_1','UMAP_2']).to_csv(out/f'UMAP_{method}.csv')
    asw=(silhouette_score(latent[eval_ix],anchors[eval_ix])+1)/2
    confident_ix=np.flatnonzero(conf)
    nn=NearestNeighbors(n_neighbors=16,n_jobs=8).fit(latent[confident_ix]).kneighbors(return_distance=False)[:,:15]
    al=anchors[confident_ix]
    purity=pd.Series((al[nn]==al[:,None]).mean(axis=1)).groupby(al).mean().mean()
    fullnn=NearestNeighbors(n_neighbors=16,n_jobs=8).fit(latent).kneighbors(return_distance=False)[:,:15]
    program_r={k:float(spearmanr(programs[k],programs[k].to_numpy()[fullnn].mean(axis=1)).statistic) for k in programs}
    values={'run_id':[],'sample_id':[]}
    for batchkey,label,treat,seed,ix,b,valid in strata:
        nbr=NearestNeighbors(n_neighbors=16,n_jobs=8).fit(latent[ix]).kneighbors(return_distance=False)[:,:15]
        local=np.column_stack([(b[nbr]==v).mean(axis=1) for v in valid])
        expected=np.repeat(1/len(valid),len(valid))
        score=float(1-np.mean([jensenshannon(v,expected,base=2)**2 for v in local]))
        values[batchkey].append(score)
        mixrows.append(dict(method=method,batch_key=batchkey,anchor=label,treatment=treat,seed=seed,score=score))
    bio=float((asw+purity)/2)
    metrics.append(dict(method=method,anchor_ASW=float(asw),anchor_neighbor_purity=float(purity),
        biological_conservation_proxy=bio,run_mixing=float(np.mean(values['run_id'])),
        sample_mixing=float(np.mean(values['sample_id'])),program_neighbor_correlation=float(np.mean(list(program_r.values()))),
        **{k+'_correlation':v for k,v in program_r.items()}))
    pd.DataFrame(metrics).to_csv(out/'integration_diagnostics.csv',index=False)
    print(metrics[-1],flush=True)
df=pd.DataFrame(metrics).set_index('method')
baseline=df.loc['uncorrected']
df['bio_guardrail_pass']=df.biological_conservation_proxy>=baseline.biological_conservation_proxy-.05
df['program_guardrail_pass']=df.program_neighbor_correlation>=baseline.program_neighbor_correlation-.05
df['run_mixing_gain']=df.run_mixing-baseline.run_mixing
df['diagnostic_score']=.6*df.biological_conservation_proxy+.4*df.run_mixing
eligible=df.drop(index='uncorrected').query('bio_guardrail_pass and program_guardrail_pass and run_mixing_gain>0')
chosen=eligible.diagnostic_score.idxmax() if len(eligible) else 'uncorrected'
# Prefer less biologically confounded run-only adjustment when diagnostics are near-tied.
conservative=eligible.loc[eligible.index.isin(['Harmony_run','scVI_run'])]
if len(conservative) and conservative.diagnostic_score.max()>=df.loc[chosen,'diagnostic_score']-.02:
    chosen=conservative.diagnostic_score.idxmax()
df.to_csv(out/'integration_diagnostics.csv')
pd.DataFrame(mixrows).to_csv(out/'mixing_details.csv',index=False)
pd.DataFrame(mixrows).groupby(['method','batch_key','seed']).score.mean().to_csv(out/'mixing_by_seed.csv')
plt.rcParams.update({'font.size':14,'axes.titlesize':17,'axes.labelsize':15})
order=np.random.default_rng(19).permutation(len(obs))
def categorical(ax,xy,field,title,size=5,legend=False):
    categories=sorted(obs[field].astype(str).unique())
    if field=='treatment':categories=['Intact','Castration','Castration+ENZA']
    palette=plt.get_cmap('tab20' if len(categories)>10 else 'tab10')
    colours=list(plt.get_cmap('tab20').colors)+list(plt.get_cmap('Dark2').colors) if len(categories)>20 else [palette(i) for i in range(len(categories))]
    cols={v:colours[i] for i,v in enumerate(categories)}
    ax.scatter(*xy[order].T,s=size,c=[cols[str(v)] for v in obs[field].iloc[order]],linewidths=0,rasterized=True)
    ax.set(title=title,xlabel='UMAP_1',ylabel='UMAP_2',xticks=[],yticks=[],aspect='equal')
    ax.spines[['top','right','bottom','left']].set_visible(False)
    if legend:
        ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color=cols[v],label=v) for v in categories],
                  loc='upper left',bbox_to_anchor=(1,1),fontsize=10,frameon=False)
    return cols
for field in ['run_id','sample_id','treatment','parent_harmony_cluster']:
    fig,axes=plt.subplots(2,4,figsize=(28,14))
    for method,ax in zip(methods,axes.flat):cols=categorical(ax,coord[method],field,names[method])
    axes.flat[-1].axis('off')
    axes.flat[-1].legend(handles=[Line2D([],[],marker='o',linestyle='',color=c,label=v) for v,c in cols.items()],
                        loc='center',frameon=False,ncol=2 if len(cols)>10 else 1,fontsize=13)
    fig.suptitle(f'Tumour-only reintegration | 43,057 cells | colour: {field}',fontsize=23)
    fig.tight_layout();fig.savefig(out/f'comparison_{field}.png',dpi=150);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(23,8))
for field,ax in zip(['run_id','sample_id','treatment'],axes):
    categorical(ax,coord[chosen],field,field,size=7,legend=True)
fig.suptitle('Provisional display: '+names[chosen]+' | 43,057 tumour-classified cells',fontsize=22)
fig.tight_layout();fig.savefig(out/'selected_umap.png',dpi=170);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(23,8))
xy=coord[chosen]
for t,ax in zip(['Intact','Castration','Castration+ENZA'],axes):
    ax.scatter(*xy[order].T,s=4,c='#edf0f3',linewidths=0)
    subset=obs.treatment.to_numpy()==t
    cats=sorted(obs.run_id.unique());colmap={v:plt.get_cmap('tab10')(i) for i,v in enumerate(cats)}
    ix=order[subset[order]]
    ax.scatter(*xy[ix].T,s=7,c=[colmap[v] for v in obs.run_id.iloc[ix]],linewidths=0)
    ax.set(title=f'{t}\n{subset.sum():,} cells',xlabel='UMAP_1',ylabel='UMAP_2',xticks=[],yticks=[],aspect='equal')
    ax.spines[['top','right','bottom','left']].set_visible(False)
fig.suptitle(names[chosen]+' | shared tumour-only UMAP split by treatment',fontsize=22)
fig.tight_layout();fig.savefig(out/'selected_by_treatment.png',dpi=170);plt.close(fig)
# Preserve all latent spaces/UMAPs, but no new clusters or cell-state calls.
a.obsm['X_umap']=coord[chosen];a.obsp.clear();a.uns.pop('neighbors',None)
a.uns['display_method']=chosen;a.write_h5ad(out/'tumor_umap_comparison.h5ad')
selection={'candidate':chosen,'display_name':names[chosen],'cells':43057,'samples':22,
 'methods':methods,'rule':'60% anchor conservation +40% within-anchor/treatment run mixing; reject >0.05 biological/program loss or no mixing gain. Prefer run-only if within 0.02 of best score.',
 'limitations':'Exploratory custom proxies, not independent ground truth or scIB; sample nested in run and treatment. Perfect mixing can remove genuine clone/mouse differences. One model seed; mixing subsampling seeds 19/29/39 are not inferential replicates.',
 'status':'Computation complete; visual review required before publication','stop':'UMAP only; no new clustering, DE, state annotation or trajectories'}
(out/'COMPLETED.json').write_text(json.dumps(selection,indent=2));print(json.dumps(selection),flush=True)

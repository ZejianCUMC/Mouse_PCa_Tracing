"""Transparent preliminary benchmark; expression anchors are NOT ground-truth labels."""
import os, json
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import jensenshannon
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(os.environ['RB1_RUN']); rng=np.random.default_rng(19)
a=sc.read_h5ad(out/'merged_counts.h5ad')
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
panels={
 'Basal_epithelial':['Krt5','Krt14','Trp63','Krt17'],
 'Luminal_epithelial':['Krt8','Krt18','Krt19','Epcam','Cd24a'],
 'Mesenchymal_fibroblast':['Col1a1','Col1a2','Col3a1','Dcn','Lum','Pdgfra'],
 'Myeloid':['Lyz2','Tyrobp','Csf1r','Adgre1','Aif1'],
 'T_NK':['Cd3d','Cd3e','Trbc2','Nkg7','Klr d1'.replace(' ','')],
 'B_cell':['Cd79a','Ms4a1','Cd79b','Cd19'],
 'Endothelial':['Pecam1','Cdh5','Kdr','Emcn','Vwf'],
 'Pericyte_smooth_muscle':['Rgs5','Pdgfrb','Acta2','Tagln','Myh11'],
 'Neuroendocrine':['Syp','Chga','Chgb','Ascl1','Insm1'],
}
missing={k:[g for g in v if g not in a.var_names] for k,v in panels.items()}
panels={k:[g for g in v if g in a.var_names] for k,v in panels.items()}
panels={k:v for k,v in panels.items() if len(v)>=2}
(out/'marker_panel_coverage.json').write_text(json.dumps({'used':panels,'missing':missing},indent=2))
markers=list(dict.fromkeys(sum(panels.values(),[])+[g for g in ['Ptprc','Vim','Zeb1','Zeb2','Snai2','Fn1','Mki67','Top2a','Malat1','Neat1','GFP','mCherry'] if g in a.var_names]))
expr=a[:,markers].to_df()
scores={}
for label,genes in panels.items():
 sc.tl.score_genes(a,genes,score_name='score_'+label,random_state=19,use_raw=False)
 scores[label]=a.obs['score_'+label].to_numpy()
scores=pd.DataFrame(scores,index=a.obs_names)
rank=np.argsort(-scores.to_numpy(),axis=1); labels=np.asarray(scores.columns)[rank[:,0]]
top=scores.to_numpy()[np.arange(a.n_obs),rank[:,0]]
margin=top-scores.to_numpy()[np.arange(a.n_obs),rank[:,1]]
supported=np.array([sum(expr.iloc[i][panels[label]]>0)>=2 for i,label in enumerate(labels)])
confident=(top>=.25)&(margin>=.25)&supported
a.obs['expression_anchor']=np.where(confident,labels,'Unresolved')
a.obs['anchor_score']=top; a.obs['anchor_margin']=margin
scores.to_csv(out/'expression_program_scores.csv')
a.obs[['expression_anchor','anchor_score','anchor_margin']].to_csv(out/'expression_anchors.csv')
methods=['uncorrected','Harmony_theta1','Harmony_theta2','Harmony_theta4','fastMNN','scVI_nb','scVI_zinb']
latents={}; coords={}; metricrows=[]; details=[]
# Fixed, method-independent anchors and evaluation cells across all methods.
anchor_ix=np.flatnonzero(confident)
eval_ix=rng.choice(anchor_ix,min(5000,len(anchor_ix)),replace=False)
for name in methods:
 frame=pd.read_csv(out/f'latent_{name}.csv',index_col=0)
 assert frame.index.is_unique and set(frame.index)==set(a.obs_names)
 latent=frame.loc[a.obs_names].to_numpy(); assert np.isfinite(latent).all()
 latents[name]=latent
 a.obsm['comparison_latent']=latent
 sc.pp.neighbors(a,use_rep='comparison_latent',n_neighbors=15,random_state=19)
 sc.tl.umap(a,random_state=19); coords[name]=a.obsm['X_umap'].copy()
 sc.tl.leiden(a,resolution=.5,key_added='clusters_'+name,flavor='igraph',n_iterations=2,directed=False,random_state=19)
 pd.DataFrame(coords[name],index=a.obs_names,columns=['UMAP1','UMAP2']).to_csv(out/f'UMAP_{name}.csv')
 label_asw=np.nan; purity=np.nan
 evlabels=a.obs.expression_anchor.iloc[eval_ix].to_numpy()
 if len(np.unique(evlabels))>1 and len(eval_ix)>10:
  label_asw=float((silhouette_score(latent[eval_ix],evlabels)+1)/2)
  nbr=NearestNeighbors(n_neighbors=min(16,len(anchor_ix))).fit(latent[anchor_ix]).kneighbors(latent[anchor_ix],return_distance=False)[:,1:]
  al=a.obs.expression_anchor.iloc[anchor_ix].to_numpy()
  percell=(al[nbr]==al[:,None]).mean(axis=1)
  purity=float(pd.Series(percell).groupby(al).mean().mean())
 mix=[]
 for (celltype,treatment),group in a.obs.loc[confident].groupby(['expression_anchor','treatment'],observed=True):
  counts=group.batch.value_counts(); allowed=counts[counts>=20].index
  group=group[group.batch.isin(allowed)]
  if len(allowed)<2 or len(group)<100:continue
  ix=a.obs_names.get_indexer(group.index); b=group.batch.to_numpy()
  nn=NearestNeighbors(n_neighbors=min(31,len(ix))).fit(latent[ix]).kneighbors(latent[ix],return_distance=False)[:,1:]
  categories=sorted(allowed); expected=np.array([(b==x).mean() for x in categories])
  local=np.column_stack([(b[nn]==x).mean(axis=1) for x in categories])
  js=np.array([jensenshannon(row,expected,base=2)**2 for row in local])
  score=float(1-js.mean()); mix.append(score)
  details.append(dict(method=name,celltype=celltype,treatment=treatment,cells=len(ix),batches=len(allowed),one_minus_JS=score))
 bio=float(np.nanmean([label_asw,purity]))
 metricrows.append(dict(method=name,anchor_label_ASW=label_asw,anchor_neighbor_purity=purity,
  within_type_treatment_batch_mixing=float(np.mean(mix)) if mix else np.nan,
  eligible_strata=len(mix),biological_conservation_proxy=bio))
metrics=pd.DataFrame(metricrows).set_index('method')
baseline=metrics.loc['uncorrected','biological_conservation_proxy']
metrics['bio_guardrail_pass']=metrics.biological_conservation_proxy>=baseline-.05
metrics['mixing_gain_vs_uncorrected']=metrics.within_type_treatment_batch_mixing-metrics.loc['uncorrected','within_type_treatment_batch_mixing']
metrics['diagnostic_score']=.6*metrics.biological_conservation_proxy+.4*metrics.within_type_treatment_batch_mixing
eligible=metrics.drop(index='uncorrected').query('bio_guardrail_pass and mixing_gain_vs_uncorrected > 0')
chosen=eligible.diagnostic_score.idxmax() if len(eligible) else 'uncorrected'
metrics.to_csv(out/'integration_benchmark.csv'); pd.DataFrame(details).to_csv(out/'integration_benchmark_strata.csv',index=False)
selection={'candidate':chosen,'status':'provisional; inspect markers, QC and old/new mixing before final acceptance',
 'rule':'60% anchor bio-conservation proxy + 40% within-type/treatment batch-mixing diagnostic; reject >0.05 bio loss or no mixing gain over uncorrected',
 'limitations':'Expression anchors are not independent ground truth. JS diagnostic is custom, not scIB/iLISI/kBET. Labels and thresholds are heuristic. No treatment correction. Library-treatment confounding remains.',
 'anchor_cells':int(confident.sum()),'all_cells':a.n_obs,'methods':methods}
(out/'integration_selection.json').write_text(json.dumps(selection,indent=2))
clusters=a.obs['clusters_'+chosen]
means=expr.groupby(clusters,observed=True).mean(); fractions=(expr>0).groupby(clusters,observed=True).mean()
means.to_csv(out/'cluster_marker_means.csv'); fractions.to_csv(out/'cluster_marker_fractions.csv')
cross=pd.crosstab(clusters,a.obs.expression_anchor); cross.to_csv(out/'cluster_anchor_evidence.csv')
mapping={}; evidence=[]
for cluster,row in cross.iterrows():
 known=row.drop(labels=['Unresolved'],errors='ignore'); total=row.sum()
 label=known.idxmax() if known.sum() else 'Unresolved'
 fraction=float(known.max()/total) if known.sum() else 0
 mapping[cluster]=label if fraction>=.5 else 'Mixed_or_unresolved'
 evidence.append(dict(cluster=cluster,candidate_type=label,anchor_fraction=fraction,preliminary_major_type=mapping[cluster],cells=int(total)))
pd.DataFrame(evidence).to_csv(out/'major_type_annotation_evidence.csv',index=False)
a.obs['major_type']=clusters.map(mapping).astype(str)
a.obs['primary_cluster']=clusters.astype(str)
tumor_compatible=a.obs.major_type.isin(['Basal_epithelial','Luminal_epithelial','Mesenchymal_fibroblast','Pericyte_smooth_muscle','Neuroendocrine'])
a.obs['putative_tumor']=np.where(a.obs.reporter_positive & tumor_compatible,'Reporter-positive putative tumor',
 np.where(a.obs.reporter_positive,'Reporter-positive: lineage/ambient review','Reporter-negative: not automatically host'))
for flag,name in [('reporter_positive','all_reporter_positive_cells'),('reporter_ge2','reporter_ge2_sensitivity_cells')]:
 a.obs.loc[a.obs[flag]].to_csv(out/f'{name}.csv')
a.obs.loc[a.obs.putative_tumor=='Reporter-positive putative tumor'].to_csv(out/'putative_tumor_cells.csv')
summary=a.obs.groupby(['primary_cluster','major_type'],observed=True).agg(cells=('reporter_positive','size'),reporter_fraction=('reporter_positive','mean'),reporter_ge2_fraction=('reporter_ge2','mean'),mean_GFP=('raw_GFP','mean'),mean_mCherry=('raw_mCherry','mean'),median_genes=('n_genes_by_counts','median'),median_mt=('pct_counts_mt','median'))
summary.to_csv(out/'cluster_reporter_QC_summary.csv')
composition=a.obs.major_type.value_counts().rename_axis('major_type').to_frame('cells')
composition['percent']=100*composition.cells/a.n_obs; composition.to_csv(out/'major_type_composition.csv')
for field in ['sample_id','batch','treatment','cohort']:
 pd.crosstab(a.obs[field],a.obs.major_type).to_csv(out/f'major_type_by_{field}.csv')
 a.obs.groupby(field,observed=True).agg(cells=('reporter_positive','size'),reporter_positive=('reporter_positive','sum'),reporter_ge2=('reporter_ge2','sum')).to_csv(out/f'reporters_by_{field}.csv')
for field in ['batch','sample_id','treatment','cohort','major_type']:
 fig,axes=plt.subplots(2,4,figsize=(28,13))
 for name,ax in zip(methods,axes.flat):
  a.obsm['X_umap']=coords[name]; sc.pl.umap(a,color=field,ax=ax,show=False,title=name,size=3)
 axes.flat[-1].axis('off'); fig.tight_layout(); fig.savefig(out/f'comparison_{field}.png',dpi=140); plt.close(fig)
a.obsm['X_umap']=coords[chosen]
fig,axes=plt.subplots(2,3,figsize=(23,13))
for field,ax in zip(['raw_GFP','raw_mCherry','primary_cluster','putative_tumor','n_genes_by_counts','pct_counts_mt'],axes.flat):
 sc.pl.umap(a,color=field,ax=ax,show=False,size=3,title=field)
fig.tight_layout();fig.savefig(out/'reporter_QC_panels.png',dpi=140);plt.close(fig)
sc.pl.dotplot(a,panels,groupby='primary_cluster',standard_scale='var',show=False)
plt.savefig(out/'major_type_dotplot.png',dpi=150,bbox_inches='tight');plt.close('all')
sc.pl.dotplot(a,['GFP','mCherry'],groupby='major_type',show=False)
plt.savefig(out/'reporter_dotplot.png',dpi=150,bbox_inches='tight');plt.close('all')
sc.tl.rank_genes_groups(a,'primary_cluster',method='wilcoxon',n_genes=50,use_raw=False)
sc.get.rank_genes_groups_df(a,group=None).to_csv(out/'cluster_top50_markers.csv',index=False)
for name in methods:
 a.obsm['X_'+name]=latents[name];a.obsm['X_umap_'+name]=coords[name]
# Restore the selected graph, since the benchmarking loop last computed scVI-zinb.
a.obsm['comparison_latent']=latents[chosen]
sc.pp.neighbors(a,use_rep='comparison_latent',n_neighbors=15,random_state=19)
a.obs.to_csv(out/'cell_annotations.csv'); a.write(out/'major_types_full.h5ad')
a[a.obs.putative_tumor=='Reporter-positive putative tumor'].write(out/'putative_tumor.h5ad')
(out/'COMPLETED.json').write_text(json.dumps(selection,indent=2))
print(metrics.to_string());print(composition.to_string());print(selection)

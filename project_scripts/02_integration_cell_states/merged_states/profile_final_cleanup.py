"""Full shared-gene descriptive profiles with within-sample enrichment support."""
import sys,json
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
from scipy import sparse
out=Path(sys.argv[1]);states=out/'states_seed42';labelpath=Path(sys.argv[2]);column=sys.argv[3];name=sys.argv[4]
dest=Path(name);dest.mkdir(parents=True,exist_ok=True)
assert not (dest/'COMPLETE.json').exists()
a=sc.read_h5ad(out/'shared_counts.h5ad');lab=pd.read_csv(labelpath,index_col=0,dtype=str)
assert lab.index.isin(a.obs_names).all();a=a[lab.index].copy();assert a.obs_names.equals(lab.index);labels=lab[column].values
meta=a.obs.copy();meta['cluster']=labels
X=a.X.tocsr();N=X.multiply(10000/meta.total_counts.values[:,None]).tocsr()
genes=a.var.gene_symbol.values
sample_data={}
for sample in meta.sample_id.unique():
 smask=meta.sample_id.values==sample
 sample_data[sample]=(smask,np.asarray(X[smask].sum(0)).ravel(),meta.total_counts.values[smask].sum(),meta.loc[smask,'model'].iloc[0])
norm_total=np.asarray(N.sum(0)).ravel();detect_total=np.asarray((X>0).sum(0)).ravel()
summary=[];alltop=[]
for cluster in sorted(np.unique(labels),key=lambda x:(len(x),x)):
 mask=labels==cluster;n=int(mask.sum());other=~mask
 norm_cluster=np.asarray(N[mask].sum(0)).ravel();detect_cluster=np.asarray((X[mask]>0).sum(0)).ravel()
 mean=norm_cluster/n;background=(norm_total-norm_cluster)/other.sum()
 pct=detect_cluster/n;pcto=(detect_total-detect_cluster)/other.sum()
 lfc=np.log2((mean+.1)/(background+.1))
 eligible=[];sample_lfcs=[];sample_models=[]
 for sample,(smask,sample_sum,sample_lib,sample_model) in sample_data.items():
  inside=smask&mask;outside=smask&other
  if inside.sum()<20 or outside.sum()<20:continue
  # Aggregate raw counts and total library sizes: sample is biological unit.
  inside_sum=np.asarray(X[inside].sum(0)).ravel();inside_lib=meta.total_counts.values[inside].sum()
  ins=inside_sum/inside_lib*1e6
  outs=(sample_sum-inside_sum)/(sample_lib-inside_lib)*1e6
  sample_lfcs.append(np.log2((ins+1)/(outs+1)));eligible.append(sample);sample_models.append(sample_model)
 samples=np.vstack(sample_lfcs) if sample_lfcs else np.empty((0,a.n_vars))
 df=pd.DataFrame({'gene_id':a.var_names,'gene':genes,'mean_norm':mean,'mean_rest':background,'pct_in':pct*100,'pct_out':pcto*100,'log2fc_mean':lfc,'eligible_samples':len(samples)})
 df['sample_median_log2fc']=np.median(samples,axis=0) if len(samples) else np.nan
 df['sample_positive_fraction']=(samples>0).mean(0) if len(samples) else np.nan
 for model in ['Npp53','Npp53Rb1']:
  subset=samples[np.array(sample_models)==model] if len(samples) else samples
  df[model+'_eligible_samples']=len(subset)
  df[model+'_median_log2fc']=np.median(subset,axis=0) if len(subset) else np.nan
 # No cell-level p-values or spurious population inference from pooled cells.
 df['marker_rank_score']=np.maximum(lfc,0)*np.maximum(pct-pcto,0)
 df=df.sort_values('marker_rank_score',ascending=False)
 df.to_csv(dest/f'cluster_{cluster}_allgenes.csv',index=False)
 endogenous=~df.gene.isin(['GFP','mCherry','Cas9','Cre','Cas9-full','Cre-full'])
 top=df[(df.pct_in>=10)&(df.log2fc_mean>=.5)&endogenous].head(30).copy();top.insert(0,'cluster',cluster);alltop.append(top)
 sub=meta.loc[mask];scount=sub.sample_id.astype(str).value_counts();rcount=sub.run_id.astype(str).value_counts()
 summary.append({'cluster':cluster,'n_cells':n,'n_samples':len(scount),'n_runs':len(rcount),'dominant_sample':scount.index[0],'dominant_sample_pct':100*scount.iloc[0]/n,'dominant_run_pct':100*rcount.iloc[0]/n,'median_genes':sub.n_genes_by_counts.median(),'median_counts':sub.total_counts.median(),'median_mt_pct':sub.pct_counts_mt.median(),'GFP_positive_pct':100*(sub.raw_GFP>0).mean(),'mCherry_positive_pct':100*(sub.raw_mCherry>0).mean(),'either_reporter_pct':100*((sub.raw_GFP>0)|(sub.raw_mCherry>0)).mean(),'eligible_samples':len(eligible),'top_enriched_genes':', '.join(top.gene.head(12))})
 print('Profiled',cluster,n,'cells;',len(eligible),'eligible samples',flush=True)
pd.DataFrame(summary).to_csv(dest/'cluster_summary.csv',index=False)
pd.concat(alltop,ignore_index=True).to_csv(dest/'top_markers.csv',index=False)
for col in ['model','treatment','run_id','sample_id','original_state']:
 pd.crosstab(meta.cluster,meta[col]).to_csv(dest/f'cluster_by_{col}.csv')
(dest/'COMPLETE.json').write_text(json.dumps({'cells':a.n_obs,'genes':a.n_vars,'clusters':len(summary),'label_file':str(labelpath),'column':column,'note':'Descriptive enrichment; sample pseudobulk support requires >=20 cells inside and outside cluster per sample. No causal or independent confirmatory tests.'},indent=2))

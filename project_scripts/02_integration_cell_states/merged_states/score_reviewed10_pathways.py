"""Current-cohort AUCell with sample-paired inference; no cell pseudoreplication."""
from pathlib import Path
import json,time,os,hashlib,shutil
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import binom
from numpy.lib.format import open_memmap
p=Path(__file__).resolve().parents[1];o=p/'analysis/merged_20260923_seed42'
parent=o/'restored_old_c18_seed42/recluster_r1_seed42/approved_merges/reviewed10'
d=parent/'pathway_statistics';d.mkdir(exist_ok=True)
assert not (d/'SCORING_COMPLETE.json').exists()
old=o/'filtered_no_c17_c18_seed42/resolution_1p0/final_cleanup/process_annotation'
start=time.time()
def log(s):
 line=f'{time.time()-start:.1f}s {s}';print(line,flush=True)
 with (d/'run.log').open('a') as h:h.write(line+'\n')
meta=pd.read_csv(parent/'final_metadata_PRIVATE.csv',index_col=0)
assert len(meta)==55133 and meta.final_cluster.nunique()==10
assert meta.loc[meta.final_cluster.eq('c9'),'putative_state'].unique().tolist()==['c9_Mes-OxidativeMetabolism']
meta.to_csv(d/'metadata_PRIVATE.csv')
for name in ['scored_signatures.csv','genesets_mouse_scored.json','signature_coverage.csv','ranking_genes.csv']:
 shutil.copy2(old/name,d/name)
sets=json.loads((d/'genesets_mouse_scored.json').read_text());sig=pd.read_csv(d/'scored_signatures.csv');names=sig.signature.tolist()
a=ad.read_h5ad(o/'shared_counts.h5ad');a=a[meta.index,a.var.gene_id.astype(str).str.startswith('ENSMUSG')].copy()
genes=pd.Index(a.var.gene_symbol);assert genes.tolist()==pd.read_csv(d/'ranking_genes.csv').gene.tolist()
x=a.X.tocsr();x.eliminate_zeros();del a
cut=int(np.ceil(.05*len(genes)));assert cut==1683
ir=[];jc=[]
for j,name in enumerate(names):
 ids=genes.get_indexer(sets[name]);assert (ids>=0).all();ir.extend(ids);jc.extend([j]*len(ids))
membership=sparse.csr_matrix((np.ones(len(ir)),(ir,jc)),shape=(len(genes),len(names)))
sizes=np.array([len(sets[n]) for n in names]);k=np.minimum(sizes,cut-1);den=k*cut-k*(k+1)/2
clusters=[f'c{i}' for i in range(1,11)];state_names={c:meta.loc[meta.final_cluster.eq(c),'putative_state'].iloc[0] for c in clusters}
record=dict(status='RUNNING',pid=os.getpid(),cells=len(meta),genes=len(genes),signatures=len(names),seeds=[42],aucMaxRank=cut,min_cells_per_arm=20,min_samples_for_test=6,test='Two-sided exact paired sign test of sample mean AUCell differences; zeros abs(delta)<=1e-8 omitted',multiple_testing='BH across all 10 states x all 1619 pathways for seed 42, including untestable p=1',comparison='Each state versus all other current tumor states within the same sample_id',inference_unit='sample_id; independence of different sample IDs assumed, animal identifiers unavailable',exploratory=True,selection_before_results='q<0.05 for seed 42, positive median delta>=0.005 and >=70% positive samples for seed 42; per-collection cap and Jaccard redundancy filter',input_sha256={str(z.relative_to(p)):hashlib.file_digest(z.open('rb'),'sha256').hexdigest() for z in [parent/'final_metadata_PRIVATE.csv',old/'genesets_mouse_scored.json']})
(d/'STARTED.json').write_text(json.dumps(record,indent=2));log(str(record))
def bh(v):
 order=np.argsort(v);ranked=np.asarray(v)[order];q=np.minimum.accumulate((ranked*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];out=np.empty(len(v));out[order]=np.minimum(q,1);return out
for seed in [42]:
 output=open_memmap(d/f'aucell_seed{seed}_PRIVATE.npy',mode='w+',dtype='float32',shape=(len(meta),len(names)))
 rng=np.random.default_rng(seed)
 for lo in range(0,len(meta),512):
  hi=min(lo+512,len(meta));rr=[];cc=[];wv=[]
  for i in range(lo,hi):
   st,en=x.indptr[i:i+2];ids=x.indices[st:en];v=x.data[st:en];ix=np.lexsort((rng.random(len(ids)),-v.astype(float)));ranks=np.arange(1,len(ids)+1)
   if seed==42 and i<24:
    if i==0:verification=np.full((len(genes),24),np.nan)
    verification[ids[ix],i]=ranks
   use=ranks<cut;rr.extend([i-lo]*int(use.sum()));cc.extend(ids[ix][use]);wv.extend(cut-ranks[use])
  w=sparse.csr_matrix((wv,(rr,cc)),shape=(hi-lo,len(genes)))
  output[lo:hi]=(w@membership).toarray()/den
  if lo%10240==0:log(f'seed {seed}: {hi}/{len(meta)}')
 output.flush();assert np.isfinite(output).all() and output.min()>=0 and output.max()<=1.000001
 rows=[];paired=[];means=[];ns=[]
 for c in clusters:
  ins=meta.final_cluster.eq(c).to_numpy();ds=[];sm=[]
  for sample in meta.sample_id.unique():
   samplemask=meta.sample_id.eq(sample).to_numpy();ii=ins&samplemask;oo=(~ins)&samplemask
   if ii.sum()<20 or oo.sum()<20:continue
   aa=output[ii].mean(0).astype(float);bb=output[oo].mean(0).astype(float);delta=aa-bb;ds.append(delta);sm.append(aa)
   model=meta.loc[ii,'model'].iloc[0];treatment=meta.loc[ii,'treatment'].iloc[0]
   for j,n in enumerate(names):paired.append(dict(seed=seed,cluster=c,sample_id=sample,model=model,treatment=treatment,signature=n,n_state=int(ii.sum()),n_other=int(oo.sum()),state_mean=aa[j],other_mean=bb[j],delta=delta[j]))
  ds=np.array(ds);n=len(ds);ns.append(n);means.append(np.mean(sm,axis=0));pos=(ds>1e-8).sum(0);neg=(ds< -1e-8).sum(0);nonzero=pos+neg
  pv=np.minimum(1,2*binom.cdf(np.minimum(pos,neg),nonzero,.5));pv[(nonzero<6)|(n<6)]=1
  for j,name in enumerate(names):rows.append(dict(seed=seed,cluster=c,state=state_names[c],signature=name,collection=sig.collection.iloc[j],eligible_samples=n,nonzero_samples=int(nonzero[j]),positive_samples=int(pos[j]),positive_fraction=float(pos[j]/n),sample_balanced_auc=float(means[-1][j]),median_delta=float(np.median(ds[:,j])),mean_delta=float(ds[:,j].mean()),p_value=float(pv[j]),testable=bool(n>=6 and nonzero[j]>=6)))
 results=pd.DataFrame(rows);results['q_value']=bh(results.p_value.to_numpy());results.to_csv(d/f'statistics_seed{seed}.csv',index=False)
 pd.DataFrame(paired).to_csv(d/f'sample_paired_scores_seed{seed}.csv',index=False)
 pd.DataFrame(means,index=clusters,columns=names).to_csv(d/f'sample_balanced_auc_seed{seed}.csv')
 if seed==42:
  pd.DataFrame(verification,index=genes,columns=meta.index[:24]).to_csv(d/'verification_rankings_PRIVATE.csv.gz',compression='gzip')
  pd.DataFrame(np.asarray(output)[:24].T,index=names,columns=meta.index[:24]).to_csv(d/'verification_python_auc_PRIVATE.csv')
 log(f'Finished seed {seed}; eligible samples={dict(zip(clusters,ns))}')
record.update(status='PASS',elapsed_seconds=time.time()-start,full_cohort_rescored=True)
(d/'SCORING_COMPLETE.json').write_text(json.dumps(record,indent=2));log('SCORING COMPLETE')

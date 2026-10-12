"""Descriptive full-gene enrichment with peer and within-sample comparisons.

No per-cell significance test or automatic biological state assignment.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='4'
from pathlib import Path
import json,time,hashlib
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

project=Path(__file__).resolve().parents[1]
out=project/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42'
dest=out/'cluster_state_inference'
dest.mkdir(exist_ok=True)
assert not (dest/'ENRICHMENT_COMPLETE.json').exists()
started=time.time()
meta=pd.read_csv(out/'metadata_PRIVATE.csv',index_col=0)
lab=pd.read_csv(out/'resolution_assignments_PRIVATE.csv',index_col=0)['leiden_r1.1_seed42'].astype(int)
assert meta.index.equals(lab.index) and len(meta)==42806 and lab.nunique()==19
source=project/'analysis/min400_20260921/tumor_umap_20260921/tumor_full_counts.h5ad'
print('Loading full counts',flush=True)
a=ad.read_h5ad(source)
assert meta.index.isin(a.obs_names).all() and a.var_names.is_unique
x=a.X[a.obs_names.get_indexer(meta.index)].tocsr().astype(np.float64)
genes=a.var_names.copy();del a
assert x.shape==(42806,25310) and np.all(x.data>=0) and np.allclose(x.data,np.rint(x.data))
assert np.all(np.asarray(x.sum(axis=1)).ravel()<=meta.total_counts.to_numpy()+1e-6)
norm=sparse.diags(1e4/meta.total_counts.to_numpy()).dot(x).tocsr()
det=x.copy();det.data[:]=1
log=norm.copy();log.data=np.log1p(log.data)
ind=sparse.csr_matrix((np.ones(len(meta)),(lab.to_numpy(),np.arange(len(meta)))),shape=(19,len(meta)))
n=np.asarray(ind.sum(axis=1)).ravel()
avg=np.asarray((sparse.diags(1/n)@ind@norm).toarray())
frac=np.asarray((sparse.diags(1/n)@ind@det).toarray())
meanlog=np.asarray((sparse.diags(1/n)@ind@log).toarray())
pd.DataFrame(avg,index=range(19),columns=genes).to_csv(dest/'cluster_mean_normalized.csv.gz')
pd.DataFrame(frac,index=range(19),columns=genes).to_csv(dest/'cluster_detection_fraction.csv.gz')
pd.DataFrame(meanlog,index=range(19),columns=genes).to_csv(dest/'cluster_mean_log1p.csv.gz')

# Broad descriptive compartments from prior marker review, not new final labels.
compartments={'epithelial':[1,2,9,12,13,15,17,18],
              'mesenchymal_or_transitional':[0,3,4,5,6,7,8,10,11],
              'mixed_unresolved':[14], 'cycling_unresolved':[16]}
groupmap={c:k for k,v in compartments.items() for c in v}
assert sorted(groupmap)==list(range(19))
samples=sorted(meta.sample_id.unique())
keys=pd.MultiIndex.from_product([samples,range(19)])
codes=keys.get_indexer(pd.MultiIndex.from_arrays([meta.sample_id,lab]))
sampleind=sparse.csr_matrix((np.ones(len(meta)),(codes,np.arange(len(meta)))),shape=(len(keys),len(meta)))
ns=np.asarray(sampleind.sum(axis=1)).ravel().reshape(len(samples),19)
sample_sum=(sampleind@norm).toarray().reshape(len(samples),19,len(genes))
pd.DataFrame(ns,index=samples,columns=range(19)).to_csv(dest/'sample_cluster_counts.csv')
exclude=genes.str.match(r'^(mt-|Rp[sl])') | genes.isin(['GFP','mCherry','Cas9-full','Cas9','Cre','WPRE'])
rows=[];tops=[];nearestrows=[]
for cl in range(19):
    same=[k for k in compartments[groupmap[cl]] if k!=cl]
    if not same:same=[k for k in range(19) if k!=cl]
    # Closest expression centroid within its broad compartment: diagnostic only.
    correlations=np.corrcoef(meanlog[:,~exclude])
    near=max(same,key=lambda k:correlations[cl,k])
    nearestrows.append(dict(cluster=cl,compartment=groupmap[cl],nearest_cluster=near,correlation=correlations[cl,near]))
    for comparison,others in [('all_other',[k for k in range(19) if k!=cl]),('peer',same),('nearest',[near])]:
        nr=n[others].sum()
        otheravg=(avg[others]*n[others,None]).sum(axis=0)/nr
        otherfrac=(frac[others]*n[others,None]).sum(axis=0)/nr
        otherlog=(meanlog[others]*n[others,None]).sum(axis=0)/nr
        lfc=np.log2((avg[cl]+.01)/(otheravg+.01))
        eligible=(ns[:,cl]>=20)&(ns[:,others].sum(axis=1)>=20)
        denom=ns[:,others].sum(axis=1)
        s_own=sample_sum[:,cl,:][eligible]/ns[eligible,cl,None]
        s_peer=sample_sum[:,others,:].sum(axis=1)[eligible]/denom[eligible,None]
        within=np.log2((s_own+.01)/(s_peer+.01))
        med=np.median(within,axis=0) if eligible.any() else np.full(len(genes),np.nan)
        consistent=(within>0).sum(axis=0) if eligible.any() else np.zeros(len(genes),dtype=int)
        tab=pd.DataFrame(dict(cluster=cl,comparison=comparison,gene=genes,n_cells=int(n[cl]),
            mean_normalized=avg[cl],mean_other=otheravg,log2FC=lfc,
            pct_in=100*frac[cl],pct_other=100*otherfrac,
            delta_pct=100*(frac[cl]-otherfrac),mean_log_difference=meanlog[cl]-otherlog,
            eligible_samples=int(eligible.sum()),samples_positive=consistent,median_within_sample_log2FC=med))
        tab['specific_enrichment']=(~exclude)&(frac[cl]>=.15)&(lfc>=.5)&(frac[cl]-otherfrac>=.1)
        # Sort only qualified genes; no p-value/false-discovery claims.
        selected=tab.loc[tab.specific_enrichment].sort_values(['mean_log_difference','log2FC'],ascending=False).head(40)
        tops.append(selected)
        tab.to_csv(dest/f'cl{cl}_{comparison}_all_genes.csv.gz',index=False)
    print('Profiled cluster',cl,flush=True)
pd.concat(tops,ignore_index=True).to_csv(dest/'top_specific_genes.csv',index=False)
pd.DataFrame(nearestrows).to_csv(dest/'nearest_peer_clusters.csv',index=False)
manifest=dict(cells=42806,clusters=19,resolution=1.1,seed=42,genes=25310,elapsed_seconds=time.time()-started,
  method='Descriptive log2 ratios of mean counts per10000 original full-library UMIs; pseudocount .01; detection contrasts; no p-values',
  specific_filter='pct_in>=15%, log2FC>=0.5, delta_pct>=10 percentage points; exclude mt/ribosomal and tracing transgenes from ranking',
  sample_support='Within-sample mean normalized expression ratio; at least20 cells on each side; report sign and median, not significance',
  compartments=compartments,
  input_sha256={f:hashlib.sha256((out/f).read_bytes()).hexdigest() for f in ['metadata_PRIVATE.csv','resolution_assignments_PRIVATE.csv','scvi_latent_PRIVATE.csv']})
(dest/'ENRICHMENT_COMPLETE.json').write_text(json.dumps(manifest,indent=2))
print('PASS: descriptive full-gene enrichment complete; no clusters or states changed.',flush=True)

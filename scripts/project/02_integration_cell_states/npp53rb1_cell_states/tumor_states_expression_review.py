"""Post-hoc expression/QC evidence for all candidate resolutions, not selection inputs."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
parser=argparse.ArgumentParser()
parser.add_argument('--source',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--reference',type=Path,required=True)
parser.add_argument('--seed',type=int,default=19)
parser.add_argument('--skip-islands',action='store_true',help='Do not reuse historical geometric island comparisons after exclusions')
args=parser.parse_args();out=args.output;dest=out/'expression_review';dest.mkdir(exist_ok=True)
assert not (dest/'COMPLETE.json').exists()
meta=pd.read_csv(out/'metadata_PRIVATE.csv',index_col=0)
assign=pd.read_csv(out/'resolution_assignments_PRIVATE.csv',index_col=0,dtype=str).loc[meta.index]
a=ad.read_h5ad(args.source/'tumor_full_counts.h5ad')
expected=json.loads((out/'PREPARED.json').read_text())['cells']
assert meta.index.is_unique and a.obs_names.is_unique and meta.index.isin(a.obs_names).all()
assert len(meta)==expected
if not a.obs_names.equals(meta.index):
    assert args.skip_islands,'Filtered cohorts require explicitly skipping historical island comparisons'
    a=a[meta.index].copy()
assert a.obs_names.equals(meta.index)
x=a.X.tocsr().astype(np.float32);genes=a.var_names
totals=meta.total_counts.to_numpy(dtype=float)
retained=np.asarray(x.sum(axis=1)).ravel()
assert (retained<=totals).all()
meta['pct_counts_ribo']=100*np.asarray(x[:,genes.str.match('^Rp[sl]')].sum(axis=1)).ravel()/totals
norm=sparse.diags(1e4/totals).dot(x).tocsr().astype(np.float32);norm.data=np.log1p(norm.data)
detect=x.copy();detect.data[:]=1
panels={
 'Luminal':['Krt8','Krt18','Krt19','Epcam','Cdh1','Foxa1'],
 'Basal':['Krt5','Krt14','Trp63','Krt17','Krt15'],
 'AR_axis':['Ar','Nkx3-1','Fkbp5','Tmprss2','Pbsn','Msmb'],
 'Mesenchymal':['Vim','Zeb1','Zeb2','Fn1','Col1a1','Col3a1','Dcn','Lum','Pdgfra'],
 'myCAF_smooth_muscle':['Acta2','Tagln','Myh11','Rgs5','Pdgfrb','Lrrc15','Postn','Fap','Aspn'],
 'Neuroendocrine':['Chga','Chgb','Syp','Ascl1','Insm1','Ncam1','Neurod1','Dll3','Ddc'],
 'Ionocyte':['Foxi1','Atp6v1g3','Aqp6','Clcnkb','Atp6v0d2'],
 'Immune':['Ptprc','Lyz2','Tyrobp','Adgre1','Cd3d','Cd79a','Nkg7'],
 'Endothelial':['Pecam1','Cdh5','Kdr','Emcn','Vwf'],
 'Erythroid_platelet':['Hba-a1','Hbb-bs','Alas2','Pf4','Ppbp'],
 'Cycling':['Mki67','Top2a','Pclaf','Birc5','Ccnb1','Cdk1'],
 'Stress':['Fos','Jun','Atf3','Ddit3','Hspa1a','Hspa1b'],
 'Reporters':['GFP','mCherry']}
reference=pd.read_csv(args.reference)
coverage_sets={**{'Npp53_'+str(name):tab.gene.tolist() for name,tab in reference.groupby('cluster',sort=False)},**panels}
pd.DataFrame([dict(program=name,gene=g,present_in_counts=g in genes)
              for name,gs in coverage_sets.items() for g in gs]).to_csv(dest/'marker_coverage.csv',index=False)
ref={name:[g for g in tab.gene if g in genes] for name,tab in reference.groupby('cluster',sort=False)}
allsets={**{'Npp53_'+k:v for k,v in ref.items()},**panels}
selected=list(dict.fromkeys(g for gs in allsets.values() for g in gs if g in genes))
sm=norm[:,genes.get_indexer(selected)].toarray()
mu=sm.mean(axis=0);sd=sm.std(axis=0);sd[sd==0]=1
z=(sm-mu)/sd
lookup={g:i for i,g in enumerate(selected)}
scores=pd.DataFrame({name:z[:,[lookup[g] for g in gs if g in lookup]].mean(axis=1)
                     for name,gs in allsets.items() if len([g for g in gs if g in lookup])>=2},index=meta.index)
scores.to_csv(dest/'signature_scores_PRIVATE.csv')
(dest/'signature_dictionary.json').write_text(json.dumps(allsets,indent=2)+'\n')
def profile(group,prefix):
    group=pd.Series(group,index=meta.index).astype(str)
    cats=sorted(group.unique(),key=lambda s:(not s.isdigit(),int(s) if s.isdigit() else s))
    index=pd.Categorical(group,categories=cats).codes
    indicator=sparse.csr_matrix((np.ones(len(meta)),(index,np.arange(len(meta)))),shape=(len(cats),len(meta)))
    n=np.asarray(indicator.sum(axis=1)).ravel()
    means=(sparse.diags(1/n)@indicator@norm).toarray()
    fractions=(sparse.diags(1/n)@indicator@detect).toarray()
    grand=np.asarray(norm.sum(axis=0)).ravel()
    denom=np.maximum(len(meta)-n,1)
    other=(grand[None,:]-means*n[:,None])/denom[:,None]
    enriched=[];canonical=[];quality=[]
    for i,cl in enumerate(cats):
        valid=(fractions[i]>=.15)&~genes.str.match(r'^(mt-|Rp[sl])')&~genes.isin(['GFP','mCherry'])
        best=np.argsort(np.where(valid,means[i]-other[i],-np.inf))[-30:][::-1]
        for j in best:enriched.append(dict(cluster=cl,n=int(n[i]),gene=genes[j],mean_log1p_10k=float(means[i,j]),
                                            fraction_positive=float(fractions[i,j]),mean_log_difference=float(means[i,j]-other[i,j])))
        for family,gs in panels.items():
            for g in gs:
                if g not in genes:continue
                j=genes.get_loc(g)
                canonical.append(dict(cluster=cl,program=family,gene=g,mean_log1p_10k=float(means[i,j]),fraction_positive=float(fractions[i,j])))
        sub=meta.loc[group==cl]
        quality.append(dict(cluster=cl,n=len(sub),samples=sub.sample_id.nunique(),runs=sub.run_id.nunique(),
                            largest_sample_fraction=float(sub.sample_id.value_counts().max()/len(sub)),
                            median_genes=float(sub.n_genes_by_counts.median()),median_UMI=float(sub.total_counts.median()),
                            median_mt=float(sub.pct_counts_mt.median()),median_ribo=float(sub.pct_counts_ribo.median()),
                            reporter_positive=int(((sub.raw_GFP>0)|(sub.raw_mCherry>0)).sum()),
                            reporter_ge2=int(((sub.raw_GFP>=2)|(sub.raw_mCherry>=2)).sum())))
    pd.DataFrame(enriched).to_csv(dest/f'{prefix}_descriptive_top30.csv',index=False)
    pd.DataFrame(canonical).to_csv(dest/f'{prefix}_canonical_markers.csv',index=False)
    pd.DataFrame(quality).to_csv(dest/f'{prefix}_QC.csv',index=False)
    scores.groupby(group).mean().to_csv(dest/f'{prefix}_signature_means.csv')
    for field in ['sample_id','run_id','treatment','parent_harmony_cluster','review_region']:
        pd.crosstab(group,meta[field]).to_csv(dest/f'{prefix}_{field}.csv')
for i in range(1,13):
    tag=f'r{i/10:.1f}';profile(assign[f'leiden_{tag}_seed{args.seed}'],tag);print('Profiled',tag,flush=True)
if args.skip_islands:
    (dest/'COMPLETE.json').write_text(json.dumps(dict(cells=len(meta),genes=len(genes),resolutions=12,seed=args.seed,
        normalization='log1p corrected counts per 10,000 original full-library UMIs',
        signature_scoring='Mean gene-wise z-scored log expression within the retained cohort; Npp53 top25 reference plus safeguard panels; not exact AddModuleScore reproduction',
        rankings='Descriptive mean-log-expression difference; no DE hypothesis test or biological replication claim',
        historical_island_comparisons=False,selection_input=False),indent=2)+'\n')
    print('Filtered-cohort full-gene expression review COMPLETE',flush=True)
    raise SystemExit(0)
profile(meta.review_region,'islands')
# Matched island summaries: same sample × parent cluster, excluding both islands.
strata=meta.sample_id.astype(str)+'|'+meta.parent_harmony_cluster.astype(str)
matchedrows=[];matchedgenes=[]
for region in ['Upper island','Far-right island']:
    flag=meta.review_region==region;weights=np.zeros(len(meta))
    for key,count in strata[flag].value_counts().items():
        control=(strata==key)&(meta.review_region=='Background');assert control.sum()>0
        weights[control]=count/flag.sum()/control.sum()
    assert np.isclose(weights.sum(),1)
    for metric in ['n_genes_by_counts','total_counts','pct_counts_mt','pct_counts_ribo']:
        v=meta[metric].to_numpy();order=np.argsort(v);wm=v[order[np.searchsorted(np.cumsum(weights[order]),.5)]]
        matchedrows.append(dict(region=region,metric=metric,island_median=float(np.median(v[flag])),matched_control_median=float(wm)))
    avg=np.asarray(norm[flag].mean(axis=0)).ravel();ctrl=np.asarray(norm.T@weights).ravel()
    frac=np.asarray(detect[flag].mean(axis=0)).ravel();cf=np.asarray(detect.T@weights).ravel()
    valid=(frac>=.1)&~genes.str.match(r'^(mt-|Rp[sl])')&~genes.isin(['GFP','mCherry'])
    for j in np.argsort(np.where(valid,avg-ctrl,-np.inf))[-40:][::-1]:
        matchedgenes.append(dict(region=region,gene=genes[j],mean_log1p_10k=float(avg[j]),matched_mean=float(ctrl[j]),
                                 fraction_positive=float(frac[j]),matched_fraction=float(cf[j]),difference=float(avg[j]-ctrl[j])))
pd.DataFrame(matchedrows).to_csv(dest/'island_matched_QC.csv',index=False)
pd.DataFrame(matchedgenes).to_csv(dest/'island_matched_descriptive_genes.csv',index=False)
(dest/'COMPLETE.json').write_text(json.dumps(dict(cells=len(meta),genes=len(genes),resolutions=12,
    normalization='log1p corrected counts per 10,000 original full-library UMIs',
    signature_scoring='Mean gene-wise z-scored log expression; Npp53 top25 reference plus safeguard panels; not exact AddModuleScore reproduction',
    rankings='Descriptive mean-log-expression difference; no DE hypothesis test or biological replication claim',
    selection_input=False),indent=2)+'\n')
print('Expression review COMPLETE',flush=True)

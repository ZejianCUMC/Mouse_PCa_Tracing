"""Descriptive candidate QC/marker evidence, separate from CHOIR/AMI selection.

Uses the saved 3,000-HVG count matrix; missing genes are unavailable, not zero.
No gene p-values, state names, model fitting or changes to candidate ranking.
"""
import gzip,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

project=Path(__file__).resolve().parents[1]
out=project/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42'
dest=out/'candidate_diagnostics'
dest.mkdir(exist_ok=True)
assert not (dest/'COMPLETE.json').exists()
meta=pd.read_csv(out/'metadata_PRIVATE.csv',index_col=0)
assign=pd.read_csv(out/'resolution_assignments_PRIVATE.csv',index_col=0,dtype=str)
screen=pd.read_csv(out/'screening_labels_PRIVATE.csv',index_col=0,dtype=str)
assert meta.index.equals(assign.index) and meta.index.equals(screen.index)
assert len(meta)==42806 and (out/'barcodes_PRIVATE.txt').read_text().splitlines()==meta.index.tolist()
genes=pd.Index((out/'features.txt').read_text().splitlines())
x=mmread(gzip.open(out/'hvg_counts_PRIVATE.mtx.gz','rb')).T.tocsr().astype(np.float32)
assert x.shape==(42806,3000) and genes.is_unique
norm=sparse.diags(1e4/meta.total_counts.to_numpy())@x
norm=norm.tocsr();norm.data=np.log1p(norm.data)
det=x.copy();det.data[:]=1
panels={'Basal':['Krt5','Krt14','Trp63','Krt17','Krt15'],
        'Luminal':['Krt8','Krt18','Krt19','Epcam','Foxa1','Ar','Nkx3-1'],
        'Mesenchymal':['Vim','Col1a1','Col3a1','Acta2','Tagln','Dcn','Fmod','Spp1','Zeb1','Fn1'],
        'Epithelial modifiers':['Lcn2','Clu','Ly6a','Tacstd2','Cd44'],
        'Cycling':['Mki67','Top2a','Pclaf','Birc5','Ccnb1','Cdk1'],
        'Stress':['Fos','Jun','Atf3','Ddit3','Hspa1a','Hspa1b'],
        'Host controls':['Ptprc','Lyz2','Tyrobp','Pecam1','Cdh5'],
        'Neuroendocrine controls':['Chga','Chgb','Syp','Ascl1','Insm1']}
coverage=[dict(program=family,gene=g,present_in_HVGs=g in genes) for family,gs in panels.items() for g in gs]
pd.DataFrame(coverage).to_csv(dest/'marker_coverage.csv',index=False)
quality=[];markers=[];enriched=[];composition=[]
total=np.asarray(norm.sum(axis=0)).ravel()
for r in np.round(np.arange(.1,1.21,.1),1):
    group=assign[f'leiden_r{r:.1f}_seed42']
    categories=sorted(group.unique(),key=int)
    codes=pd.Categorical(group,categories=categories).codes
    indicator=sparse.csr_matrix((np.ones(len(meta)),(codes,np.arange(len(meta)))),shape=(len(categories),len(meta)))
    n=np.asarray(indicator.sum(axis=1)).ravel()
    means=(sparse.diags(1/n)@indicator@norm).toarray()
    fractions=(sparse.diags(1/n)@indicator@det).toarray()
    other=(total[None,:]-means*n[:,None])/(len(meta)-n[:,None])
    for i,cl in enumerate(categories):
        mask=group.eq(cl);sub=meta.loc[mask]
        quality.append(dict(resolution=float(r),cluster=cl,cells=len(sub),samples=sub.sample_id.nunique(),
            runs=sub.run_id.nunique(),largest_sample_fraction=sub.sample_id.value_counts().iloc[0]/len(sub),
            dominant_sample=sub.sample_id.value_counts().index[0],largest_run_fraction=sub.run_id.value_counts().iloc[0]/len(sub),
            median_genes=sub.n_genes_by_counts.median(),median_UMI=sub.total_counts.median(),median_mt=sub.pct_counts_mt.median(),
            reporter_positive=int(((sub.raw_GFP>0)|(sub.raw_mCherry>0)).sum()),
            historical_cl2_fraction=screen.loc[mask,'cluster'].eq('2').mean()))
        valid=(fractions[i]>=.15)&~genes.str.match(r'^(mt-|Rp[sl])')
        indices=np.flatnonzero(valid)
        best=indices[np.argsort((means[i]-other[i])[indices])[-30:][::-1]]
        for j in best:
            enriched.append(dict(resolution=float(r),cluster=cl,gene=genes[j],mean_log1p_10k=means[i,j],
                fraction_positive=fractions[i,j],mean_log_difference=means[i,j]-other[i,j]))
        for family,gs in panels.items():
            for gene in gs:
                if gene in genes:
                    j=genes.get_loc(gene)
                    markers.append(dict(resolution=float(r),cluster=cl,program=family,gene=gene,
                        mean_log1p_10k=means[i,j],fraction_positive=fractions[i,j]))
        for (sample,run,treatment),count in sub.groupby(['sample_id','run_id','treatment']).size().items():
            composition.append(dict(resolution=float(r),cluster=cl,sample_id=sample,run_id=run,treatment=treatment,cells=int(count)))
pd.DataFrame(quality).to_csv(dest/'candidate_QC.csv',index=False)
pd.DataFrame(markers).to_csv(dest/'candidate_HVG_markers.csv',index=False)
pd.DataFrame(enriched).to_csv(dest/'candidate_HVG_top30.csv',index=False)
pd.DataFrame(composition).to_csv(dest/'candidate_sample_composition.csv',index=False)
(dest/'COMPLETE.json').write_text(json.dumps(dict(cells=len(meta),resolutions=12,seed=42,
    features=3000,role='Post-hoc interpretation only; not used for AMI rank or statistical tests',
    limitation='HVG-only descriptive expression. Genes outside HVGs are unavailable, not negative. No DE p-values or treatment-effect inference.'),indent=2))
assert pd.DataFrame(quality).groupby('resolution').cells.sum().eq(42806).all()
print('PASS: twelve candidate QC/composition profiles and HVG-only marker evidence reconciled.')

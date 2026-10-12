"""Read-only-source, descriptive island audit. No filtering, fitting, clustering or DE tests."""
import argparse, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser()
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
root=args.out; dest=root/'island_review_20260922';dest.mkdir(exist_ok=True)
started=time.time()
meta=pd.read_csv(root/'metadata.csv',index_col=0)
xy=pd.read_csv(root/'UMAP_scVI_run_sample.csv',index_col=0).loc[meta.index]
assert len(meta)==43057 and meta.index.is_unique
# Gates lie in the largest empty coordinate gaps, fixed before expression review.
labels=np.full(len(meta),'Background',dtype=object)
labels[xy.UMAP_2.to_numpy()>12.36]='Upper'
labels[xy.UMAP_1.to_numpy()>8]='Far-right'
assert (labels=='Upper').sum()==210 and (labels=='Far-right').sum()==48
meta['island_review']=labels
meta['genes_per_1000_UMI']=1000*meta.n_genes_by_counts/meta.total_counts
print('Loading existing count matrix',flush=True)
a=ad.read_h5ad(root/'tumor_full_counts.h5ad')
assert a.obs_names.equals(meta.index)
counts=a.X.tocsr().astype(np.float32)
assert np.allclose(counts.data,np.round(counts.data))
retained_totals=np.asarray(counts.sum(axis=1)).ravel()
# Preparation preserved original totals, then removed genes seen in <3 cells.
# Normalize with original library sizes and explicitly audit the missing UMIs.
totals=meta.total_counts.to_numpy(dtype=float)
lost=totals-retained_totals
assert np.all(lost>=0) and np.allclose(totals,a.obs.total_counts)
count_audit={'matrix_genes':a.n_vars,'cells_with_rare_gene_UMIs_omitted':int((lost>0).sum()),
             'max_UMIs_omitted':float(lost.max()),'max_fraction_omitted':float((lost/totals).max()),
             'normalization_denominator':'Original pre-gene-filter total UMI counts from metadata'}
print(count_audit,flush=True)
genes=a.var_names
groups={
 'Basal':['Krt5','Krt14','Krt17','Trp63'],
 'Luminal':['Epcam','Krt8','Krt18','Krt19'],
 'Mesenchymal':['Vim','Zeb1','Col1a1','Col3a1','Dcn','Lum'],
 'Smooth_muscle':['Acta2','Tagln','Myh11','Rgs5','Pdgfrb'],
 'Immune':['Ptprc','Lyz2','Tyrobp','Adgre1','Cd3d','Cd79a'],
 'Endothelial':['Pecam1','Cdh5','Kdr','Emcn'],
 'Erythroid_platelet':['Hba-a1','Hbb-bs','Alas2','Pf4','Ppbp'],
 'Neuroendocrine':['Chga','Chgb','Syp','Ascl1'],
 'Stress':['Fos','Jun','Atf3','Ddit3','Hspa1a','Hspa1b'],
 'Cycling':['Mki67','Top2a'],
 'Reporters':['GFP','mCherry']}
gene_list=[g for gs in groups.values() for g in gs if g in genes]
norm=sparse.diags(1e4/np.maximum(totals,1)).dot(counts).tocsr()
norm.data=np.log1p(norm.data)
detected=counts.copy();detected.data[:]=1
rib=genes.str.match(r'^Rp[sl]')
meta['pct_counts_ribo']=100*np.asarray(counts[:,rib].sum(axis=1)).ravel()/totals
sel=genes.get_indexer(gene_list)
small=counts[:,sel].toarray()
smallnorm=norm[:,sel].toarray()
metrics=['n_genes_by_counts','total_counts','pct_counts_mt','pct_counts_ribo','genes_per_1000_UMI','raw_GFP','raw_mCherry']
strata=meta.sample_id.astype(str)+'|'+meta.parent_harmony_cluster.astype(str)
results={};quality=[];markers=[];enrichments=[];by_sample=[]
weights_by_island={}
for name in ['Upper','Far-right']:
    mask=labels==name;idx=np.flatnonzero(mask); n=int(mask.sum())
    weights=np.zeros(len(meta)); cell_percentiles={m:[] for m in metrics}
    for key,nkey in strata[mask].value_counts().items():
        control=(strata.to_numpy()==key)&(labels=='Background')
        assert control.sum()>0
        weights[control]=nkey/n/control.sum()
        for m in metrics:
            ref=np.sort(meta.loc[control,m].to_numpy())
            values=meta.loc[mask&(strata.to_numpy()==key),m].to_numpy()
            cell_percentiles[m].extend((np.searchsorted(ref,values,side='left')+np.searchsorted(ref,values,side='right'))/2/len(ref))
    assert np.isclose(weights.sum(),1)
    weights_by_island[name]=weights
    row={'n':n,'samples':meta.loc[mask,'sample_id'].value_counts().to_dict(),
         'runs':meta.loc[mask,'run_id'].value_counts().to_dict(),
         'treatments':meta.loc[mask,'treatment'].value_counts().to_dict(),
         'parent_clusters':meta.loc[mask,'parent_harmony_cluster'].astype(str).value_counts().to_dict(),
         'anchors':meta.loc[mask,'expression_anchor'].value_counts().to_dict(),
         'reporter_positive':int(((meta.raw_GFP>0)|(meta.raw_mCherry>0))[mask].sum()),
         'reporter_ge2':int(((meta.raw_GFP>=2)|(meta.raw_mCherry>=2))[mask].sum()),
         'qc':{},'matched_control_cells':int((weights>0).sum())}
    for m in metrics:
        v=meta[m].to_numpy();order=np.argsort(v);cum=np.cumsum(weights[order])
        wmedian=float(v[order[np.searchsorted(cum,.5)]])
        q=dict(island=name,metric=m,median=float(np.median(v[mask])),q25=float(np.quantile(v[mask],.25)),q75=float(np.quantile(v[mask],.75)),
               matched_control_median=wmedian,median_within_stratum_percentile=float(np.median(cell_percentiles[m])))
        quality.append(q);row['qc'][m]=q
    for sample, subset in meta.loc[mask].groupby('sample_id',observed=True):
        by_sample.append(dict(island=name,sample_id=sample,n=len(subset),run_id=subset.run_id.iloc[0],treatment=subset.treatment.iloc[0],
                              reporter_positive=int(((subset.raw_GFP>0)|(subset.raw_mCherry>0)).sum()),
                              median_genes=float(subset.n_genes_by_counts.median()),median_UMI=float(subset.total_counts.median()),median_mt=float(subset.pct_counts_mt.median())))
    avg=np.asarray(norm[mask].mean(axis=0)).ravel()
    ctrl=np.asarray(norm.T.dot(weights)).ravel()
    pct=np.asarray(detected[mask].mean(axis=0)).ravel()
    ctrlpct=np.asarray(detected.T.dot(weights)).ravel()
    for family, panel in groups.items():
        for g in panel:
            if g not in genes:continue
            j=genes.get_loc(g)
            markers.append(dict(island=name,program=family,gene=g,fraction_positive=float(pct[j]),matched_control_fraction=float(ctrlpct[j]),
                                mean_log1p_10k=float(avg[j]),matched_control_mean=float(ctrl[j]),difference=float(avg[j]-ctrl[j])))
    eligible=(pct>=.1)&~genes.str.match(r'^(mt-|Rp[sl])')&~genes.isin(['GFP','mCherry'])
    top=np.argsort(np.where(eligible,avg-ctrl,-np.inf))[-40:][::-1]
    for j in top:
        enrichments.append(dict(island=name,gene=genes[j],fraction_positive=float(pct[j]),matched_control_fraction=float(ctrlpct[j]),
                                mean_log1p_10k=float(avg[j]),matched_control_mean=float(ctrl[j]),difference=float(avg[j]-ctrl[j])))
    # Coexpression screen is descriptive, NOT a doublet classifier.
    gidx={g:i for i,g in enumerate(gene_list)}
    def n_detect(panel):return (small[:,[gidx[g] for g in panel if g in gidx]]>0).sum(axis=1)
    ep=n_detect(groups['Basal']+groups['Luminal'])>=2
    imm=n_detect(['Ptprc','Lyz2','Tyrobp','Adgre1','Cd3d','Cd79a'])>=2
    end=n_detect(groups['Endothelial'])>=2
    row['coexpression_screen']={'epithelial_ge2':int(ep[mask].sum()),'immune_ge2':int(imm[mask].sum()),'endothelial_ge2':int(end[mask].sum()),
                                'epithelial_and_immune':int((ep&imm)[mask].sum()),'epithelial_and_endothelial':int((ep&end)[mask].sum())}
    results[name]=row
pd.DataFrame(quality).to_csv(dest/'quality_summary.csv',index=False)
pd.DataFrame(markers).to_csv(dest/'marker_summary.csv',index=False)
pd.DataFrame(enrichments).to_csv(dest/'descriptive_gene_enrichment.csv',index=False)
pd.DataFrame(by_sample).to_csv(dest/'island_by_sample.csv',index=False)
meta.loc[labels!='Background'].to_csv(dest/'island_cell_metadata_PRIVATE.csv')
print('Count/marker review complete; checking original 30D neighborhoods',flush=True)
b=ad.read_h5ad(root/'tumor_umap_comparison.h5ad')
assert b.obs_names.equals(meta.index)
nnrows=[]
rng=np.random.default_rng(19)
for method in ['uncorrected','Harmony_run','scVI_run','scVI_sample','scVI_run_sample']:
    latent=np.asarray(b.obsm['X_'+method])
    assert latent.shape==(43057,30) and np.isfinite(latent).all()
    nn=NearestNeighbors(n_neighbors=16,algorithm='brute',n_jobs=4).fit(latent)
    for name in ['Upper','Far-right']:
        mask=labels==name;idx=np.flatnonzero(mask)
        dist,inds=nn.kneighbors(latent[idx])
        clean=[[(d,j) for d,j in zip(ds,js) if j!=i][:15] for i,ds,js in zip(idx,dist,inds)]
        near=np.array([[j for d,j in x] for x in clean]);d15=np.array([x[-1][0] for x in clean])
        within=(labels[near]==name).mean(axis=1)
        controlidx=rng.choice(len(meta),size=500,replace=True,p=weights_by_island[name])
        cd,ci=nn.kneighbors(latent[controlidx])
        control_d15=np.array([[d for d,j in zip(ds,js) if j!=i][:15][-1] for i,ds,js in zip(controlidx,cd,ci)])
        nnrows.append(dict(island=name,method=method,mean_same_island_fraction_k15=float(within.mean()),median_same_island_fraction_k15=float(np.median(within)),
                           median_distance_k15=float(np.median(d15)),matched_control_median_distance_k15=float(np.median(control_d15)),
                           distance_ratio=float(np.median(d15)/np.median(control_d15))))
    print('Neighbors checked:',method,flush=True)
pd.DataFrame(nnrows).to_csv(dest/'latent_neighborhood_summary.csv',index=False)
palette={'Upper':'#c03aaa','Far-right':'#148c66'}
plt.rcParams.update({'font.size':12,'axes.titlesize':14})
fig,axes=plt.subplots(1,3,figsize=(16,5.5))
for ax,method,title in zip(axes,['scVI_run_sample','Harmony_run','uncorrected'],['Current scVI run + sample','Harmony run θ1','Uncorrected PCA → UMAP']):
    pos=pd.read_csv(root/('UMAP_'+method+'.csv'),index_col=0).loc[meta.index].to_numpy()
    ax.scatter(*pos.T,s=.6,c='#d9d9d9',rasterized=True)
    for name in ['Upper','Far-right']:
        part=pos[labels==name];ax.scatter(*part.T,s=9,c=palette[name],label=f'{name} (n={len(part)})',rasterized=True)
    ax.set(xlabel='UMAP_1',ylabel='UMAP_2',title=title);ax.set_aspect('equal',adjustable='datalim')
axes[0].axhline(12.36,ls=':',color=palette['Upper'],lw=.7)
axes[0].axvline(8,ls=':',color=palette['Far-right'],lw=.7)
axes[0].legend(loc='lower left',fontsize=10,markerscale=2)
fig.suptitle('Same cells highlighted in three independently computed embeddings; orientations are not comparable',fontsize=14)
fig.tight_layout();fig.savefig(dest/'island_locations.png',dpi=180);plt.close(fig)
# Predetermined markers, retaining axes/ticks and matched-comparison denominators.
plotgenes=['Krt5','Krt14','Epcam','Krt8','Krt18','Krt19','Vim','Zeb1','Col1a1','Dcn','Acta2','Tagln','Myh11','Ptprc','Lyz2','Tyrobp','Pecam1','Cdh5','Hba-a1','Pf4','Chga','Syp','Fos','Ddit3','Mki67','GFP','mCherry']
plotgenes=[g for g in plotgenes if g in genes]
fig,axes=plt.subplots(2,1,figsize=(16,8),gridspec_kw={'height_ratios':[1,1.3]})
qdf=pd.DataFrame(quality)
for i,m in enumerate(['n_genes_by_counts','total_counts','pct_counts_mt','pct_counts_ribo']):
    for off,name in [(-.18,'Upper'),(.18,'Far-right')]:
        row=qdf[(qdf.island==name)&(qdf.metric==m)].iloc[0]
        axes[0].bar(i+off,row['median']/row.matched_control_median,width=.32,color=palette[name],label=name if i==0 else None)
axes[0].axhline(1,color='black',ls='--',lw=1)
axes[0].set(xticks=range(4),xticklabels=['Detected genes','UMI counts','Mitochondrial %','Ribosomal %'],ylabel='Island / matched-control median',title='Quality relative to sample × parent-cluster matched background')
axes[0].legend()
mdf=pd.DataFrame(markers)
yy=0
for name in ['Upper','Far-right']:
    data=mdf[mdf.island==name].set_index('gene').loc[plotgenes]
    for suffix,frac,expr in [('',data.fraction_positive,data.mean_log1p_10k),(' matched background',data.matched_control_fraction,data.matched_control_mean)]:
        pts=axes[1].scatter(range(len(plotgenes)),np.repeat(yy,len(plotgenes)),s=frac.to_numpy()*180,c=expr,cmap='Blues',vmin=0,vmax=5,edgecolors='#888888',linewidths=.3)
        yy+=1
axes[1].set(xticks=range(len(plotgenes)),xticklabels=plotgenes,yticks=range(4),yticklabels=['Upper','Upper matched background','Far-right','Far-right matched background'],ylim=(3.6,-.6),title='Predefined marker panel: area = fraction detected; colour = mean log1p(counts per 10,000)')
axes[1].tick_params(axis='x',rotation=65)
for val in [.25,.5,1]:axes[1].scatter([],[],s=180*val,c='grey',label=f'{val:.0%}')
axes[1].legend(title='Fraction detected',bbox_to_anchor=(1.01,.9),loc='upper left',fontsize=9)
fig.colorbar(pts,ax=axes[1],pad=.13,label='Mean log expression',fraction=.02)
fig.tight_layout();fig.savefig(dest/'island_evidence.png',dpi=180);plt.close(fig)
payload={'gate_basis':'Largest empty gaps in current saved scVI run+sample UMAP; upper UMAP_2 >12.36; far-right UMAP_1 >8; fixed before inspecting genes',
 'cohort_cells':43057,'count_audit':count_audit,'islands':results,'matched_control':'Non-island cells matched by sample_id × parent_harmony_cluster; strata weighted by island composition',
 'methods':'Descriptive only; no hypothesis tests, model refitting, new clustering, cell removal or malignancy confirmation. Historical doublet filtering retained; no new doublet classifier.',
 'elapsed_seconds':round(time.time()-started,2)}
(dest/'review_summary.json').write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps(payload,indent=2),flush=True)

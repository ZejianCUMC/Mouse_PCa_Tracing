"""Display the AMI-leading partition, not an automatic final state assignment."""
import json, os, tempfile, hashlib, shutil, subprocess, time
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='rb1-candidate-review-'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import anndata as ad
from scipy import sparse, io

out=Path(__file__).resolve().parents[1]/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42'
assert (out/'COMPUTATION_VERIFIED.json').exists()
started=time.time()
revision=out/'figure_revision_20260922'
revision.mkdir(exist_ok=True)
previous=revision/'previous'
if not previous.exists():
    previous.mkdir()
    for name in ['candidate_r1.1_embeddings.png','candidate_r1.1_markers.png','CANDIDATE_REVIEW_PLOTS.json','PUBLICATION_VERIFIED.json']:
        if (out/name).exists():shutil.copy2(out/name,previous/name)
    project=Path(__file__).resolve().parents[1]
    shutil.copy2(project/'Npp53Rb1_tracing_report_part1_preprocessing_tumor_states.html',previous/'report_before_figure_revision.html')
protected=['scvi_latent_PRIVATE.csv','umap_PRIVATE.csv','tsne_PRIVATE.csv','resolution_assignments_PRIVATE.csv','AMI_ranking.csv','top5_evidence.csv']
before={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected}
rank=pd.read_csv(out/'top5_evidence.csv')
assert rank.iloc[0].resolution==1.1 and rank.iloc[0].defended_clusters==19
labels=pd.read_csv(out/'resolution_assignments_PRIVATE.csv',index_col=0)['leiden_r1.1_seed42']
states=pd.read_csv(out/'cluster_state_inference/state_summary.csv')
state_by_cl={int(r.original_cluster[2:]):r.state for r in states.itertuples()}
rank_by_cl={int(r.original_cluster[2:]):r.rank for r in states.itertuples()}
qc=pd.read_csv(out/'expression_review/r1.1_QC.csv').set_index('cluster')
assert len(labels)==42806 and labels.nunique()==19
plt.rcParams.update({'font.size':14,'axes.titlesize':18})
palette=plt.get_cmap('tab20').colors
fig,axes=plt.subplots(1,2,figsize=(16,7),layout='constrained')
order=np.random.default_rng(42).permutation(len(labels))
meta=pd.read_csv(out/'metadata_PRIVATE.csv',index_col=0)
emb=pd.read_csv(out/'umap_PRIVATE.csv',index_col=0)
assert emb.index.equals(labels.index) and meta.index.equals(labels.index) and emb.shape==(42806,2)
xy=emb.to_numpy()
treatments={'Intact':'#278c74','Castration':'#e99b27','Castration+ENZA':'#9569b7'}
assert meta.treatment.value_counts().to_dict()=={'Intact':23840,'Castration+ENZA':13683,'Castration':5283}
for i,ax in enumerate(axes):
    colors=[palette[rank_by_cl[int(c)]-1] for c in labels.iloc[order]] if i==0 else meta.treatment.iloc[order].map(treatments).tolist()
    ax.scatter(xy[order,0],xy[order,1],c=colors,s=5,alpha=.85,linewidths=0,rasterized=True)
    if i==0:
        for cl in range(19):
            med=np.median(xy[labels.to_numpy()==cl],axis=0)
            ax.text(*med,f'c{rank_by_cl[cl]}',ha='center',va='center',fontsize=13,weight='bold',
                    path_effects=[pe.withStroke(linewidth=3,foreground='white')])
    else:
        handles=[Line2D([],[],marker='o',linestyle='',color=c,markersize=8,label=f'{t} (n={meta.treatment.eq(t).sum():,})') for t,c in treatments.items()]
        ax.legend(handles=handles,loc='upper right',frameon=True,fontsize=12)
    ax.set(title='UMAP · inferred state ID' if i==0 else 'UMAP · treatment',xlabel='UMAP_1',ylabel='UMAP_2',xticks=[],yticks=[])
    ax.spines[['top','right']].set_visible(False)
fig.suptitle('Gene-informed putative states • Leiden 1.1, seed 42 • 42,806 cells',fontsize=19)
fig.savefig(out/'candidate_r1.1_embeddings.png',dpi=180);plt.close(fig)

genes=['Krt8','Krt18','Epcam','Lcn2','Clu','Krt19','Nupr1','Krt5','Krt14','Col17a1',
       'Mecom','Ccnd2','Krt13','Krt10','Bnip3','Ndrg1','Mki67','Top2a','Txnip',
       'Lyz2','Tyrobp','C1qa','Igfbp2','Spp1','Krt7','Vim','Arhgap22','Frmd4a',
       'Crabp1','Id3','Cxcl12','Ccl2','Acta2','Tagln','Il33','Dcn','Inhba',
       'Col1a1','Col3a1','Col5a2','Ecrg4','GFP','mCherry']
# Actual Seurat DotPlot, not a recolouring of mean(log1p) summaries.
source=Path(__file__).resolve().parents[1]/'analysis/min400_20260921/tumor_umap_20260921/tumor_full_counts.h5ad'
a=ad.read_h5ad(source,backed='r')
positions=a.obs_names.get_indexer(labels.index)
assert (positions>=0).all() and set(genes)<=set(a.var_names)
x=a[:,genes].X.tocsr()[positions,:].astype(float)
a.file.close()
assert x.shape==(42806,len(genes)) and np.all(x.data>=0)
io.mmwrite(revision/'marker_counts.mtx',x.T)
pd.DataFrame({'gene':genes}).to_csv(revision/'marker_genes.csv',index=False)
pd.DataFrame({'cell':labels.index,'cluster':labels.to_numpy(),'total_counts':meta.total_counts.to_numpy()}).to_csv(revision/'marker_metadata_PRIVATE.csv',index=False)
rscript=Path(__file__).with_name('plot_filtered_dotplot_seurat.R')
with (revision/'seurat_plot.log').open('w') as handle:
    subprocess.run(['/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript',str(rscript),str(out)],stdout=handle,stderr=subprocess.STDOUT,check=True,timeout=600)
# Independent numerical reconciliation of actual Seurat plot data.
plotdata=pd.read_csv(revision/'seurat_dotplot_data.csv')
assert len(plotdata)==19*len(genes) and set(plotdata['features.plot'])==set(genes)
norm=sparse.diags(1e4/meta.total_counts.to_numpy()).dot(x).tocsr()
avg=np.vstack([np.asarray(norm[labels.to_numpy()==cl].mean(axis=0)).ravel() for cl in range(19)])
pct=np.vstack([np.asarray((x[labels.to_numpy()==cl]>0).mean(axis=0)).ravel()*100 for cl in range(19)])
logmean=np.log1p(avg)
scaled=np.clip((logmean-logmean.mean(axis=0))/logmean.std(axis=0,ddof=1),-2.5,2.5)
for row in plotdata.to_dict('records'):
    state=row['id'].split(' (n=')[0]
    cl=next(c for c,s in state_by_cl.items() if s==state);j=genes.index(row['features.plot'])
    assert np.isclose(row['avg.exp'],avg[cl,j],rtol=1e-10,atol=1e-10)
    assert np.isclose(row['pct.exp'],pct[cl,j],rtol=1e-10,atol=1e-10)
    assert np.isclose(row['avg.exp.scaled'],scaled[cl,j],rtol=1e-10,atol=1e-10)
assert before=={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected}
(out/'CANDIDATE_REVIEW_PLOTS.json').write_text(json.dumps(dict(cells=42806,resolution=1.1,seed=42,
    clusters=19,final_state_labels=False,inferred_state_labels=True,expression='Seurat DotPlot: log1p(mean normalized counts), gene-wise sample-SD scaling, clipped [-2.5,2.5]; original full-library denominators',
    panels=['UMAP cluster','UMAP treatment'],genes=genes,dotplot='Seurat DotPlot; lightgrey to #08306b; dot.scale7; scale.by radius',
    protected_input_sha256=before,elapsed_seconds=time.time()-started,
    plots=['candidate_r1.1_embeddings.png','candidate_r1.1_markers.png']),indent=2))
print('PASS: reordered state/treatment UMAP; actual Seurat marker dotplot including GFP/mCherry; independent mean/detection/scaling checks; scientific inputs unchanged.')

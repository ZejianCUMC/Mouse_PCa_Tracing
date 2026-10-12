"""Display completed descriptive island results without filtering cells."""
import argparse
from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
p=argparse.ArgumentParser();p.add_argument('output',type=Path);out=p.parse_args().output
dest=out/'expression_review'
meta=pd.read_csv(out/'metadata_PRIVATE.csv',index_col=0)
xy=pd.read_csv(out/'umap_unchanged_PRIVATE.csv',index_col=0).loc[meta.index]
marker=pd.read_csv(dest/'islands_canonical_markers.csv')
qc=pd.read_csv(dest/'island_matched_QC.csv')
regions=['Upper island','Far-right island'];colours=['#c43a9a','#12826c']
plt.rcParams.update({'font.size':13,'axes.titlesize':16})
fig=plt.figure(figsize=(17,10));gs=fig.add_gridspec(2,2,height_ratios=[1.25,1],width_ratios=[1,1.25],hspace=.42,wspace=.3)
ax=fig.add_subplot(gs[0,0]);ax.scatter(xy.UMAP_1,xy.UMAP_2,s=1,c='#dddddd',rasterized=True)
for region,colour in zip(regions,colours):
    pos=xy.loc[meta.review_region==region];ax.scatter(pos.UMAP_1,pos.UMAP_2,s=15,c=colour,label=f'{region} (n={len(pos)})')
ax.set(xlabel='UMAP_1',ylabel='UMAP_2',title='Review regions on unchanged scVI UMAP',aspect='equal')
ax.legend(loc='lower left',fontsize=11,markerscale=1.5)
ax=fig.add_subplot(gs[0,1])
metrics=['n_genes_by_counts','total_counts','pct_counts_mt','pct_counts_ribo']
for offset,region,colour in zip([-.18,.18],regions,colours):
    data=qc[qc.region==region].set_index('metric').loc[metrics]
    ratios=(data.island_median/data.matched_control_median).to_numpy()
    bars=ax.bar(np.arange(4)+offset,ratios,width=.33,color=colour,label=region)
    for bar,value in zip(bars,ratios):ax.text(bar.get_x()+bar.get_width()/2,value+.05,f'{value:.2f}×',ha='center',fontsize=12)
ax.axhline(1,c='black',ls='--',lw=1)
ax.set(xticks=range(4),xticklabels=['Detected\ngenes','UMI\ncounts','Mitochondrial\nRNA %','Ribosomal\nRNA %'],
       ylim=(0,2.65),ylabel='Island / matched-control median',title='QC against matched non-island cells')
ax.legend(fontsize=11,loc='upper left')
ax=fig.add_subplot(gs[1,:])
genes=['GFP','mCherry','Krt8','Krt18','Krt14','Vim','Acta2','Tagln','Ptprc','Pecam1','Foxi1','Atp6v1g3','Atp6v0d2','Aqp6','Clcnkb','Chga','Chgb','Syp','Insm1','Mki67']
for y,region in enumerate(regions):
    data=marker[marker.cluster==region].set_index('gene').loc[genes]
    points=ax.scatter(range(len(genes)),np.repeat(y,len(genes)),s=data.fraction_positive.clip(0,1)*300,
                      c=data.mean_log1p_10k,cmap='Blues',norm=Normalize(0,4),edgecolors='#888888',linewidths=.4)
ax.set(xticks=range(len(genes)),xticklabels=genes,yticks=[0,1],yticklabels=regions,ylim=(1.6,-.6),
       title='Marker evidence: dot area = cells expressing; colour = mean log-normalized expression')
ax.tick_params(axis='x',rotation=55)
for value in [.25,.5,1]:ax.scatter([],[],s=value*300,c='#777777',label=f'{value:.0%}')
ax.legend(title='Cells expressing',bbox_to_anchor=(1.005,1.1),loc='upper left',fontsize=11)
fig.colorbar(points,ax=ax,location='bottom',pad=.46,fraction=.09,shrink=.35,label='Mean log1p(counts per 10,000)')
fig.suptitle('Small-island review | all 258 cells retained; neither region independently proves malignancy',fontsize=18,y=.985)
fig.text(.5,.01,'Matched controls: same sample × parent Harmony cluster; excludes both islands. Descriptive comparison, not a doublet or malignancy test.',ha='center',fontsize=11)
fig.subplots_adjust(top=.9,bottom=.1,left=.11,right=.86)
fig.savefig(dest/'island_evidence.png',dpi=180)
print('Island evidence figure saved.')

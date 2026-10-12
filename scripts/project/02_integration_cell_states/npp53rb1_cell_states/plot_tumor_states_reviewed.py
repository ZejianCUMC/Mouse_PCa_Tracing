"""Plot completed state evidence using an explicitly reviewed cluster mapping.

This script neither chooses a resolution nor creates publication acceptance.
"""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

p=argparse.ArgumentParser();p.add_argument('output',type=Path)
out=p.parse_args().output
review=json.loads((out/'reviewed_annotation.json').read_text())
assert review['reviewed'] is True
assert (out/'scSHC_COMPLETE.json').exists()
r=float(review['resolution']);tag=f'r{r:.1f}'
mapping=review['clusters']
meta=pd.read_csv(out/'metadata_PRIVATE.csv',index_col=0)
xy=pd.read_csv(out/'umap_unchanged_PRIVATE.csv',index_col=0).loc[meta.index]
assign=pd.read_csv(out/'resolution_assignments_PRIVATE.csv',index_col=0,dtype=str).loc[meta.index]
cl=assign[f'leiden_{tag}_seed19']
assert set(cl)==set(mapping) and len(meta)==43057
labels=cl.map({k:v['label'] for k,v in mapping.items()})
order=list(dict.fromkeys(v['label'] for v in mapping.values()))
colours=dict(zip(order,plt.get_cmap('tab20').colors[:len(order)]))
assert len(order)<=20
plt.rcParams.update({'font.size':13,'axes.titlesize':16})

# Resolution evidence: show all tested resolutions, not only the winner.
rank=pd.read_csv(out/'AMI_ranking.csv').sort_values('resolution')
sweep=pd.read_csv(out/'resolution_sweep.csv')
tested=pd.read_csv(out/'scSHC_summary.csv').sort_values('resolution')
fig,axes=plt.subplots(1,3,figsize=(18,5.6),layout='constrained')
ax=axes[0]
for col,lab,alpha in [('AMI','Seed 19 (primary)',1),('AMI_seed29','Seed 29',.55),('AMI_seed39','Seed 39',.55)]:
    ax.plot(rank.resolution,rank[col],marker='o',label=lab,alpha=alpha)
ax.axvline(r,c='black',ls=':',label=f'Reviewed r={r:.1f}')
ax.set(xlabel='Leiden resolution',ylabel='Adjusted mutual information',title='Agreement with raw CHOIR')
ax.legend(fontsize=10)
ax=axes[1];xx=np.arange(len(tested))
ax.bar(xx-.18,tested.input_clusters,.36,label='Input Leiden clusters',color='#5b8db8')
ax.bar(xx+.18,tested.defended_clusters,.36,label='After sc-SHC',color='#dd9856')
for i,row in enumerate(tested.itertuples()):
    ax.text(i-.18,row.input_clusters+.2,str(row.input_clusters),ha='center')
    ax.text(i+.18,row.defended_clusters+.2,str(row.defended_clusters),ha='center')
ax.set(xticks=xx,xticklabels=[f'{x:.1f}' for x in tested.resolution],xlabel='Candidate resolution',ylabel='Number of clusters',title='Within-run sc-SHC checks')
ax.set_ylim(0,max(tested.input_clusters)*1.25);ax.legend(fontsize=10)
ax=axes[2]
ari_cols=[c for c in sweep if 'ARI' in c or 'ari' in c]
assert len(ari_cols)==3,ari_cols
for c in ari_cols:ax.plot(sweep.resolution,sweep[c],marker='o',label=c.replace('_',' '))
ax.axvline(r,c='black',ls=':');ax.set(xlabel='Leiden resolution',ylabel='Adjusted Rand index',ylim=(0,1.02),title='Seed sensitivity (same graph)')
ax.legend(fontsize=9)
fig.suptitle('Resolution review | statistical checks are not independent biological truth',fontsize=18)
fig.savefig(out/'state_resolution_review.png',dpi=180);plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(18,8))
for ax,values,title in [(axes[0],cl,f'Candidate clusters | Leiden r={r:.1f}'),(axes[1],labels,'Expression programs / QC flags')]:
    levels=sorted(set(cl),key=int) if ax is axes[0] else order
    palette={k:colours[mapping[k]['label']] for k in levels} if ax is axes[0] else colours
    for level in levels:
        pos=xy.loc[values==level]
        ax.scatter(pos.UMAP_1,pos.UMAP_2,s=4.5,c=[palette[level]],linewidths=0,rasterized=True,label=level)
        if ax is axes[0]:
            centre=pos.median();ax.text(centre.UMAP_1,centre.UMAP_2,level,fontsize=12,ha='center',va='center',bbox=dict(boxstyle='round,pad=.15',fc='white',ec='none',alpha=.82))
    ax.set(xlabel='UMAP_1',ylabel='UMAP_2',title=title,aspect='equal',xticks=[],yticks=[])
axes[1].legend(bbox_to_anchor=(1.01,1),loc='upper left',fontsize=10,markerscale=3,frameon=False)
fig.suptitle('Unchanged scVI NB run+sample UMAP | 43,057 cells retained',fontsize=19)
fig.subplots_adjust(left=.06,right=.77,wspace=.22,top=.89,bottom=.09)
fig.savefig(out/'state_umap.png',dpi=180,bbox_inches='tight');plt.close(fig)

# Combine descriptive marker aggregates with exact cluster-size weights.
markers=pd.read_csv(out/'candidate_markers'/f'{tag}_display_marker_summary.csv',dtype={'cluster':str})
markers=markers.drop_duplicates(['cluster','gene'])
sizes=cl.value_counts();markers['n']=markers.cluster.map(sizes)
markers['label']=markers.cluster.map({k:v['label'] for k,v in mapping.items()})
genes=review['dotplot_genes']
assert set(genes)<=set(markers.gene)
rows=[]
validated=pd.read_csv(out/f'scSHC_{tag}_labels_PRIVATE.csv',index_col=0,dtype=str).loc[meta.index]
for label in order:
    d=markers[markers.label==label].copy()
    for gene in genes:
        q=d[d.gene==gene]
        rows.append(dict(label=label,gene=gene,mean=np.average(q.mean_log1p_10k,weights=q.n),fraction=np.average(q.fraction_positive,weights=q.n)))
dots=pd.DataFrame(rows);dots.to_csv(out/'reviewed_state_marker_means.csv',index=False)
fig,ax=plt.subplots(figsize=(max(16,len(genes)*.36),max(6,len(order)*.48+3)))
norm=Normalize(0,max(1,float(dots['mean'].max())))
for y,label in enumerate(order):
    d=dots[dots.label==label].set_index('gene').loc[genes]
    points=ax.scatter(range(len(genes)),np.repeat(y,len(genes)),s=d.fraction.clip(0,1)*230,c=d['mean'],cmap='Blues',norm=norm,edgecolors='#888',linewidths=.3)
ax.set(xticks=range(len(genes)),xticklabels=genes,yticks=range(len(order)),yticklabels=order,ylim=(len(order)-.4,-.6),xlabel='Marker gene',ylabel='Putative state / review category',title='Marker evidence | dot area = cells expressing; colour = mean log-normalized expression')
ax.tick_params(axis='x',rotation=60,labelsize=11)
for f in [.25,.5,1]:ax.scatter([],[],s=f*230,c='#777',label=f'{f:.0%}')
ax.legend(title='Cells expressing',bbox_to_anchor=(1.01,1),loc='upper left')
fig.colorbar(points,ax=ax,location='bottom',pad=.34,fraction=.035,shrink=.4,label='Mean log1p(counts per 10,000)')
fig.subplots_adjust(left=.23,right=.87,top=.9,bottom=.14)
fig.savefig(out/'state_marker_dotplot.png',dpi=180,bbox_inches='tight');plt.close(fig)

rows=[]
for label in order:
    mask=labels==label;d=meta.loc[mask]
    clusters=sorted(set(cl[mask]),key=int)
    notes=list(dict.fromkeys(mapping[k]['evidence'] for k in clusters))
    sc_groups=', '.join(x.replace('new','G') for x in sorted(set(validated.loc[mask,'scSHC_cluster'])))
    rows.append(dict(putative_state_or_flag=label,cluster_IDs=', '.join(clusters),scSHC_group=sc_groups,cells=len(d),samples=d.sample_id.nunique(),runs=d.run_id.nunique(),reporter_positive=int(d.reporter_positive.sum()),largest_sample_percent=round(d.sample_id.value_counts().max()/len(d)*100,1),median_detected_genes=round(d.n_genes_by_counts.median(),1),evidence='; '.join(notes)))
summary=pd.DataFrame(rows);assert summary.cells.sum()==43057
assert summary.reporter_positive.sum()==39683
summary.to_csv(out/'reviewed_state_summary.csv',index=False)
pd.DataFrame({'cell':meta.index,'cluster':cl.to_numpy(),'reviewed_label':labels.to_numpy()}).to_csv(out/'reviewed_labels_PRIVATE.csv',index=False)
# Samples, not pooled treatment cells, are the units displayed here. This is
# descriptive composition; no treatment hypothesis test is claimed.
composition=pd.crosstab(meta.sample_id,labels).reindex(columns=order,fill_value=0)
sample_meta=meta.groupby('sample_id').agg(treatment=('treatment','first'),run_id=('run_id','first'))
assert (meta.groupby('sample_id').treatment.nunique()==1).all()
sample_meta['treatment_order']=sample_meta.treatment.map({'Intact':0,'Castration':1,'Castration+ENZA':2})
sample_meta['cells']=composition.sum(axis=1)
sample_meta=sample_meta.sort_values(['treatment_order','cells'],ascending=[True,False])
composition=composition.loc[sample_meta.index]
sample_meta.drop(columns='treatment_order').join(composition).to_csv(out/'reviewed_state_sample_counts.csv')
proportions=composition.div(composition.sum(axis=1),axis=0)
fig,ax=plt.subplots(figsize=(17,8));bottom=np.zeros(len(proportions));xx=np.arange(len(proportions))
for label in order:
    vals=100*proportions[label].to_numpy();ax.bar(xx,vals,bottom=bottom,color=colours[label],label=label,width=.85);bottom+=vals
for treatment in ['Intact','Castration','Castration+ENZA']:
    positions=np.flatnonzero(sample_meta.treatment.to_numpy()==treatment)
    ax.text(positions.mean(),109,treatment,ha='center',fontsize=14,fontweight='bold')
    if positions[0]>0:ax.axvline(positions[0]-.5,c='#777',lw=.7)
for i,n in enumerate(sample_meta.cells):ax.text(i,101,f'{n:,}',ha='center',fontsize=9,rotation=90,va='bottom')
ax.set(xticks=xx,xticklabels=proportions.index,ylim=(0,100),xlabel='Sample ID (grouped by treatment)',ylabel='Retained tumour-cohort cells (%)')
ax.tick_params(axis='x',rotation=60);ax.legend(bbox_to_anchor=(1.01,1),loc='upper left',fontsize=10,frameon=False)
fig.suptitle('Per-sample composition | descriptive; no treatment-effect test',fontsize=19,y=.99)
fig.subplots_adjust(left=.07,right=.77,bottom=.18,top=.8)
fig.savefig(out/'state_sample_composition.png',dpi=180,bbox_inches='tight');plt.close(fig)
print('Reviewed evidence plotted; publication acceptance remains a separate manual gate.')

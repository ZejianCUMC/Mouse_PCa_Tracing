"""Compact joint panels; fixed coordinates, shared state palette and limits."""
import os,tempfile
os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='rb1-split-mpl-'))
from pathlib import Path
import json,hashlib,shutil
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
p=Path(__file__).resolve().parents[1]
out=p/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/exclude_current_c17_c19_seed42'
meta=pd.read_csv(out/'metadata_PRIVATE.csv',index_col=0)
lab=pd.read_csv(out/'state_labels_PRIVATE.csv',index_col=0).inferred_state
states=sorted(lab.unique(),key=lambda x:int(x.split('_')[0][1:]))
assert len(lab)==40141 and len(states)==16 and lab.index.equals(meta.index)
colours={s:plt.get_cmap('tab20').colors[int(s.split('_')[0][1:])-1] for s in states}
plt.rcParams.update({'font.size':12,'axes.titlesize':15,'axes.labelsize':12})
treat=['Intact','Castration','Castration+ENZA'];runs=sorted(meta.run_id.unique())
spec=[('diffmap_PRIVATE.csv','final_diffmap_treatment.png',[('All cells',np.ones(len(meta),dtype=bool))]+[(t,meta.treatment.eq(t).to_numpy()) for t in treat],1,4,4.2,'DC'),
      ('umap_PRIVATE.csv','final_umap_treatment.png',[(t,meta.treatment.eq(t).to_numpy()) for t in treat],1,3,6.5,'UMAP'),
      ('umap_PRIVATE.csv','final_umap_run.png',[(r,meta.run_id.eq(r).to_numpy()) for r in runs],2,3,9.0,'UMAP')]
manifest=[]
for src,name,groups,nr,nc,height,prefix in spec:
 coords=pd.read_csv(out/src,index_col=0);assert coords.index.equals(meta.index)
 xy=coords.to_numpy();assert np.isfinite(xy).all()
 fig,axes=plt.subplots(nr,nc,figsize=(18,height),squeeze=False)
 full_names=prefix=='UMAP'
 legend_labels=[s.replace('_',' · ',1) if full_names else s.split('_')[0] for s in states]
 handles=[Line2D([],[],marker='o',linestyle='',color=colours[s],markersize=6,label=label) for s,label in zip(states,legend_labels)]
 # A dedicated bottom band keeps all four rows of full names outside axes.
 fig.subplots_adjust(left=.045,right=.99,top=.91 if full_names else .88,bottom=1.85/height if full_names else .23,wspace=.24,hspace=.55)
 rng=np.random.default_rng(42)
 lo=xy.min(axis=0);hi=xy.max(axis=0);pad=(hi-lo)*.04
 counts=[]
 for ax,(group,mask) in zip(axes.flat,groups):
  idx=np.flatnonzero(mask);idx=rng.permutation(idx);counts.append(int(len(idx)))
  ax.scatter(xy[idx,0],xy[idx,1],c=lab.iloc[idx].map(colours).tolist(),s=3,alpha=.85,linewidths=0,rasterized=True)
  ax.set(title=f'{group} · n={len(idx):,}',xlabel=prefix+'_1',ylabel=prefix+'_2',xticks=[],yticks=[],
         xlim=(lo[0]-pad[0],hi[0]+pad[0]),ylim=(lo[1]-pad[1],hi[1]+pad[1]))
  ax.spines[['top','right']].set_visible(False)
 legend=fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.02 if full_names else .005),ncol=4 if full_names else 8,frameon=False,fontsize=11,columnspacing=1.6,labelspacing=.65 if full_names else .5)
 fig.canvas.draw();renderer=fig.canvas.get_renderer();b=legend.get_window_extent(renderer)
 assert all(not b.overlaps(ax.get_window_extent(renderer)) for ax in axes.flat)
 assert all(not b.overlaps(ax.get_tightbbox(renderer)) for ax in axes.flat), 'Legend overlaps axis text'
 assert fig.bbox.contains(b.x0,b.y0) and fig.bbox.contains(b.x1,b.y1)
 boxes=[t.get_window_extent(renderer) for t in legend.get_texts()]
 assert all(not left.overlaps(right) for i,left in enumerate(boxes) for right in boxes[i+1:]), 'Legend labels overlap'
 fig.savefig(out/name,dpi=180);plt.close(fig)
 manifest.append(dict(file=name,width_inches=18,height_inches=height,panels=len(groups),groups=[g for g,_ in groups],cell_counts=counts,coordinate_sha256=hashlib.sha256((out/src).read_bytes()).hexdigest(),colour='putative_state',same_limits=True,full_state_names=full_names,legend_labels=legend_labels,legend_overlap_and_clipping_checks=True))
# The newer full-name legend request supersedes the former total-height target.
ratio=sum(x['height_inches']/x['width_inches'] for x in manifest)/(12/23)
(out/'SPLIT_PLOTS_COMPLETE.json').write_text(json.dumps(dict(cells=40141,states=16,plots=manifest,combined_height_relative_to_original_fig3_8=ratio,colour_mapping={s:list(c) for s,c in colours.items()}),indent=2))
print('PASS: compact Diffmap merged/treatment, UMAP treatment and run panels; common state colours; combined height ratio',ratio)

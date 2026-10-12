"""Render current reviewed UMAP without modifying cells, names or coordinates."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def render(a, destination):
    states=sorted(a.obs.putative_state.unique(),key=lambda s:int(s.split('_')[0][1:]))
    palette={s:matplotlib.colors.to_hex(plt.get_cmap('tab20')(i)) for i,s in enumerate(states)}
    treatment={'Intact':'#4c78a8','Castration':'#f58518','Castration+ENZA':'#54a24b'}
    runs={s:matplotlib.colors.to_hex(plt.get_cmap('tab20')(i)) for i,s in enumerate(sorted(a.obs.run_id.unique()))}
    xy=a.obsm['X_umap'];order=np.random.default_rng(42).permutation(len(a))
    fig,axes=plt.subplots(1,3,figsize=(23,11))
    for ax,col,colors in zip(axes,['putative_state','treatment','run_id'],[palette,treatment,runs]):
        ax.scatter(xy[order,0],xy[order,1],c=[colors[v] for v in a.obs[col].iloc[order]],s=3,alpha=.85,linewidths=0,rasterized=True)
        ax.set_xlabel('UMAP_1',fontsize=13);ax.set_ylabel('UMAP_2',fontsize=13)
        ax.set_xticks([]);ax.set_yticks([])
        if col=='putative_state':
            for state in colors:
                ax.text(*np.median(xy[a.obs[col].eq(state)],axis=0),state.split('_')[0],fontsize=10,ha='center',weight='bold',bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=1))
    for ax,colors,title,ncols in zip(axes,[palette,treatment,runs],['Putative cell state','Treatment','Run'],[1,1,2]):
        handles=[Line2D([],[],marker='o',linestyle='',color=v,markersize=7,label=k) for k,v in colors.items()]
        ax.legend(handles=handles,title=title,title_fontsize=14,fontsize=14,ncol=ncols,loc='upper center',bbox_to_anchor=(.5,-.12),frameon=False,columnspacing=1.8,handletextpad=.5,labelspacing=.65)
    fig.subplots_adjust(left=.04,right=.99,top=.98,bottom=.47,wspace=.22)
    fig.savefig(destination,dpi=170,bbox_inches='tight');plt.close(fig)

if __name__=='__main__':
    import anndata as ad
    p=Path(__file__).resolve().parents[1]
    d=p/'analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42/approved_merges'
    render(ad.read_h5ad(d/'named_graph_PRIVATE.h5ad'),d/'reclustered_joint_umap.png')

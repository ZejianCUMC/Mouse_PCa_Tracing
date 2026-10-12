#!/usr/bin/env python
"""Genes-only trajectory heatmap: drop ribosomal + non-protein-coding/predicted genes,
show top 25 UP + 25 DOWN by Spearman rho vs corrected pseudotime."""
import warnings; warnings.filterwarnings('ignore')
import os, re, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

ROOT='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
OUT=os.environ['OUTDIR']; PT=os.environ['PTCSV']; NB=20
STATE_ORDER=['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR','Luminal-Persistent','Luminal-Squamous',
             'Basal-like','Partial-EMT','Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
SCOL={'Luminal-Distal':'#1f5fa8','Luminal-Proximal':'#9467bd','Luminal-ARdep':'#d23030','Luminal-CR':'#f1903d',
 'Luminal-Persistent':'#117c70','Luminal-Squamous':'#f78ec7','Basal-like':'#8c564b','Partial-EMT':'#bcbd22',
 'Mes-Quiescent':'#c49c94','Mes-Collagen':'#5b5b5b','Mes-Stemness':'#8c1a1a','Mes-myCAF':'#4daf4a'}
RIBO=re.compile(r'^(?:Rp[sl]\d|Rplp\d|Mrp[sl]\d|Rpsa$|Fau$|Uba52$|Rack1$)')
NONC=re.compile(r'(?:^Gm\d|Rik$|^mt-|^Mir\d|^Sno|^Snhg|^Malat1$|^Neat1$|^Xist$|^Rn\d|^A[CFL]\d{5,}|^Rpph|^Vaultrc|-ps\d?$)')

pt=pd.read_csv(PT); pt['k']=pt['sample']+'|'+pt['cellBC']
k2pt=dict(zip(pt['k'],pt['pt_rank'])); k2st=dict(zip(pt['k'],pt['state']))
a=sc.read_h5ad(f'{ROOT}/data/npp53_merge/npp53_edited.h5ad')
a.obs['bc']=[re.search(r'[ACGTN]{14,18}',x).group(0) for x in a.obs_names]
a.obs['k']=a.obs['sample_id'].astype(str)+'|'+a.obs['bc']
a.obs['pt']=a.obs['k'].map(k2pt); a.obs['st']=a.obs['k'].map(k2st)
a=a[a.obs['pt'].notna()].copy()
ptv=a.obs['pt'].values.astype(float); bins=np.clip((ptv*NB).astype(int),0,NB-1); stv=a.obs['st'].values

g=pd.read_csv(f'{OUT}/evopath_genes_ranked.csv')
mask=~(g.gene.str.contains(RIBO)|g.gene.str.contains(NONC))
gf=g[mask].copy()
print(f'genes: {len(g)} -> {len(gf)} after dropping ribosomal + non-coding',flush=True)
up=list(gf.sort_values('rho_cohort',ascending=False).head(25)['gene'])
dn=list(gf.sort_values('rho_cohort').head(25)['gene'])
sel=up+dn
gsel=[x for x in sel if x in set(a.var_names)]

# per-cell expression for selected genes
X=a[:,gsel].X; X=X.toarray() if hasattr(X,'toarray') else np.asarray(X)
M=np.full((len(gsel),NB),np.nan)
for r in range(len(gsel)):
    v=X[:,r]
    for b in range(NB):
        mm=bins==b
        if mm.sum()>=10: M[r,b]=np.nanmean(v[mm])
Mz=(M-np.nanmean(M,1,keepdims=True))/(np.nanstd(M,1,keepdims=True)+1e-9)
peak=np.nanargmax(Mz,1); order=np.argsort(peak); Mz=Mz[order]; labs=[gsel[i] for i in order]

comp=np.zeros((len(STATE_ORDER),NB))
for b in range(NB):
    mm=bins==b
    if mm.sum():
        vc=pd.Series(stv[mm]).value_counts(normalize=True)
        for i,s in enumerate(STATE_ORDER): comp[i,b]=vc.get(s,0)

cellh=0.26; feats=labs
CID={'Luminal-Distal':'c1','Luminal-Proximal':'c2','Luminal-ARdep':'c3','Luminal-CR':'c4','Luminal-Persistent':'c5','Luminal-Squamous':'c6','Basal-like':'c7','Partial-EMT':'c8','Mes-Quiescent':'c9','Mes-Collagen':'c10','Mes-Stemness':'c11','Mes-myCAF':'c12'}
TIN,GAP,SIN,BOT=0.55,0.12,0.85,1.55
heat=len(feats)*cellh; H=TIN+heat+GAP+SIN+BOT
L,R=0.06,0.87; W=R-L
fig=plt.figure(figsize=(12,H))
axH=fig.add_axes([L,1-(TIN+heat)/H,W,heat/H])
axC=fig.add_axes([L,BOT/H,W,SIN/H])
im=axH.imshow(Mz,aspect='auto',cmap='RdBu_r',vmin=-1.8,vmax=1.8)
axH.set_xticks([]); axH.set_yticks(range(len(feats))); axH.set_yticklabels(feats,fontsize=8); axH.yaxis.tick_right()
fig.text(L,1-(TIN*0.42)/H,'Genes following the tumor-evolution trajectory (top 25 up / 25 down; cohort)',fontsize=13,fontweight='bold',color='#111827',ha='left',va='center')
bottom=np.zeros(NB)
for i,st in enumerate(STATE_ORDER):
    axC.bar(range(NB),comp[i],bottom=bottom,width=1.0,color=SCOL[st],edgecolor='none'); bottom+=comp[i]
axC.set_xlim(-0.5,NB-0.5); axC.set_ylim(0,1); axC.set_yticks([]); axC.set_xticks([0,NB-1]); axC.set_xticklabels(['early (initiator)','late (terminal)'],fontsize=9)
axC.set_xlabel('Corrected PhyloVelo pseudotime  \u2192',fontsize=10)
fig.legend(handles=[Patch(color=SCOL[st],label=f'{CID[st]}_{st}') for st in STATE_ORDER],loc='lower center',ncol=6,fontsize=8,frameon=False,bbox_to_anchor=(0.5,0.60/H))
cax=fig.add_axes([0.42,0.25/H,0.16,0.10/H]); plt.colorbar(im,cax=cax,orientation='horizontal').set_label('row z-score',fontsize=8)
for fmt in ['png','pdf']: fig.savefig(f'{OUT}/evopath_GENES.{fmt}',dpi=170,bbox_inches='tight',facecolor='white')
print('saved evopath_GENES.png rows=',len(gsel));
print('UP:',up); print('DOWN:',dn); print('DONE',flush=True)

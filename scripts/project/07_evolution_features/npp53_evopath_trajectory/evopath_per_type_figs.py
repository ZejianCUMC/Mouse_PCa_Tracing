#!/usr/bin/env python
"""Four separate trajectory heatmaps (genes / pathways / regulons / FSU modules), each
= features x pseudotime-bin z-score (cohort), rows ordered by peak, with state-composition strip."""
import warnings; warnings.filterwarnings('ignore')
import os, re, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

ROOT='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
AN=f'{ROOT}/analysis'; OUT=os.environ['OUTDIR']; PT=os.environ['PTCSV']; NB=20
STATE_ORDER=['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR','Luminal-Persistent','Luminal-Squamous',
             'Basal-like','Partial-EMT','Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
SCOL={'Luminal-Distal':'#1f5fa8','Luminal-Proximal':'#9467bd','Luminal-ARdep':'#d23030','Luminal-CR':'#f1903d',
 'Luminal-Persistent':'#117c70','Luminal-Squamous':'#f78ec7','Basal-like':'#8c564b','Partial-EMT':'#bcbd22',
 'Mes-Quiescent':'#c49c94','Mes-Collagen':'#5b5b5b','Mes-Stemness':'#8c1a1a','Mes-myCAF':'#4daf4a'}
TITLE={'GENE':'Genes','PATH':'Hallmark / signature pathways (AUCell)','REG':'Regulons (Npp53 pyscenic AUCell)','FSU':'FSU Hotspot modules'}
ACCENT={'GENE':'#111827','PATH':'#2563eb','REG':'#b45309','FSU':'#059669'}

pt=pd.read_csv(PT); pt['k']=pt['sample']+'|'+pt['cellBC']
k2pt=dict(zip(pt['k'],pt['pt_rank'])); k2st=dict(zip(pt['k'],pt['state']))
a=sc.read_h5ad(f'{ROOT}/data/npp53_merge/npp53_edited.h5ad')
a.obs['bc']=[re.search(r'[ACGTN]{14,18}',x).group(0) for x in a.obs_names]
a.obs['k']=a.obs['sample_id'].astype(str)+'|'+a.obs['bc']
a.obs['pt']=a.obs['k'].map(k2pt); a.obs['st']=a.obs['k'].map(k2st)
a=a[a.obs['pt'].notna()].copy()
ptv=a.obs['pt'].values.astype(float); bins=np.clip((ptv*NB).astype(int),0,NB-1); stv=a.obs['st'].values

rs=pd.read_csv(f'{OUT}/evopath_sets_ranked.csv')
NE_EXCLUDE={'PATH:Aggarwal_NEPC_UP','PATH:Aggarwal_NEPC_DN','PATH:Dehm_NE','PATH:Tang_CRPC_NE','PATH:Tang_CRPC_SCL','PATH:NE_FHCRC','PATH:Lundberg_NE','PATH:Beltran_NEPC_UP'}
rs=rs[~rs.feature.isin(NE_EXCLUDE)].copy()  # remove neuroendocrine sigs (no NE in system)
reg=pd.read_csv(f'{AN}/pyscenic_out/auc_mtx.csv',index_col=0)
pth=pd.read_csv(f'{AN}/aucell_arhi/aucell_scores.csv',index_col=0)
fsu_df=pd.read_csv(f'{AN}/hotspot_canonical/cohort_genes_by_module_Intact.csv'); present=set(a.var_names)
for m in sorted(fsu_df.Module.unique()):
    gg=[x for x in fsu_df[fsu_df.Module==m]['Gene'] if x in present]; sc.tl.score_genes(a,gg,score_name=f'FSU:M{m}',use_raw=False)

# state composition per bin (shared)
comp=np.zeros((len(STATE_ORDER),NB))
for b in range(NB):
    mm=bins==b
    if mm.sum():
        vc=pd.Series(stv[mm]).value_counts(normalize=True)
        for i,s in enumerate(STATE_ORDER): comp[i,b]=vc.get(s,0)

def per_cell(feat):
    typ,nm=feat.split(':',1)
    if typ=='REG': return reg[nm].reindex(a.obs_names).values.astype(float)
    if typ=='PATH': return pth[nm].reindex(a.obs_names).values.astype(float)
    if typ=='FSU': return a.obs[feat].values.astype(float)
    ix=list(a.var_names).index(nm); col=a[:,nm].X; col=col.toarray().ravel() if hasattr(col,'toarray') else np.asarray(col).ravel()
    return col

def draw(typ, feats, labels, fname, cellh=0.30, sort_by_peak=True):
    M=np.full((len(feats),NB),np.nan)
    for r,f in enumerate(feats):
        v=per_cell(f)
        for b in range(NB):
            mm=bins==b
            if mm.sum()>=10: M[r,b]=np.nanmean(v[mm])
    Mz=(M-np.nanmean(M,1,keepdims=True))/(np.nanstd(M,1,keepdims=True)+1e-9)
    if sort_by_peak:
        peak=np.nanargmax(Mz,1); order=np.argsort(peak); Mz=Mz[order]; labs=[labels[i] for i in order]
    else:
        labs=list(labels)
    CID={'Luminal-Distal':'c1','Luminal-Proximal':'c2','Luminal-ARdep':'c3','Luminal-CR':'c4','Luminal-Persistent':'c5','Luminal-Squamous':'c6','Basal-like':'c7','Partial-EMT':'c8','Mes-Quiescent':'c9','Mes-Collagen':'c10','Mes-Stemness':'c11','Mes-myCAF':'c12'}
    TIN,GAP,SIN,BOT=0.55,0.12,0.85,1.55
    heat=len(labs)*cellh; H=TIN+heat+GAP+SIN+BOT
    L,R=0.06,0.87; W=R-L
    fig=plt.figure(figsize=(12,H))
    axH=fig.add_axes([L,1-(TIN+heat)/H,W,heat/H])
    axC=fig.add_axes([L,BOT/H,W,SIN/H])
    im=axH.imshow(Mz,aspect='auto',cmap='RdBu_r',vmin=-1.8,vmax=1.8)
    axH.set_xticks([]); axH.set_yticks(range(len(labs))); axH.set_yticklabels(labs,fontsize=8); axH.yaxis.tick_right()
    fig.text(L,1-(TIN*0.42)/H,f'{TITLE[typ]} following the tumor-evolution trajectory (corrected pseudotime, cohort)',fontsize=11.5,fontweight='bold',color=ACCENT[typ],ha='left',va='center')
    bottom=np.zeros(NB)
    for i,s in enumerate(STATE_ORDER):
        axC.bar(range(NB),comp[i],bottom=bottom,width=1.0,color=SCOL[s],edgecolor='none'); bottom+=comp[i]
    axC.set_xlim(-0.5,NB-0.5); axC.set_ylim(0,1); axC.set_yticks([]); axC.set_xticks([0,NB-1]); axC.set_xticklabels(['early (initiator)','late (terminal)'],fontsize=9)
    axC.set_xlabel('Corrected PhyloVelo pseudotime  \u2192',fontsize=10)
    fig.legend(handles=[Patch(color=SCOL[s],label=f'{CID[s]}_{s}') for s in STATE_ORDER],loc='lower center',ncol=6,fontsize=8,frameon=False,bbox_to_anchor=(0.5,0.60/H))
    cax=fig.add_axes([0.42,0.25/H,0.16,0.10/H]); plt.colorbar(im,cax=cax,orientation='horizontal').set_label('row z-score',fontsize=8)
    for fmt in ['png','pdf']: fig.savefig(f'{OUT}/{fname}.{fmt}',dpi=170,bbox_inches='tight',facecolor='white')
    plt.close(fig); print(f'saved {fname}.png ({len(feats)} rows)',flush=True)

# pathways: top 15 up + 15 down
pu=list(rs[rs.type=='PATH'].sort_values('rho_cohort',ascending=False).head(15)['feature']); pd_=list(rs[rs.type=='PATH'].sort_values('rho_cohort').head(15)['feature'])
draw('PATH',pu+pd_,[x.split(':',1)[1].replace('HALLMARK_','') for x in pu+pd_],'evopath_PATHWAYS')
# regulons: top 15 up + 15 down
ru=list(rs[rs.type=='REG'].sort_values('rho_cohort',ascending=False).head(15)['feature']); rd=list(rs[rs.type=='REG'].sort_values('rho_cohort').head(15)['feature'])
draw('REG',ru+rd,[x.split(':',1)[1] for x in ru+rd],'evopath_REGULONS')
# FSU: pass modules in canonical M1..M5 order (skip peak-based reordering) so the
# row order matches the per-arm Fig 10.3.14 and is easy to read across figures.
fu_ranked = list(rs[rs.type=='FSU']['feature'])
fu = sorted(fu_ranked, key=lambda x: int(x.split(':M',1)[1]))
draw('FSU',fu,[x.split(':',1)[1] for x in fu],'evopath_MODULES',cellh=0.55,sort_by_peak=False)
print('DONE',flush=True)

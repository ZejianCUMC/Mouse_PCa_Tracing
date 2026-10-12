#!/usr/bin/env python
"""Npp53 data-driven COHORT-WIDE consensus fate map (single panel, all 25 samples pooled).
Replacement candidate for report Fig 10.3.8 (Step 6 integrated view). Same method as Fig 10.3.5:
direction = pooled PhyloVelo pseudotime (Mann-Whitney) requiring cross-tumor agreement.
Edge types: solid single arrow=confident; solid double (<->)=strong link/unresolved; dashed arrow=weak lean/bridge;
headless dashed(- - -)=below-threshold/connectivity."""
import os, json
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerPatch
from scipy.stats import mannwhitneyu

BASE='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
IN_JSON='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12/biology_review/datadriven_fatemap/consensus_fate_per_arm.json'   # step0 copy: read the canonical JSON in place
H5AD=f'{BASE}/data/npp53_merge/npp53_edited.h5ad'
PT=os.environ['PTCSV']; OUTDIR=os.environ.get('OUTDIR','.')
MINC=20
STATE_COLORS={'Luminal-Distal':'#1f5fa8','Luminal-Proximal':'#9467bd','Luminal-ARdep':'#d23030','Luminal-CR':'#f1903d',
    'Luminal-Persistent':'#117c70','Luminal-Squamous':'#f78ec7','Basal-like':'#8c564b','Partial-EMT':'#bcbd22',
    'Mes-Quiescent':'#c49c94','Mes-Collagen':'#5b5b5b','Mes-Stemness':'#8c1a1a','Mes-myCAF':'#4daf4a'}
LB={'Luminal-Distal':'Lum-\nDistal','Luminal-Proximal':'Lum-\nProximal','Luminal-ARdep':'Lum-\nARdep','Luminal-CR':'Lum-\nCR',
    'Luminal-Persistent':'Lum-\nPersistent','Luminal-Squamous':'Lum-\nSquamous','Basal-like':'Basal-\nlike','Partial-EMT':'Partial-\nEMT',
    'Mes-Quiescent':'Mes-\nQuiescent','Mes-Collagen':'Mes-\nCollagen','Mes-Stemness':'Mes-\nStemness','Mes-myCAF':'Mes-\nmyCAF'}
SI={'Luminal-Distal':1,'Luminal-Proximal':2,'Luminal-ARdep':3,'Luminal-CR':4,'Luminal-Persistent':5,'Luminal-Squamous':6,
    'Basal-like':7,'Partial-EMT':8,'Mes-Quiescent':9,'Mes-Collagen':10,'Mes-Stemness':11,'Mes-myCAF':12}
SH={s:f'c{SI[s]}' for s in SI}
POS={'Luminal-Distal':(0.50,0.84),'Luminal-Proximal':(0.24,0.72),'Luminal-Squamous':(0.12,0.54),'Basal-like':(0.50,0.55),
    'Luminal-ARdep':(0.76,0.72),'Luminal-CR':(0.88,0.54),'Luminal-Persistent':(0.76,0.48),'Partial-EMT':(0.32,0.32),
    'Mes-Quiescent':(0.68,0.32),'Mes-Collagen':(0.24,0.14),'Mes-Stemness':(0.50,0.10),'Mes-myCAF':(0.76,0.14)}
MES={'Partial-EMT','Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF'}
FALLBACK={'Luminal-Squamous':['Luminal-Proximal','Basal-like'],'Luminal-Persistent':['Basal-like'],
    'Partial-EMT':['Basal-like','Luminal-Persistent'],'Mes-Quiescent':['Partial-EMT','Basal-like'],
    'Mes-Collagen':['Partial-EMT','Mes-Quiescent','Basal-like'],'Mes-myCAF':['Mes-Quiescent','Partial-EMT','Basal-like'],
    'Luminal-CR':['Luminal-Distal'],'Luminal-ARdep':['Luminal-Distal'],
    'Mes-Stemness':['Mes-Quiescent','Luminal-CR','Partial-EMT','Basal-like','Luminal-Distal']}
CANON=[('Luminal-Distal','Basal-like'),('Luminal-Distal','Luminal-Proximal'),('Luminal-Proximal','Luminal-Squamous'),
    ('Luminal-Proximal','Basal-like'),('Luminal-Squamous','Basal-like'),('Basal-like','Luminal-Persistent'),
    ('Luminal-Distal','Luminal-Persistent'),('Basal-like','Partial-EMT'),('Luminal-Persistent','Partial-EMT'),
    ('Partial-EMT','Mes-Quiescent'),('Partial-EMT','Mes-Collagen'),('Mes-Quiescent','Mes-myCAF'),
    ('Luminal-Distal','Luminal-ARdep'),('Luminal-Distal','Luminal-CR'),('Mes-Quiescent','Mes-Stemness')]

data=json.loads(open(IN_JSON).read())
# cohort terminal frequency = sum over arms / total samples
N_TOT=sum(data['n_by_arm'].values())
in_term={s:sum(data['in_terminal_by_arm'][a].get(s,0) for a in data['n_by_arm']) for s in STATE_COLORS}
cooccur_all={}
for a in data['n_by_arm']:
    for k,v in data['cooccur_by_arm'][a].items():
        cooccur_all[k]=cooccur_all.get(k,0)+v
def cf(a,b):
    key='|'.join(sorted([a,b],key=list(STATE_COLORS.keys()).index))
    return cooccur_all.get(key,0)/max(1,N_TOT)
adobj=ad.read_h5ad(H5AD,backed='r'); ob=adobj.obs[['sample_id','putative_state']].copy()
ct=ob.groupby(['sample_id','putative_state']).size().reset_index(name='n'); mx=ct.groupby('putative_state')['n'].max().to_dict()
faded={s for s in STATE_COLORS if mx.get(s,0)<MINC}
pdf=pd.read_csv(PT)   # all arms
g=pdf.groupby(['sample','state']).agg(n=('pt_rank','size'),m=('pt_rank','mean')).reset_index()
mpiv=g.pivot(index='sample',columns='state',values='m'); npiv=g.pivot(index='sample',columns='state',values='n')
cells={s:pdf[pdf.state==s]['pt_rank'].values for s in STATE_COLORS}
nsamp={s:(pdf[pdf.state==s].groupby('sample').size()>20).sum() for s in STATE_COLORS}
def both(a,b):
    if a not in npiv or b not in npiv: return []
    return [s for s in npiv.index if npiv.loc[s,a]>20 and npiv.loc[s,b]>20]
def pooled(a,b):
    da,db=cells[a],cells[b]
    if len(da)<5 or len(db)<5 or nsamp[a]<2 or nsamp[b]<2: return None
    U,pv=mannwhitneyu(da,db,alternative='two-sided');auc=U/(len(da)*len(db));eff=abs(auc-0.5)
    direc=(a,b) if np.median(da)<=np.median(db) else (b,a)
    bo=both(a,b); npair=len(bo)
    agree=(sum(1 for s in bo if ((mpiv.loc[s,a]<mpiv.loc[s,b])==(direc==(a,b))))/npair) if npair else None
    strong=pv<0.01 and eff>=0.10
    if strong and npair>=2 and agree<=0.5: return (direc,'bidir')
    if strong and (npair<2 or agree>0.5): return (direc,'solid')
    if pv<0.05 and eff>=0.05 and (npair<2 or agree>=0.5): return (direc,'lean')
    return (direc,'headless')
def nrad(s): return 0.032+0.052*(in_term.get(s,0)/max(1,N_TOT))
def circle(ax,s,r,fd):
    x,y=POS[s]; c=STATE_COLORS[s]; ec=c if not fd else '#bbb'; lw=2.8 if not fd else 1.4; al=1 if not fd else 0.35
    ax.add_patch(Circle((x,y),r,facecolor='white',edgecolor=ec,linewidth=lw,zorder=5,alpha=al))
    ax.text(x,y+r*0.14,LB[s],ha='center',va='center',fontsize=8.4,color='#111' if not fd else '#aaa',fontweight='bold',zorder=6,linespacing=0.9,alpha=1 if not fd else 0.55)
    by=y-r-0.016
    ax.add_patch(FancyBboxPatch((x-0.019,by-0.015),0.038,0.028,boxstyle='round,pad=0.001,rounding_size=0.004',facecolor='white',edgecolor='#666' if not fd else '#ccc',linewidth=0.9,zorder=6,alpha=1 if not fd else 0.5))
    ax.text(x,by-0.001,f'{SI[s]}',ha='center',va='center',fontsize=10,fontweight='bold',color='#111' if not fd else '#bbb',zorder=7)
EDGES=[]
def edge(ax,src,dst,w,color,ls,rad,rs,rd,double=False,headless=False):
    EDGES.append(dict(src=src,dst=dst,w=round(w,3),color=color,double=double,headless=headless,dashed=ls!='-'))
    x0,y0=POS[src];x1,y1=POS[dst];rs+=0.004;rd+=0.009;dx,dy=x1-x0,y1-y0;L=(dx**2+dy**2)**.5
    if L<1e-6:return
    style='<|-|>' if double else ('-' if headless else '-|>')
    ax.add_patch(FancyArrowPatch((x0+dx/L*rs,y0+dy/L*rs),(x1-dx/L*rd,y1-dy/L*rd),arrowstyle=style,color=color,linewidth=1.0+3.2*w,linestyle=ls,mutation_scale=18,connectionstyle=f'arc3,rad={rad}',zorder=3))
node_r={s:nrad(s) for s in STATE_COLORS}; avoid=list(STATE_COLORS)
def cfree(src,dst):
    def bez(rad,n=15):
        x0,y0=POS[src];x1,y1=POS[dst];dx,dy=x1-x0,y1-y0;mx_,my_=(x0+x1)/2,(y0+y1)/2;cx,cy=mx_+rad*dy,my_-rad*dx
        return [((1-t/n)**2*x0+2*(t/n)*(1-t/n)*cx+(t/n)**2*x1,(1-t/n)**2*y0+2*(t/n)*(1-t/n)*cy+(t/n)**2*y1) for t in range(1,n)]
    for rad in [0.0,0.22,-0.22,0.40,-0.40,0.55,-0.55]:
        ok=True
        for nm in avoid:
            if nm in (src,dst):continue
            nx,ny=POS[nm]
            if any(((px-nx)**2+(py-ny)**2)**.5<0.08 for px,py in bez(rad)):ok=False;break
        if ok:return rad
    return 0.55

fig,ax=plt.subplots(figsize=(11,11)); ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_aspect('equal');ax.axis('off')
ax.text(0.50,0.955,'Tumor initiation',ha='center',va='center',fontsize=11,fontweight='bold',color='#111',zorder=6)
r1=nrad('Luminal-Distal'); ax.add_patch(FancyArrowPatch((0.50,0.84+r1+0.028),(0.50,0.84+r1+0.006),arrowstyle='-|>',color='#111',lw=2.0,mutation_scale=14,zorder=3))
ax.add_patch(FancyBboxPatch((0.02,0.965),0.96,0.045,boxstyle='round,pad=0.002,rounding_size=0.01',facecolor='#0f2747',edgecolor='none',alpha=0.12,zorder=1))
ax.text(0.50,0.9875,f'Cohort-wide  (all {N_TOT} samples pooled)',ha='center',va='center',fontsize=15,fontweight='bold',color='#0f2747')
bridge_drawn=False; has_incoming=set()
for a_,b_ in CANON:
    if a_ in faded or b_ in faded: continue
    has_incoming.add(b_)
    pr=pooled(a_,b_); w=max(cf(a_,b_),0.12); db=False
    if (a_,b_)==('Luminal-Distal','Basal-like'): w=max(w,0.30)
    if pr is None or pr[1]=='headless':
        edge(ax,a_,b_,min(w,1),'#9a9a9a',(0,(4,3)),cfree(a_,b_),node_r[a_],node_r[b_],headless=True); continue
    direc,lvl=pr
    if lvl=='solid': col,ls='#111','-'
    elif lvl=='bidir': col,ls,db='#111','-',True
    else: col,ls='#7a86a0',(0,(4,3))
    edge(ax,direc[0],direc[1],min(w,1),col,ls,0.10,node_r[direc[0]],node_r[direc[1]],double=db)
    if (direc[0] not in MES and direc[1] in MES) or (direc[1] not in MES and direc[0] in MES): bridge_drawn=True
present_mes=[s for s in MES if s not in faded]
if not bridge_drawn and present_mes:
    target=max(present_mes,key=lambda s: len(cells[s])); up='Luminal-Distal'
    for c in ['Basal-like','Luminal-Persistent','Partial-EMT','Luminal-Distal']:
        if c not in faded and c not in MES: up=c;break
    edge(ax,up,target,0.20,'#7a86a0',(0,(4,3)),cfree(up,target),node_r[up],node_r[target]); has_incoming.add(target)
for s in STATE_COLORS:
    if s in faded or s=='Luminal-Distal' or s in has_incoming: continue
    up='Luminal-Distal'
    for c in FALLBACK.get(s,['Luminal-Distal']):
        if c not in faded: up=c;break
    edge(ax,up,s,0.10,'#c0c0c0',(0,(4,3)),cfree(up,s),node_r[up],node_r[s],headless=True)
# right-side 'Trajectory' axis (matches Appendix Fig A1)
axx=0.985
ax.add_patch(FancyArrowPatch((axx,0.87),(axx,0.05),arrowstyle='-|>',color='black',lw=1.6,mutation_scale=12,zorder=4))
ax.text(axx+0.012,0.46,'Trajectory',ha='center',va='center',fontsize=10.5,color='black',fontweight='bold',rotation=270)
for s in STATE_COLORS: circle(ax,s,node_r[s],fd=(s in faded))
fig.suptitle('Npp53 cohort-wide consensus fate map — data-driven arrow direction (CARTA + PhyloVelo pseudotime)',fontsize=13.5,fontweight='bold',color='#0f2747',y=0.965)
class _HArrow(HandlerPatch):
    def create_artists(self, legend, oh, xd, yd, w, h, fs, trans):
        y=h/2.-yd
        p=FancyArrowPatch((-xd,y),(-xd+w,y),arrowstyle=oh._astyle,color=oh._acol,lw=oh._alw,linestyle=oh._als,mutation_scale=9)
        p.set_transform(trans); return [p]
def _legarrow(color,ls,astyle,lw):
    h=FancyArrowPatch((0,0),(1,0),arrowstyle=astyle,color=color,lw=lw,linestyle=ls,mutation_scale=9)
    h._astyle=astyle; h._acol=color; h._als=ls; h._alw=lw; return h
leg=[_legarrow('#111','-','-|>',2.6),
     _legarrow('#111','-','<|-|>',2.6),
     _legarrow('#7a86a0',(0,(4,3)),'-|>',2.2),
     _legarrow('#9a9a9a',(0,(4,3)),'-',1.8)]
labels=['confident direction — single arrow (pooled MWU p<0.01, effect≥0.10, tumors agree)',
        'strong link, direction unresolved — double arrow ↔ (tumors disagree)',
        'weak lean — dashed single arrow (p<0.05) / guaranteed mesenchymal bridge',
        'headless dashed - - - (below threshold / unresolved / connectivity; no arrowhead)']
ax.legend(handles=leg,labels=labels,handler_map={FancyArrowPatch:_HArrow()},loc='lower center',bbox_to_anchor=(0.5,-0.10),ncol=2,fontsize=8.5,frameon=False,handlelength=3.4,handletextpad=0.7,columnspacing=2.4)
fig.text(0.5,0.005,'Circle size ∝ CARTA terminal-progenitor frequency (cohort) · number in box = state index (c1–c12) · faded = <20 cells in every sample.',ha='center',fontsize=8.5,color='#333',style='italic')
plt.subplots_adjust(top=0.92,bottom=0.10,left=0.02,right=0.98)
for fmt in ['png','pdf']:
    fig.savefig(f'{OUTDIR}/npp53_fate_map_datadriven_cohort.{fmt}',dpi=200,bbox_inches='tight',facecolor='white')
print(f'saved: {OUTDIR}/npp53_fate_map_datadriven_cohort.png  (N_TOT={N_TOT}, faded={sorted(SH[s] for s in faded)})')

import json as _j; open(f'{OUTDIR}/fig10_3_6_edges_rebuilt.json','w').write(_j.dumps(EDGES,indent=1))

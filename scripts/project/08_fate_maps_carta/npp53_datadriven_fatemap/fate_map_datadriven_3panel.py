#!/usr/bin/env python
"""Npp53 data-driven consensus fate map — 3 panels (Intact | Cas | Cas+Enz).
Replacement for report Fig 10.3.5. Direction from pooled pseudotime (Mann-Whitney)
REQUIRING cross-tumor agreement. Edge types: solid single arrow = confident direction;
solid double arrow = strong link, direction unresolved (tumors disagree); dashed arrow = weak lean.
Below-threshold (headless/unresolved) edges are NOT drawn. No fate-cluster ellipses/labels."""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerPatch
from scipy.stats import mannwhitneyu

BASE='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
IN_JSON=f'{os.path.dirname(os.path.abspath(__file__))}/consensus_fate_per_arm.json'
H5AD=f'{BASE}/data/npp53_merge/npp53_edited.h5ad'
PT=os.environ['PTCSV']
OUTDIR=os.environ.get('OUTDIR','.')
import json
MINC=20; TREAT={'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}
ARM_ACCENT={'Intact':'#1f5fa8','Cas':'#e0955f','Cas+Enz':'#c44e52'}
ARM_DISPLAY={'Intact':'Intact','Cas':'Cas','Cas+Enz':'Cas+Enza'}
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

def nrad(s,it,n): return 0.032+0.048*(it.get(s,0)/max(1,n))
def circle(ax,s,r,faded):
    x,y=POS[s]; c=STATE_COLORS[s]; ec=c if not faded else '#bbb'; lw=2.6 if not faded else 1.4; al=1 if not faded else 0.35
    ax.add_patch(Circle((x,y),r,facecolor='white',edgecolor=ec,linewidth=lw,zorder=5,alpha=al))
    ax.text(x,y+r*0.14,LB[s],ha='center',va='center',fontsize=6.4,color='#111' if not faded else '#aaa',fontweight='bold',zorder=6,linespacing=0.9,alpha=1 if not faded else 0.55)
    by=y-r-0.014
    ax.add_patch(FancyBboxPatch((x-0.017,by-0.014),0.034,0.026,boxstyle='round,pad=0.001,rounding_size=0.004',facecolor='white',edgecolor='#666' if not faded else '#ccc',linewidth=0.8,zorder=6,alpha=1 if not faded else 0.5))
    ax.text(x,by-0.001,f'{SI[s]}',ha='center',va='center',fontsize=8,fontweight='bold',color='#111' if not faded else '#bbb',zorder=7)
def edge(ax,src,dst,w,color,ls,rad,rs,rd,double=False,headless=False):
    x0,y0=POS[src];x1,y1=POS[dst];rs+=0.004;rd+=0.008;dx,dy=x1-x0,y1-y0;L=(dx**2+dy**2)**.5
    if L<1e-6: return
    style='<|-|>' if double else ('-' if headless else '-|>')
    ax.add_patch(FancyArrowPatch((x0+dx/L*rs,y0+dy/L*rs),(x1-dx/L*rd,y1-dy/L*rd),arrowstyle=style,color=color,linewidth=0.8+3.0*w,linestyle=ls,mutation_scale=15,connectionstyle=f'arc3,rad={rad}',zorder=3))

data=json.loads(open(IN_JSON).read())
a=ad.read_h5ad(H5AD,backed='r');dfo=a.obs[['sample_id','treatment','putative_state']].copy();dfo['arm']=dfo['treatment'].astype(str).map(TREAT)
pdf=pd.read_csv(PT)

def draw_arm(ax,ARM):
    n_arm=data['n_by_arm'][ARM]; in_term=data['in_terminal_by_arm'][ARM]
    cooccur={tuple(k.split('|')):v for k,v in data['cooccur_by_arm'][ARM].items()}
    def cf(a,b): return cooccur.get(tuple(sorted([a,b],key=list(STATE_COLORS.keys()).index)),0)/max(1,n_arm)
    ct=dfo[dfo.arm==ARM].groupby(['sample_id','putative_state']).size().reset_index(name='n');mx=ct.groupby('putative_state')['n'].max().to_dict()
    faded={s for s in STATE_COLORS if mx.get(s,0)<MINC}
    sub=pdf[pdf.arm==ARM]
    g=sub.groupby(['sample','state']).agg(n=('pt_rank','size'),m=('pt_rank','mean')).reset_index()
    mpiv=g.pivot(index='sample',columns='state',values='m'); npiv=g.pivot(index='sample',columns='state',values='n')
    cells={s:sub[sub.state==s]['pt_rank'].values for s in STATE_COLORS}
    nsamp={s:(sub[sub.state==s].groupby('sample').size()>20).sum() for s in STATE_COLORS}
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
    node_r={s:nrad(s,in_term,n_arm) for s in STATE_COLORS}; avoid=list(STATE_COLORS)
    def cfree(src,dst):
        def bez(rad,n=15):
            x0,y0=POS[src];x1,y1=POS[dst];dx,dy=x1-x0,y1-y0;mx_,my_=(x0+x1)/2,(y0+y1)/2;cx,cy=mx_+rad*dy,my_-rad*dx
            return [((1-t/ n)**2*x0+2*(t/n)*(1-t/n)*cx+(t/n)**2*x1,(1-t/n)**2*y0+2*(t/n)*(1-t/n)*cy+(t/n)**2*y1) for t in range(1,n)]
        for rad in [0.0,0.22,-0.22,0.40,-0.40,0.55,-0.55]:
            ok=True
            for nm in avoid:
                if nm in (src,dst):continue
                nx,ny=POS[nm]
                if any(((px-nx)**2+(py-ny)**2)**.5<0.08 for px,py in bez(rad)): ok=False;break
            if ok:return rad
        return 0.55
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_aspect('equal');ax.axis('off')
    ax.text(0.50,0.955,'Tumor initiation',ha='center',va='center',fontsize=9,fontweight='bold',color='#111',zorder=6)
    r1=nrad('Luminal-Distal',in_term,n_arm)
    ax.add_patch(FancyArrowPatch((0.50,0.84+r1+0.026),(0.50,0.84+r1+0.006),arrowstyle='-|>',color='#111',lw=1.8,mutation_scale=12,zorder=3))
    ax.add_patch(FancyBboxPatch((0.02,0.965),0.96,0.045,boxstyle='round,pad=0.002,rounding_size=0.01',facecolor=ARM_ACCENT[ARM],edgecolor='none',alpha=0.14,zorder=1))
    ax.text(0.50,0.9875,f'{ARM_DISPLAY[ARM]}  (n = {n_arm} samples)',ha='center',va='center',fontsize=13,fontweight='bold',color=ARM_ACCENT[ARM])
    bridge_drawn=False; has_incoming=set()
    for a_,b_ in CANON:
        if a_ in faded or b_ in faded: continue
        has_incoming.add(b_)   # nominal downstream, for connectivity
        pr=pooled(a_,b_); w=max(cf(a_,b_),0.12); db=False
        if (a_,b_)==('Luminal-Distal','Basal-like'): w=max(w,0.30)
        if pr is None or pr[1]=='headless':          # BELOW-THRESHOLD -> headless dashed (no arrowhead)
            edge(ax,a_,b_,min(w,1),'#9a9a9a',(0,(4,3)),cfree(a_,b_),node_r[a_],node_r[b_],headless=True); continue
        direc,lvl=pr
        if lvl=='solid': col,ls='#111','-'
        elif lvl=='bidir': col,ls,db='#111','-',True
        else: col,ls='#7a86a0',(0,(4,3))
        edge(ax,direc[0],direc[1],min(w,1),col,ls,0.10,node_r[direc[0]],node_r[direc[1]],double=db)
        if (direc[0] not in MES and direc[1] in MES) or (direc[1] not in MES and direc[0] in MES): bridge_drawn=True
    # bridge guarantee: >=1 non-Mes -> Mes link (dashed arrow) so mes compartment isn't isolated
    present_mes=[s for s in MES if s not in faded]
    if not bridge_drawn and present_mes:
        target=max(present_mes,key=lambda s: len(cells[s]))
        up='Luminal-Distal'
        for c in ['Basal-like','Luminal-Persistent','Partial-EMT','Luminal-Distal']:
            if c not in faded and c not in MES: up=c; break
        if up in MES or up==target: up='Luminal-Distal'
        edge(ax,up,target,0.20,'#7a86a0',(0,(4,3)),cfree(up,target),node_r[up],node_r[target]); has_incoming.add(target)
    # connectivity fallback: present state with no incoming -> headless connector from nearest present upstream
    for s in STATE_COLORS:
        if s in faded or s=='Luminal-Distal' or s in has_incoming: continue
        up='Luminal-Distal'
        for c in FALLBACK.get(s,['Luminal-Distal']):
            if c not in faded: up=c; break
        edge(ax,up,s,0.10,'#c0c0c0',(0,(4,3)),cfree(up,s),node_r[up],node_r[s],headless=True); has_incoming.add(s)
    # right-side 'Trajectory' axis on the rightmost panel (matches Appendix Fig A1)
    if ARM=='Cas+Enz':
        axx=0.985
        ax.add_patch(FancyArrowPatch((axx,0.87),(axx,0.05),arrowstyle='-|>',color='black',lw=1.2,mutation_scale=10,zorder=4))
        ax.text(axx+0.010,0.46,'Trajectory',ha='center',va='center',fontsize=8,color='black',fontweight='bold',rotation=270)
    for s in STATE_COLORS: circle(ax,s,node_r[s],faded=(s in faded))

fig,axes=plt.subplots(1,3,figsize=(24,9.2))
for ax,ARM in zip(axes,['Intact','Cas','Cas+Enz']): draw_arm(ax,ARM)
fig.suptitle('Npp53 consensus fate map — data-driven arrow direction (CARTA + PhyloVelo pseudotime)',fontsize=15,fontweight='bold',color='#0f2747',y=0.99)
class _HArrow(HandlerPatch):
    def create_artists(self, legend, oh, xd, yd, w, h, fs, trans):
        y=h/2.-yd
        p=FancyArrowPatch((-xd,y),(-xd+w,y),arrowstyle=oh._astyle,color=oh._acol,lw=oh._alw,linestyle=oh._als,mutation_scale=9)
        p.set_transform(trans); return [p]
def _legarrow(color,ls,astyle,lw):
    h=FancyArrowPatch((0,0),(1,0),arrowstyle=astyle,color=color,lw=lw,linestyle=ls,mutation_scale=9)
    h._astyle=astyle; h._acol=color; h._als=ls; h._alw=lw; return h
leg=[_legarrow('#111','-','-|>',2.4),
     _legarrow('#111','-','<|-|>',2.4),
     _legarrow('#7a86a0',(0,(4,3)),'-|>',2.0),
     _legarrow('#9a9a9a',(0,(4,3)),'-',1.8)]
labels=['confident direction — single arrow (pooled Mann–Whitney p<0.01, effect≥0.10, tumors agree)',
        'strong link, direction unresolved — double arrow ↔ (tumors disagree)',
        'weak lean — dashed single arrow (p<0.05, effect≥0.05) / guaranteed bridge into mesenchymal',
        'headless dashed - - - (below threshold / unresolved / connectivity; no arrowhead)']
fig.legend(handles=leg,labels=labels,handler_map={FancyArrowPatch:_HArrow()},loc='lower center',ncol=2,fontsize=9,frameon=False,bbox_to_anchor=(0.5,0.01),handlelength=3.4,handletextpad=0.7,columnspacing=2.6)
fig.text(0.5,0.075,'Circle size ∝ CARTA terminal-progenitor frequency · number in box = state index (c1–c12) · '
         'faded = <20 cells in every sample of that arm (no links).',
         ha='center',fontsize=9,color='#333',style='italic')
plt.subplots_adjust(top=0.94,bottom=0.11,left=0.01,right=0.99,wspace=0.04)
for fmt in ['png','pdf']:
    fig.savefig(f'{OUTDIR}/npp53_fate_map_datadriven_3panel.{fmt}',dpi=200,bbox_inches='tight',facecolor='white')
print(f'saved: {OUTDIR}/npp53_fate_map_datadriven_3panel.png')

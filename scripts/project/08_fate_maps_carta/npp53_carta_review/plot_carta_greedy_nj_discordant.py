#!/usr/bin/env python
"""Fig 10.3.4 — CARTA DAG stability: Greedy vs NJ at same k*, for the discordant
samples. Reuses draw_carta_map() from plot_carta_dag_per_sample.py (which already
suppresses cell states with <=20 cells per sample, like Fig 10.3.2)."""
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import ast, pandas as pd
M='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
CF=Path(f'{M}/analysis/evolution_features_v2/carta_fine12')
SRC=CF/'biology_review/plot_carta_dag_per_sample.py'
src=SRC.read_text(); cut=src.index('\nintact = sorted(VERD')   # exec only the defs, not the main loop
ns={}; exec(compile(src[:cut], str(SRC), 'exec'), ns)
draw_carta_map=ns['draw_carta_map']; VERD=ns['VERD']; K=ns['K_FROM_CONSENSUS']
RES_G=CF/'results_fine12_greedy'; RES_N=CF/'results_fine12_nj'
def rp(resdir,s,k):
    f=resdir/s/f'{s}_progenitors_{k}.txt'
    if not f.exists(): return []
    o=[]
    for ln in f.read_text().splitlines():
        ln=ln.strip()
        if ln:
            try: o.append(frozenset(ast.literal_eval(ln)))
            except Exception: pass
    return o
tsv=pd.read_csv(CF/'biology_review/cross_solver_stability.tsv',sep='\t')
import json as _j; disc_set=set(_j.load(open(CF/'biology_review/discordant_at_adopted_k_gt20.json'))['discordant'])
order=sorted(disc_set)
disc=[s for s in order if s in disc_set]
print("discordant samples:",disc)
fig=plt.figure(figsize=(20,13.6)); gs=GridSpec(3,4,figure=fig,hspace=0.55,wspace=0.15)
for i,s in enumerate(disc):
    k=int(K[s]); rr=VERD[VERD['sample']==s].iloc[0]; arm=rr['arm']; nl=int(rr['n_leaves_tree'])
    for j,(solver,rd) in enumerate([('Greedy',RES_G),('NJ',RES_N)]):
        ax=fig.add_subplot(gs[i//2,(i%2)*2+j])
        draw_carta_map(ax, s, arm, k, rp(rd,s,k), 0, nl, 0)
        ax.text(0.5,1.05,solver,transform=ax.transAxes,ha='center',va='bottom',fontsize=12.5,fontweight='bold',color='#0f2747')
fig.suptitle('CARTA DAG stability: Greedy vs NJ at same k*  ·  5 discordant samples  ·  >20-cell states only',
             fontsize=15,fontweight='bold',y=0.985)
fig.subplots_adjust(top=0.90,bottom=0.02,left=0.02,right=0.99)
for fmt in ['png','pdf']:
    fig.savefig(CF/f'biology_review/carta_greedy_vs_nj_discordant.{fmt}',dpi=170,bbox_inches='tight',facecolor='white')
print("saved carta_greedy_vs_nj_discordant.png")

#!/usr/bin/env python
"""Build combined per-cell pseudotime for all samples: best-run-per-sample by unique-value count.
Sources: v8 CSV (umap), old scaled-newick (root->tip), JZ201 SWEEPBEST (if rescued).
Writes combined_pseudotime_percell.csv (sample,cellBC,pt_raw,pt_rank,run) + per-state ordering."""
import os, glob, numpy as np, pandas as pd
from ete3 import Tree
from scipy.stats import rankdata
ROOT='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
OLD=f'{ROOT}/lineage_tracing/cassiopeia_phylovelo_out'
V8=f'{ROOT}/lineage_tracing/cassiopeia_trees_v8_rebuilt/phylovelo_out'
LBL=f'{ROOT}/analysis/evolution_features_v2/carta_fine12/inputs_fine12'
SCRATCH=os.environ['SWEEP_OUT']
scfp=pd.read_csv(f'{ROOT}/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv')
am=scfp.groupby('sample')['treatment'].first().to_dict()
ARM={'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}
SO=['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR','Luminal-Persistent','Luminal-Squamous',
    'Basal-like','Partial-EMT','Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
SH={s:f'c{i+1}' for i,s in enumerate(SO)}; NUNIQ_MIN=20

def v8pt(s):
    f=f'{V8}/{s}_phylovelo_pseudotime.csv'
    return pd.read_csv(f,index_col=0)['pseudotime'] if os.path.exists(f) else None
def oldpt(s):
    f=f'{OLD}/{s}_pseudotime_scaled.newick'
    if not os.path.exists(f): return None
    t=Tree(f,format=1); return pd.Series({l.name:t.get_distance(l) for l in t.get_leaves()})
def sweepbest(s):
    g=sorted(glob.glob(f'{SCRATCH}/{s}_pseudotime_SWEEPBEST_*.csv'))
    return pd.read_csv(g[-1],index_col=0)['pseudotime'] if g else None

recs=[]; chosen={}
for s in sorted(am):
    cands=[('v8',v8pt(s)),('old',oldpt(s)),('sweep',sweepbest(s))]
    valid=[(tag,pt) for tag,pt in cands if pt is not None and pt.round(4).nunique()>=NUNIQ_MIN]
    if not valid: chosen[s]='DEAD'; continue
    tag,pt=max(valid,key=lambda x:x[1].round(4).nunique()); chosen[s]=tag
    lf=f'{LBL}/{s}_labels_fine12.tsv'
    if not os.path.exists(lf): continue
    lab=pd.read_csv(lf,sep='\t').set_index('cellBC')['cell_state']
    j=pd.DataFrame({'pt_raw':pt,'state':lab}).dropna()
    j['pt_rank']=rankdata(j['pt_raw'])/len(j); j['sample']=s; j['arm']=ARM.get(am.get(s)); j['run']=tag
    recs.append(j.reset_index().rename(columns={'index':'cellBC'}))
df=pd.concat(recs)
df.to_csv(f'{SCRATCH}/combined_pseudotime_percell.csv',index=False)
print('chosen run per sample:',chosen)
print('DEAD:',[s for s,v in chosen.items() if v=='DEAD'])
print(f'\nwrote combined_pseudotime_percell.csv: {df["sample"].nunique()} samples, {len(df)} cells')
def order(g): return ' < '.join(f"{SH[s]}({v:.2f})" for s,v in g.sort_values().items())
print('\nSHARED:',order(df.groupby('state')['pt_rank'].mean()))
for a in ['Intact','Cas','Cas+Enz']:
    sub=df[df['arm']==a]; print(f'{a:8s}(n={sub["sample"].nunique()}):',order(sub.groupby('state')['pt_rank'].mean()))

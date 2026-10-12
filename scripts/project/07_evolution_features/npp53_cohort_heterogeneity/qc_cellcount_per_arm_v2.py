#!/usr/bin/env python
"""v2 of QC: linear-scale Panel B + add cell-state SAMPLE-SPECIFICITY classification.

New question (Zejian 2026-09-02): which states are SHARED across most samples vs
SAMPLE-SPECIFIC (dominated by few samples)? Motivates dropping arm-level mean comparisons.
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

M = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
SCFP = M/'1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT = M/'1-mouse_based/Npp53_tracing/analysis/evolution_features_v2'

STATES_DISPLAY = {
    'Luminal-Distal':'Lum-Distal','Luminal-Proximal':'Lum-Proximal','Luminal-ARdep':'Lum-ARdep',
    'Luminal-CR':'Lum-CR','Luminal-Persistent':'Lum-Persistent','Luminal-Squamous':'Lum-Squamous',
    'Basal-like':'Basal-like','Partial-EMT':'Partial-EMT','Mes-Quiescent':'Mes-Quiescent',
    'Mes-Collagen':'Mes-Collagen','Mes-Stemness':'Mes-Stemness','Mes-myCAF':'Mes-myCAF',
}
STATE_ORDER = list(STATES_DISPLAY.keys())
STATE_ORDER_DISPLAY = list(STATES_DISPLAY.values())
STATE_COLORS = {
    'Lum-Distal':'#26a69a','Lum-Proximal':'#ec407a','Lum-ARdep':'#7cb342','Lum-CR':'#e53935',
    'Lum-Persistent':'#5c6bc0','Lum-Squamous':'#8d6e63','Basal-like':'#4a7fbf','Partial-EMT':'#ffca28',
    'Mes-Quiescent':'#bcaaa4','Mes-Collagen':'#78909c','Mes-Stemness':'#ce93d8','Mes-myCAF':'#a1887f',
}
ARM_ORDER = ['Intact', 'Castration', 'Castration + Enzalutamide']
ARM_SHORT = {'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}
ARM_COLORS = {'Intact':'#4a7fbf','Castration':'#e0955f','Castration + Enzalutamide':'#c44e52'}

# Rebuild count matrix
scfp = pd.read_csv(SCFP)
sample_arm = scfp.groupby('sample')['treatment'].first().reset_index()
sample_arm['arm_order'] = sample_arm['treatment'].map({a:i for i,a in enumerate(ARM_ORDER)})
sample_arm = sample_arm.sort_values(['arm_order','sample']).reset_index(drop=True)
long = scfp.groupby(['sample','treatment','putative_state']).size().reset_index(name='n_cells')
mat = long.pivot_table(index='sample', columns='putative_state', values='n_cells', aggfunc='sum').fillna(0).astype(int)
mat = mat.reindex(sample_arm['sample'].tolist())
mat = mat.reindex(columns=STATE_ORDER, fill_value=0)
sample_to_arm = sample_arm.set_index('sample')['treatment'].to_dict()

# ------- Sample-specificity analysis -------
# For each state:
#   n_samples_present_thresh = # samples with >= N cells (N=1, 5, 10, 30)
#   Gini / entropy of cell distribution across samples (concentration)
#   top-sample-fraction (max cells in one sample / total)
def gini(x):
    x = np.asarray(x, dtype=float)
    if x.sum() == 0: return 0
    x = np.sort(x)
    n = len(x)
    return (2 * np.sum((np.arange(1, n+1)) * x) - (n+1) * x.sum()) / (n * x.sum())

univ_rows = []
for st in STATE_ORDER:
    col = mat[st].values
    tot = col.sum()
    n_samp = len(col)
    n_ge1  = int((col >= 1).sum())
    n_ge5  = int((col >= 5).sum())
    n_gt20 = int((col > 20).sum())
    n_ge30 = int((col >= 30).sum())
    max_frac = float(col.max() / tot) if tot > 0 else 0.0
    top3_frac = float(np.sort(col)[::-1][:3].sum() / tot) if tot > 0 else 0.0
    g = float(gini(col))
    # Classification (Zejian's request):
    if n_gt20 >= 20:      cls = 'ubiquitous'         # ≥20/25 samples >20 cells
    elif n_gt20 >= 12:    cls = 'shared_majority'    # ≥12/25 samples >20 cells (majority)
    elif n_gt20 >= 5:     cls = 'multi_sample'       # 5-11 samples
    elif n_gt20 >= 2:     cls = 'sample_specific'    # 2-4 samples
    else:                 cls = 'single_sample'      # ≤1 sample
    univ_rows.append(dict(state=st, state_display=STATES_DISPLAY[st],
                          total_cells=int(tot), n_samples_ge1=n_ge1, n_samples_ge5=n_ge5,
                          n_samples_gt20=n_gt20, n_samples_ge30=n_ge30,
                          max_sample_frac=max_frac, top3_sample_frac=top3_frac,
                          gini_across_samples=g, classification=cls))
univ = pd.DataFrame(univ_rows)
univ.to_csv(OUT/'cellstate_universality_v10.csv', index=False)
print(f'saved: cellstate_universality_v10.csv')
print(univ[['state_display','total_cells','n_samples_gt20','max_sample_frac','gini_across_samples','classification']].to_string(index=False))

# ------- FIGURE (v2): linear Panel B + Panel C = universality ranking -------
fig = plt.figure(figsize=(15, 13))
gs = gridspec.GridSpec(3, 1, height_ratios=[2.1, 0.9, 1.1], hspace=0.55, top=0.955, bottom=0.05, left=0.07, right=0.97)

# Panel A: per-sample × per-state dot plot (unchanged from v1)
axA = fig.add_subplot(gs[0, 0])
for i, s in enumerate(mat.index):
    for j, st in enumerate(STATE_ORDER):
        n = mat.loc[s, st]
        if n == 0: continue
        size = min(500, 20 * np.log10(n+1)**2)
        axA.scatter(j, i, s=size, c=STATE_COLORS[STATES_DISPLAY[st]], edgecolors='black', linewidth=0.5, alpha=0.85)
        axA.text(j, i, str(n), ha='center', va='center', fontsize=7, color='white' if n>50 else 'black',
                 fontweight='bold' if n<10 else 'normal')
axA.set_xticks(np.arange(len(STATE_ORDER))); axA.set_xticklabels(STATE_ORDER_DISPLAY, rotation=45, ha='right', fontsize=10)
sample_labels = [f'{s} ({ARM_SHORT[sample_to_arm[s]]})' for s in mat.index]
axA.set_yticks(np.arange(len(mat.index))); axA.set_yticklabels(sample_labels, fontsize=9)
axA.set_xlim(-0.7, len(STATE_ORDER)-0.3); axA.set_ylim(-0.7, len(mat.index)-0.3)
axA.set_xlabel('Cell state', fontsize=11); axA.set_ylabel('Sample (arm)', fontsize=11)
axA.set_title('A · Per-sample × per-state cell counts (v10 sc_fp; Cycling + Ionocyte excluded)',
              fontsize=12, fontweight='bold', color='#0f2747', loc='left')
axA.grid(True, alpha=0.2, linestyle=':')
# arm-band background shading
prev = None; boundaries = [0]
for i, s in enumerate(mat.index):
    a = sample_to_arm[s]
    if prev and a != prev: boundaries.append(i)
    prev = a
boundaries.append(len(mat.index))
for start, end in zip(boundaries[:-1], boundaries[1:]):
    arm = sample_to_arm[mat.index[start]]
    axA.axhspan(start-0.5, end-0.5, alpha=0.05, color=ARM_COLORS[arm])
axA.invert_yaxis()

# Panel B: LINEAR-scale mean cells per sample (Zejian: no log)
axB = fig.add_subplot(gs[1, 0])
x = np.arange(len(STATE_ORDER)); width = 0.27
for i, arm in enumerate(ARM_ORDER):
    samples_in_arm = sample_arm[sample_arm['treatment']==arm]['sample'].tolist()
    sub = mat.loc[samples_in_arm]
    means = [sub[s].mean() for s in STATE_ORDER]
    mins  = [sub[s].min()  for s in STATE_ORDER]
    maxs  = [sub[s].max()  for s in STATE_ORDER]
    yerr = np.array([[m-mn for m,mn in zip(means,mins)], [mx-m for m,mx in zip(means,maxs)]])
    axB.bar(x + (i-1)*width, means, width, yerr=yerr, capsize=3,
            color=ARM_COLORS[arm], alpha=0.75, edgecolor='black', linewidth=0.4,
            label=f'{ARM_SHORT[arm]} (n_samples={len(samples_in_arm)})')
axB.axhline(20, color='red', lw=0.7, linestyle='--', alpha=0.6)
axB.text(len(STATE_ORDER)-0.3, 25, '>20 cells threshold', color='red', fontsize=9, ha='right', va='bottom')
axB.set_xticks(x); axB.set_xticklabels(STATE_ORDER_DISPLAY, rotation=45, ha='right', fontsize=10)
axB.set_ylabel('Cells per sample (mean; whiskers = min-max) — LINEAR scale', fontsize=10)
axB.set_title('B · Per-state cell abundance by treatment arm — LINEAR y (mean ± range)',
              fontsize=12, fontweight='bold', color='#0f2747', loc='left')
axB.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.005,1.0), frameon=False)
axB.grid(axis='y', alpha=0.3)

# Panel C: universality — bar of n_samples_ge10 per state, colored by classification
CLASS_COLORS = {'ubiquitous':'#047857','shared_majority':'#22c55e','multi_sample':'#f59e0b',
                'sample_specific':'#c44e52','single_sample':'#7f1d1d'}
axC = fig.add_subplot(gs[2, 0])
univ_sorted = univ.sort_values('n_samples_gt20', ascending=False).reset_index(drop=True)
bars = axC.bar(range(len(univ_sorted)), univ_sorted['n_samples_gt20'],
               color=[CLASS_COLORS[c] for c in univ_sorted['classification']],
               edgecolor='black', linewidth=0.6)
for i, r in univ_sorted.iterrows():
    axC.text(i, r['n_samples_gt20'] + 0.5, f'{int(r["n_samples_gt20"])}/25', ha='center', fontsize=9)
axC.axhline(20, color='#047857', ls='--', lw=0.6, alpha=0.5)
axC.axhline(12, color='#22c55e', ls='--', lw=0.6, alpha=0.5)
axC.axhline(5,  color='#f59e0b', ls='--', lw=0.6, alpha=0.5)
axC.set_xticks(range(len(univ_sorted))); axC.set_xticklabels(univ_sorted['state_display'], rotation=45, ha='right', fontsize=10)
axC.set_ylabel('N samples with >20 cells (of 25)', fontsize=10)
axC.set_ylim(0, 27)
axC.set_title('C · Cell-state universality across 25 samples (bar height = n_samples >20 cells; color = class)',
              fontsize=11, fontweight='bold', color='#0f2747', loc='left')
# legend for classes
handles = [plt.Rectangle((0,0),1,1,color=c) for c in CLASS_COLORS.values()]
labels  = [f'{cls} ({(univ["classification"]==cls).sum()} states)' for cls in CLASS_COLORS]
axC.legend(handles, labels, fontsize=9, loc='upper left', bbox_to_anchor=(1.005,1.0), frameon=False)
axC.grid(axis='y', alpha=0.3)

fig.suptitle('Cell-state coverage QC — per-sample counts (A) + per-arm abundance linear (B) + universality across samples (C)',
             fontsize=13, fontweight='bold', y=0.995, color='#0f2747')
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'Fig_9_1_QC_cellcount.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'\n[saved] {OUT}/Fig_9_1_QC_cellcount.{{pdf,png}} (v2 — linear B + universality panel C)')

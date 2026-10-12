#!/usr/bin/env python
"""Per-arm variant of evopath_per_type_figs.py.

Same pipeline as evopath_MODULES / evopath_GENES (Fig 10.3.10-13) but split by
treatment arm into a 3-column layout (Intact | Cas | Cas+Enz).

Design choices held constant vs the cohort version:
  * pt_rank = "Corrected PhyloVelo pseudotime" from combined_pseudotime_percell.csv.
  * Bins = clip((pt_rank_within_arm * NB).astype(int), 0, NB-1) → equal-cell bins.
  * NB = 20 bins per arm.
  * Canonical c1..c12 composition-strip stack (bottom-to-top), same SCOL palette.
  * X-axis label: only "early (initiator)" and "late (terminal)" ticks + arrow.
  * Row order = SHARED across arms, computed once on the cohort as
    argsort(argmax_over_cohort_bins(row_z)) so that within-arm z-scored rows
    stack in the same order as Fig 10.3.13 (M5 → M4 → M1 → M3 → M2).

Emits evopath_MODULES_perarm.png/pdf. Copy to Fig 10.3.14 slot upstream."""
import warnings; warnings.filterwarnings('ignore')
import os, re, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

ROOT = '/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
AN   = f'{ROOT}/analysis'
HERE = f'{AN}/evolution_features_v2/carta_fine12/biology_review/evopath_trajectory_features'
PT   = f'{HERE}/combined_pseudotime_percell.csv'
OUT  = HERE
NB   = 20
ARMS = ['Intact', 'Cas', 'Cas+Enz']
ARM_ACCENT = {'Intact': '#1f5fa8', 'Cas': '#e0955f', 'Cas+Enz': '#c44e52'}

STATE_ORDER = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
               'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
               'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
SCOL = {'Luminal-Distal':'#1f5fa8','Luminal-Proximal':'#9467bd','Luminal-ARdep':'#d23030',
        'Luminal-CR':'#f1903d','Luminal-Persistent':'#117c70','Luminal-Squamous':'#f78ec7',
        'Basal-like':'#8c564b','Partial-EMT':'#bcbd22','Mes-Quiescent':'#c49c94',
        'Mes-Collagen':'#5b5b5b','Mes-Stemness':'#8c1a1a','Mes-myCAF':'#4daf4a'}
CID = {'Luminal-Distal':'c1','Luminal-Proximal':'c2','Luminal-ARdep':'c3','Luminal-CR':'c4',
       'Luminal-Persistent':'c5','Luminal-Squamous':'c6','Basal-like':'c7','Partial-EMT':'c8',
       'Mes-Quiescent':'c9','Mes-Collagen':'c10','Mes-Stemness':'c11','Mes-myCAF':'c12'}

# ── Load per-cell pseudotime + arm + h5ad ──────────────────────────────────
pt = pd.read_csv(PT)
pt['k'] = pt['sample'] + '|' + pt['cellBC']
k2pt   = dict(zip(pt['k'], pt['pt_rank']))
k2st   = dict(zip(pt['k'], pt['state']))
k2arm  = dict(zip(pt['k'], pt['arm']))

a = sc.read_h5ad(f'{ROOT}/data/npp53_merge/npp53_edited.h5ad')
a.obs['bc'] = [re.search(r'[ACGTN]{14,18}', x).group(0) for x in a.obs_names]
a.obs['k']  = a.obs['sample_id'].astype(str) + '|' + a.obs['bc']
a.obs['pt']  = a.obs['k'].map(k2pt)
a.obs['st']  = a.obs['k'].map(k2st)
a.obs['arm'] = a.obs['k'].map(k2arm)
a = a[a.obs['pt'].notna() & a.obs['arm'].notna()].copy()
print(f'  {a.shape[0]:,} cells with pt_rank + arm')

# ── FSU module scores (score_genes, use_raw=False) ─────────────────────────
fsu_df = pd.read_csv(f'{AN}/hotspot_canonical/cohort_genes_by_module_Intact.csv')
present = set(a.var_names)
mods = sorted(fsu_df.Module.unique())
for m in mods:
    gg = [x for x in fsu_df[fsu_df.Module == m]['Gene'] if x in present]
    sc.tl.score_genes(a, gg, score_name=f'FSU:M{m}', use_raw=False)
    print(f'  M{m}: {len(gg)} genes scored')

# Row order = canonical M1..M5 (matches Fig 10.3.13 which also uses M1..M5 now).
DESIRED_ORDER_TOP_TO_BOTTOM = [1, 2, 3, 4, 5]
row_order = np.array([mods.index(m) for m in DESIRED_ORDER_TOP_TO_BOTTOM])
mods_ordered = DESIRED_ORDER_TOP_TO_BOTTOM
labs_ordered = [f'M{m}' for m in mods_ordered]
print(f'  row order (canonical M1..M5): {labs_ordered}')

# ── Per-arm heatmap + composition ──────────────────────────────────────────
def per_arm_matrices(arm_name):
    mask = (a.obs['arm'] == arm_name).values
    if mask.sum() == 0:
        return None, None, 0
    ptv = a.obs['pt'].values[mask].astype(float)
    stv = a.obs['st'].values[mask]
    # Re-rank within arm so bins are equal-cell for THIS arm
    ranks_within = pd.Series(ptv).rank(method='average').values
    pt_rank_arm = (ranks_within - 1) / (len(ranks_within) - 1) if len(ranks_within) > 1 else np.zeros_like(ptv)
    bins = np.clip((pt_rank_arm * NB).astype(int), 0, NB - 1)
    # Module × bin (mean score per bin), row-z within arm
    M = np.full((len(mods), NB), np.nan)
    for r, m in enumerate(mods):
        v = a.obs[f'FSU:M{m}'].values[mask].astype(float)
        for b in range(NB):
            mm = bins == b
            if mm.sum() >= 5:
                M[r, b] = np.nanmean(v[mm])
    Mz = (M - np.nanmean(M, 1, keepdims=True)) / (np.nanstd(M, 1, keepdims=True) + 1e-9)
    Mz = Mz[row_order]
    # Composition per bin (canonical stack)
    comp = np.zeros((len(STATE_ORDER), NB))
    for b in range(NB):
        mm = bins == b
        if mm.sum():
            vc = pd.Series(stv[mm]).value_counts(normalize=True)
            for i, s in enumerate(STATE_ORDER):
                comp[i, b] = vc.get(s, 0)
    return Mz, comp, int(mask.sum())

arm_data = {arm: per_arm_matrices(arm) for arm in ARMS}
for arm, (_, _, n) in arm_data.items():
    print(f'  {arm}: n = {n:,} cells')

# ── Figure — 3-col per-arm (matches evopath draw() layout, one column per arm) ─
n_rows = len(labs_ordered)
cellh  = 0.55
TIN, GAP, SIN, BOT = 0.70, 0.14, 0.85, 1.55
heat_h = n_rows * cellh
H = TIN + heat_h + GAP + SIN + BOT
L, R = 0.04, 0.98
W = R - L
# Wider inter-column gap so per-arm "early (initiator)" / "late (terminal)"
# tick labels don't collide with the neighbouring arm's labels.
col_w   = W / len(ARMS) * 0.86
col_gap = W / len(ARMS) * 0.14

fig = plt.figure(figsize=(20, H))

im_last = None
for j, arm in enumerate(ARMS):
    Mz, comp, n_arm = arm_data[arm]
    left = L + j * (col_w + col_gap)
    axH = fig.add_axes([left, 1 - (TIN + heat_h) / H, col_w, heat_h / H])
    axC = fig.add_axes([left, BOT / H, col_w, SIN / H])
    if Mz is None:
        axH.set_axis_off(); axC.set_axis_off(); continue
    im = axH.imshow(Mz, aspect='auto', cmap='RdBu_r', vmin=-1.8, vmax=1.8)
    im_last = im
    axH.set_xticks([])
    if j == 0:
        axH.set_yticks(range(n_rows)); axH.set_yticklabels(labs_ordered, fontsize=10)
    else:
        axH.set_yticks([])
    axH.set_title(f'{arm}   (n = {n_arm:,} cells, {NB} pt-rank bins)',
                  fontsize=12, fontweight='bold', color=ARM_ACCENT[arm], pad=6)

    # Composition strip
    bottom = np.zeros(NB)
    for i, s in enumerate(STATE_ORDER):
        axC.bar(range(NB), comp[i], bottom=bottom, width=1.0,
                color=SCOL[s], edgecolor='none'); bottom += comp[i]
    axC.set_xlim(-0.5, NB - 0.5); axC.set_ylim(0, 1); axC.set_yticks([])
    axC.set_xticks([0, NB - 1])
    # Left-align the "early" label and right-align the "late" label at the ends
    # so they hug the panel edges instead of centering (which caused inter-arm collision).
    axC.set_xticklabels(['early\n(initiator)', 'late\n(terminal)'], fontsize=8)
    for tick, ha in zip(axC.get_xticklabels(), ['left', 'right']):
        tick.set_horizontalalignment(ha)
    axC.set_xlabel('Corrected PhyloVelo pseudotime  →', fontsize=10)

# Shared title
fig.text(L, 1 - (TIN * 0.35) / H,
         'FSU Hotspot modules following the tumor-evolution trajectory  ·  per treatment arm  '
         '(canonical M1..M5 row order, same as Fig 10.3.13)',
         fontsize=11.5, fontweight='bold', color='#059669', ha='left', va='center')

# Shared legend (canonical c1..c12) + colorbar
fig.legend(handles=[Patch(color=SCOL[s], label=f'{CID[s]}_{s}') for s in STATE_ORDER],
           loc='lower center', ncol=6, fontsize=8, frameon=False,
           bbox_to_anchor=(0.5, 0.60 / H))
if im_last is not None:
    cax = fig.add_axes([0.42, 0.25 / H, 0.16, 0.10 / H])
    plt.colorbar(im_last, cax=cax, orientation='horizontal').set_label(
        'row z-score (within-arm)', fontsize=8)

for fmt in ['png', 'pdf']:
    fig.savefig(f'{OUT}/evopath_MODULES_perarm.{fmt}',
                dpi=170, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'saved: {OUT}/evopath_MODULES_perarm.png')

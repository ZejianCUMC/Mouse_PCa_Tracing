#!/usr/bin/env python
"""§10 evolution — per-arm UMAP + diffmap (DC2 vs DC3) in BIG-PANEL style,
matching Fig 6.4b (figures/npp53_plots/make_sec6_4b_grid.py).

Two 1x3 figures (Intact | Cas | Cas+Enz), each 24 x 8.5 in — big enough to read.
Uses the SAME pre-computed embeddings that Fig 6.4b consumed:
  data/npp53_merge/umap_no_cycling_no_ionocyte.csv     (2-col UMAP)
  data/npp53_merge/diffmap_no_cycling_no_ionocyte.csv  (DC2 DC3 DC4 DC5)
Cells filtered to the 18,162-cell scfp v10 set (Ionocyte + Cycling excluded).
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

M   = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
UMAP_CSV = M/'data/npp53_merge/umap_no_cycling_no_ionocyte.csv'
DM_CSV   = M/'data/npp53_merge/diffmap_no_cycling_no_ionocyte.csv'
OBS_CSV  = M/'data/npp53_merge/obs_no_cycling_no_ionocyte.csv'
SCFP     = M/'analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT      = M/'analysis/evolution_features_v2'
OUT.mkdir(exist_ok=True)

# 12-state palette + display names — identical to Fig 6.4b
STATE_MAP = {
    'Luminal-Distal':'c1_Lum-Distal','Luminal-Proximal':'c2_Lum-Proximal','Luminal-ARdep':'c3_Lum-ARdep',
    'Luminal-CR':'c4_Lum-CR','Luminal-Persistent':'c5_Lum-Persistent','Luminal-Squamous':'c6_Lum-Squamous',
    'Basal-like':'c7_Basal-like','Partial-EMT':'c8_Partial-EMT','Mes-Quiescent':'c9_Mes-Quiescent',
    'Mes-Collagen':'c10_Mes-Collagen','Mes-Stemness':'c11_Mes-Stemness','Mes-myCAF':'c12_Mes-myCAF',
}
STATE_DISPLAY = list(STATE_MAP.values())
STATE_COLORS = {
    'c1_Lum-Distal':'#26a69a','c2_Lum-Proximal':'#ec407a','c3_Lum-ARdep':'#7cb342',
    'c4_Lum-CR':'#e53935','c5_Lum-Persistent':'#5c6bc0','c6_Lum-Squamous':'#8d6e63',  # c4 red — distinct from c8 yellow
    'c7_Basal-like':'#4a7fbf','c8_Partial-EMT':'#ffca28','c9_Mes-Quiescent':'#bcaaa4',
    'c10_Mes-Collagen':'#78909c','c11_Mes-Stemness':'#ce93d8','c12_Mes-myCAF':'#a1887f',
}
ARM_ORDER = ['Intact', 'Castration', 'Castration + Enzalutamide']
ARM_SHORT = {'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}

# --- load embeddings + obs -------------------------------------------------
print('[load] pre-computed embeddings + obs ...')
umap = pd.read_csv(UMAP_CSV, index_col=0)
dm   = pd.read_csv(DM_CSV,   index_col=0)
obs  = pd.read_csv(OBS_CSV,  index_col=0)
print(f'  umap: {umap.shape}   diffmap: {dm.shape}   obs: {obs.shape}')

# align obs & embeddings on the umap index (that's the canonical "no-cyc no-ion" set)
obs  = obs.loc[umap.index.intersection(obs.index)].copy()
umap = umap.loc[obs.index]; dm = dm.loc[obs.index]
print(f'  aligned to {len(obs):,} cells')

# --- filter to the 18,162 scfp v10 cells (Ion+Cyc excluded) ----------------
scfp = pd.read_csv(SCFP)
scfp['bc16']  = scfp['cellBC'].astype(str)
scfp['key']   = scfp['sample'].astype(str) + '::' + scfp['bc16']
obs['bc16']   = obs.index.str.split('_').str[1].str.replace('-1','', regex=False)
# obs already carries sample_id (JZ*) and treatment
obs['key']    = obs['sample_id'].astype(str) + '::' + obs['bc16']
mask          = obs['key'].isin(set(scfp['key']))
obs           = obs[mask].copy()
umap          = umap.loc[obs.index]; dm = dm.loc[obs.index]
print(f'  after scfp filter: {len(obs):,} cells')

# state display + treatment
obs['state_disp'] = obs['putative_state'].map(STATE_MAP)
obs = obs[obs['state_disp'].notna()].copy()
umap = umap.loc[obs.index]; dm = dm.loc[obs.index]
print(f'  final tumor cells: {len(obs):,}')
print(f'  per arm: {obs.treatment.value_counts().to_dict()}')

# --- plotting helper (big-panel style) -------------------------------------
def facet_big(coords, xcol_ycol_labels, title, fname):
    """coords: (N, 2) numpy array; xcol_ycol_labels: (xlab, ylab)."""
    xlab, ylab = xcol_ycol_labels
    fig, axes = plt.subplots(1, 3, figsize=(24, 8.5))
    all_x, all_y = coords[:,0], coords[:,1]
    pad_x = 0.03 * (all_x.max() - all_x.min())
    pad_y = 0.03 * (all_y.max() - all_y.min())
    xlim = (all_x.min()-pad_x, all_x.max()+pad_x)
    ylim = (all_y.min()-pad_y, all_y.max()+pad_y)
    for ax, arm in zip(axes, ARM_ORDER):
        # background: all cells in light grey (context)
        ax.scatter(all_x, all_y, c='#e2e8f0', s=2.5, linewidths=0, rasterized=True)
        # foreground: this arm's cells colored by state
        arm_mask = (obs['treatment'] == arm).values
        for st in STATE_DISPLAY:
            ms = arm_mask & (obs['state_disp'] == st).values
            if ms.sum() == 0: continue
            ax.scatter(coords[ms,0], coords[ms,1], c=STATE_COLORS[st],
                       s=2.5, alpha=0.9, edgecolors='none', rasterized=True)
        n_cells = int(arm_mask.sum()); n_mice = obs.loc[arm_mask,'sample_id'].nunique()
        ax.set_title(f'{ARM_SHORT[arm]}   ({n_mice} mice · {n_cells:,} cells)',
                     fontsize=18, fontweight='bold', color='#0f2747')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_aspect('equal', 'datalim')
        for s in ax.spines.values(): s.set_visible(False)
        ax.set_xlabel(xlab, fontsize=13); ax.set_ylabel(ylab, fontsize=13)
        # cell-count annotation, box style like Fig 6.4b
        ax.text(0.02, 0.02, f'n = {n_cells:,}', transform=ax.transAxes,
                fontsize=13, ha='left', va='bottom', color='#333',
                bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=3))
    fig.suptitle(title, fontsize=17, fontweight='bold', color='#0f2747', y=1.02)
    # reserve right whitespace for legend BEFORE placing it (figure coords)
    fig.tight_layout(rect=[0, 0, 0.86, 1.0])
    handles = [Patch(facecolor=STATE_COLORS[st], edgecolor='none', label=st) for st in STATE_DISPLAY]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.87, 0.5),
               fontsize=14, title='Tumor state', title_fontsize=15,
               frameon=False, borderpad=1, labelspacing=0.7)
    for fmt in ['pdf','png']:
        fig.savefig(OUT/f'{fname}.{fmt}', dpi=170, bbox_inches='tight')
    plt.close(fig)
    print(f'  [saved] {OUT/(fname+".png")}')

# --- Fig A · UMAP big-panel per arm ----------------------------------------
print('[fig A] UMAP big-panel per arm ...')
facet_big(umap[['umapsub_1','umapsub_2']].values,
          xcol_ycol_labels=('UMAP-1','UMAP-2'),
          title='§10 · scVI UMAP (No-Cycling + No-Ionocyte tumor cells) split by treatment arm',
          fname='UMAP_per_arm_BIG')

# --- Fig B · Diffmap DC2 vs DC3 big-panel per arm --------------------------
print('[fig B] Diffmap DC2 vs DC3 big-panel per arm ...')
facet_big(dm[['DC2','DC3']].values,
          xcol_ycol_labels=('DC2','DC3'),
          title='§10 · scVI Diffmap (DC2, DC3 — DC0 skipped as technical) split by treatment arm',
          fname='DIFFMAP_DC2vsDC3_per_arm_BIG')

print('done.')

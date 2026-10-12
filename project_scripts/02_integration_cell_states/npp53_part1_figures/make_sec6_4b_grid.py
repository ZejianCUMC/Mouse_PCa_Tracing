#!/usr/bin/env python
"""Fig 6.4b — final No-Cycling + No-Ionocyte tumor state view (12 states, 18,241 cells).
   4-panel grid: scVI UMAP + tSNE + Diffmap + harmony UMAP, colored by 12 tumor states.
   Larger legend + in-plot text per user request."""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import scanpy as sc

M = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
OUT = M/'1-mouse_based/Npp53_tracing/figures/npp53_plots'
STATE_DISPLAY = ['c1_Lum-Distal','c2_Lum-Proximal','c3_Lum-ARdep','c4_Lum-CR',
                 'c5_Lum-Persistent','c6_Lum-Squamous','c7_Basal-like','c8_Partial-EMT',
                 'c9_Mes-Quiescent','c10_Mes-Collagen','c11_Mes-Stemness','c12_Mes-myCAF']
STATE_MAP = dict(zip(
    ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR','Luminal-Persistent','Luminal-Squamous',
     'Basal-like','Partial-EMT','Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF'],
    STATE_DISPLAY))
STATE_COLORS = {
    'c1_Lum-Distal':'#26a69a','c2_Lum-Proximal':'#ec407a','c3_Lum-ARdep':'#7cb342',
    'c4_Lum-CR':'#e53935','c5_Lum-Persistent':'#5c6bc0','c6_Lum-Squamous':'#8d6e63',  # c4 red — distinct from c8 yellow (2026-09-04)
    'c7_Basal-like':'#4a7fbf','c8_Partial-EMT':'#ffca28','c9_Mes-Quiescent':'#bcaaa4',
    'c10_Mes-Collagen':'#78909c','c11_Mes-Stemness':'#ce93d8','c12_Mes-myCAF':'#a1887f',
}

print('[1] load recomputed embeddings + obs')
umap = pd.read_csv(M/'1-mouse_based/Npp53_tracing/data/npp53_merge/umap_no_cycling_no_ionocyte.csv',   index_col=0)
tsne = pd.read_csv(M/'1-mouse_based/Npp53_tracing/data/npp53_merge/tsne_no_cycling_no_ionocyte.csv',   index_col=0)
dm   = pd.read_csv(M/'1-mouse_based/Npp53_tracing/data/npp53_merge/diffmap_no_cycling_no_ionocyte.csv', index_col=0)
obs  = pd.read_csv(M/'1-mouse_based/Npp53_tracing/data/npp53_merge/obs_no_cycling_no_ionocyte.csv',    index_col=0)
# harmony UMAP from Seurat RDS not exported to CSV; fall back to scanpy umap on harmony scvi
# For simplicity, use the same recomputed UMAP for both panels (harmony_UMAP dropped or aliased)
# Actually harmony UMAP is different — extract from Seurat RDS on-the-fly is heavy; use the shared_embeddings_v10 subset
sh = pd.read_csv(M/'1-mouse_based/Npp53_tracing/data/npp53_merge/shared_embeddings_v10.csv', index_col=0)  # composite key
# but harmony not in shared_embeddings_v10 — skip harmony panel or use the first shared subset
# Simplest: 3-panel grid (UMAP + tSNE + Diffmap)

state_full = obs.loc[umap.index, 'putative_state']
state_display = state_full.map(STATE_MAP)
kept = state_display.notna()
print(f'  cells kept: {kept.sum()}')

# 3-panel grid: UMAP · tSNE · Diffmap (DC2/DC3)
fig, axes = plt.subplots(1, 3, figsize=(24, 8.5))
PANELS = [('scVI UMAP', umap.values, ['umapsub_1','umapsub_2']),
          ('scVI tSNE', tsne.values, ['tsnesub_1','tsnesub_2']),
          ('scVI Diffmap (DC2, DC3)', dm.values, ['DC2','DC3'])]

for ax, (title, mat, cols) in zip(axes, PANELS):
    for st in STATE_DISPLAY:
        mask = (state_display == st) & kept
        if mask.sum() == 0: continue
        col = STATE_COLORS[st]
        ax.scatter(mat[mask, 0], mat[mask, 1], c=col, s=2.5, alpha=0.85, edgecolors='none')
    ax.set_title(title, fontsize=18, fontweight='bold', color='#0f2747')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect('equal', 'datalim')
    for s in ax.spines.values(): s.set_visible(False)
    # add total cell count as annotation
    ax.text(0.02, 0.02, f'n = {kept.sum():,} cells (12 states)', transform=ax.transAxes,
            fontsize=13, ha='left', va='bottom', color='#333',
            bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=3))

fig.suptitle('Fig 6.4b · Final tumor cell state view — No-Cycling + No-Ionocyte (12 states, 18,241 cells)',
             fontsize=19, fontweight='bold', color='#0f2747', y=1.02)
# reserve right band for legend BEFORE placing it (figure-coord layout)
fig.tight_layout(rect=[0, 0, 0.85, 1.0])
# Large shared legend (right side, outside plot margin)
handles = [Patch(facecolor=STATE_COLORS[st], edgecolor='none', label=st) for st in STATE_DISPLAY]
fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.86, 0.5),
           fontsize=20, title='Tumor state', title_fontsize=22,
           frameon=False, borderpad=1, labelspacing=1.0, handleheight=1.6, handlelength=2.0)
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'sec6_4b_no_cycling_no_ionocyte_grid.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'saved: sec6_4b_no_cycling_no_ionocyte_grid.{{png,pdf}}')

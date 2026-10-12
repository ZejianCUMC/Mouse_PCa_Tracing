#!/usr/bin/env python
"""Rerender pySCENIC top-5-per-state clustered heatmap with fixed layout.

Input : top5_regulons_per_state_z_heatmap.csv (pre-computed 56 × 12 z-matrix)
Output: top5_regulons_per_state_heatmap.{png,pdf}   (in-place, no _v2 suffix)

Fixes vs v1:
  - figsize taller (11 × 15.5) so 56 y-labels breathe
  - y-tick rotation 0 (horizontal) — no more crooked labels
  - dendrogram_ratio (0.10, 0.07) — extra top space so title clears the top dendrogram
  - cbar moved to (0.020, 0.945, 0.030, 0.045) with 'AUCell z' label — no title/dendrogram collision
  - title placed via fig.text at (0.55, 0.995) centered, tight_layout not needed after cbar/dendro tuning
  - font sizes bumped: y=9 (readable), x=10.5
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

OUT = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/downstream')
htm = pd.read_csv(OUT/'top5_regulons_per_state_z_heatmap.csv', index_col=0)
print(f'[load] heatmap matrix: {htm.shape}  regulons × states')
print(f'       columns: {list(htm.columns)}')

# rename state labels to match the abbreviated convention used throughout the Part 1 report
# (c3_Luminal-ARdep in the raw sc_fitness_plasticity CSV → c3_Lum-ARdep in the report; etc.)
STATE_RENAME = {
    'c1_Luminal-Distal'     : 'c1_Lum-Distal',
    'c2_Luminal-Proximal'   : 'c2_Lum-Proximal',
    'c3_Luminal-ARdep'       : 'c3_Lum-ARdep',
    'c4_Luminal-CR'         : 'c4_Lum-CR',
    'c5_Luminal-Persistent' : 'c5_Lum-Persistent',
    'c6_Luminal-Squamous'   : 'c6_Lum-Squamous',
    # basal / partial / mes labels are already in report convention
}
htm = htm.rename(columns=STATE_RENAME)
print(f'[rename] columns → {list(htm.columns)}')

vmax = 3
htm_clip = htm.clip(lower=-vmax, upper=vmax)

g = sns.clustermap(
    htm_clip,
    cmap='RdBu_r', vmin=-vmax, vmax=vmax, center=0,
    figsize=(11.0, 15.5),                           # taller so 56 y-labels breathe
    method='average', metric='correlation',         # keep clustering approach
    row_cluster=True, col_cluster=True,
    dendrogram_ratio=(0.11, 0.075),                 # a touch more room for both dendros
    cbar_kws={'label':'AUCell z-score (clip ±3)',
              'orientation':'horizontal'},
    cbar_pos=(0.09, 0.965, 0.28, 0.014),            # wider horizontal bar TOP-LEFT, above dendrogram — fits full label
    xticklabels=True, yticklabels=True,
    linewidths=0, linecolor='none',
)

# horizontal y-labels, larger font
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=9)
# tilted x-labels
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10.5)
g.ax_heatmap.set_xlabel(''); g.ax_heatmap.set_ylabel('')

# colorbar cosmetics
g.ax_cbar.tick_params(labelsize=8)
g.ax_cbar.xaxis.set_label_position('top')
# left-anchor the colorbar label so long text extends RIGHT of the axis, not off the left edge
g.ax_cbar.set_xlabel('AUCell z-score  (state-wise, clip ±3)',
                     fontsize=9, color='#0f2747', labelpad=4,
                     loc='left')

# title as fig.text — sits ABOVE the top dendrogram + colorbar row, centered on the heatmap column
g.fig.text(0.58, 0.988,
           'pySCENIC · top-5 regulons per state × 12 tumor states  '
           f'({htm.shape[0]} regulons union) — clustered both axes',
           ha='center', va='top',
           fontsize=12, fontweight='bold', color='#0f2747')

# a tiny bit of top padding — leave the current subplot geometry from clustermap in place otherwise
g.fig.subplots_adjust(top=0.945)

for fmt in ['png','pdf']:
    g.fig.savefig(OUT/f'top5_regulons_per_state_heatmap.{fmt}',
                  dpi=180, bbox_inches='tight', facecolor='white')
plt.close(g.fig)
print(f'[saved] {OUT/"top5_regulons_per_state_heatmap.png"}')
print(f'[saved] {OUT/"top5_regulons_per_state_heatmap.pdf"}')

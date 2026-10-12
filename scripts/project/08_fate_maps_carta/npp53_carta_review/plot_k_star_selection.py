#!/usr/bin/env python
"""K* selection figure for §10.3 · Fig 10.3A.
Left panel: histogram of chosen k* across 25 samples (grouped by arm).
Right panel: k* vs n_states scatter, colored by arm; k* is bounded by n_states.
Shows the chord-elbow method result.
"""
import pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12')
OUT  = BASE / 'biology_review'
df = pd.read_csv(BASE / 'best_k_per_sample_v2_chord.csv')
print(df[['sample','arm','n_states','k_star_chord']].head())

ARM_ORDER  = ['Intact','Cas','Cas+Enz']
ARM_DISPLAY = {'Intact':'Intact','Cas':'Cas','Cas+Enz':'Cas+Enza'}
ARM_COLORS = {'Intact':'#1f5fa8','Cas':'#e0955f','Cas+Enz':'#c44e52'}

fig, ax = plt.subplots(figsize=(9, 6))  # single-panel: histogram only
k_range = sorted(df['k_star_chord'].unique())
bottom = {k: 0 for k in k_range}
for arm in ARM_ORDER:
    sub = df[df['arm']==arm]
    counts = sub['k_star_chord'].value_counts().reindex(k_range, fill_value=0)
    ax.bar(k_range, counts.values, bottom=[bottom[k] for k in k_range],
           color=ARM_COLORS[arm], edgecolor='white', linewidth=1.2,
           label=f'{ARM_DISPLAY[arm]} (n={len(sub)})', width=0.7)
    for k in k_range: bottom[k] += counts[k]

for k in k_range:
    total = bottom[k]
    if total > 0:
        ax.text(k, total + 0.15, f'n={total}', ha='center', fontsize=10, fontweight='bold')

ax.set_xticks(k_range)
ax.set_xlabel('k* (chord-elbow selection)', fontsize=11)
ax.set_ylabel('Number of samples', fontsize=11)
ax.set_title('(a) k* distribution across 25 samples\n(stacked by treatment arm)',
             fontsize=12, fontweight='bold', color='#0f2747')
# Panel (a) legend: kept inside since the histogram doesn't have data at the top-right
# where the legend sits — this is the ONLY legend in the figure (panel (b) has none;
# the outside-right legend on panel (b) is the master legend showing the same items
# once). To avoid two-legend duplication reported in QC, DROP the panel-(a) legend
# and rely on panel-(b)'s outside-right legend as the single figure-wide legend.
# ax.legend(...)   # <-- intentionally omitted
ax.legend(fontsize=10, loc='upper right', frameon=False, title='Treatment arm')
ax.grid(True, axis='y', alpha=0.3)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_ylim(0, max(bottom.values()) + 1.5)



fig.suptitle('CARTA k* selection (chord-elbow method on obj / n_leaves curves)',
             fontsize=13.5, fontweight='bold', color='#0f2747', y=1.02)
plt.subplots_adjust(top=0.88, bottom=0.10, left=0.06, right=0.86, wspace=0.30)
for fmt in ['png','pdf']:
    fig.savefig(OUT/f'npp53_carta_k_star_selection.{fmt}', dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'saved: {OUT}/npp53_carta_k_star_selection.png')

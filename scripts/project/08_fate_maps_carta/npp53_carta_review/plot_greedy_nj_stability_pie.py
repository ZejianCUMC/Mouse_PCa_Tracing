#!/usr/bin/env python
"""CARTA cross-solver stability pie — terminal-progenitor set, 2 categories only."""
import pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12/biology_review')
TSV  = BASE / 'cross_solver_stability.tsv'
OUT  = BASE

df = pd.read_csv(TSV, sep='\t')

def categorize(j):
    return 'Identical' if j >= 0.999 else 'Difference'

CAT_ORDER = ['Identical','Difference']
CAT_COLORS = {'Identical':'#2ecc71','Difference':'#e74c3c'}

df['cat_terminal'] = df['jaccard_terminal'].apply(categorize)

fig, ax = plt.subplots(1, 1, figsize=(7, 6))
import json as _j; _dd=_j.load(open(BASE/'discordant_at_adopted_k_gt20.json'))
counts = pd.Series({'Identical':len(_dd['concordant']),'Difference':len(_dd['discordant'])}).reindex(CAT_ORDER, fill_value=0)
total = counts.sum()
labels = [f'{c}\nn={counts[c]} ({counts[c]/total*100:.0f}%)' for c in CAT_ORDER]
wedges, texts = ax.pie(
    counts.values, labels=labels, startangle=90,
    colors=[CAT_COLORS[c] for c in CAT_ORDER],
    wedgeprops=dict(linewidth=1.5, edgecolor='white'),
    textprops=dict(fontsize=12, fontweight='bold', color='#111'),
)
ax.set_title('Full CARTA DAG at each sample\u2019s adopted k*  ·  >20-cell states',
             fontsize=13, fontweight='bold', color='#0f2747', pad=12)

fig.suptitle('CARTA cross-solver stability  ·  Greedy vs NJ trees (n = 25)',
             fontsize=14, fontweight='bold', color='#0f2747', y=1.0)
fig.text(0.5, -0.01,
         'Identical = Greedy and NJ produce the same DAG (>20-cell states) at the sample\u2019s adopted k*;  '
         'Difference = the two maps differ',
         ha='center', fontsize=10, color='#555', style='italic')

for fmt in ['png','pdf']:
    fig.savefig(OUT/f'npp53_carta_greedy_vs_nj_pie.{fmt}', dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'saved: {OUT}/npp53_carta_greedy_vs_nj_pie.png')
print(dict(df['cat_terminal'].value_counts()))

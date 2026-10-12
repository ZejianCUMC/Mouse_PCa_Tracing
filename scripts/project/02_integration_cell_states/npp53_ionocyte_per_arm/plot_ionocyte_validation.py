#!/usr/bin/env python
# =============================================================================
# Ionocyte validation figure for report §7.2
# =============================================================================
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from ete3 import Tree
from collections import Counter

M   = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
OUT = M/'1-mouse_based/Npp53_tracing/analysis/evolution_features_v2'

sc_fp = pd.read_csv(M/'1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity.csv')

# Panel a: per-sample Ionocyte counts colored by arm
ion_counts = sc_fp[sc_fp['putative_state']=='Ionocyte'].groupby(['sample','treatment']).size().reset_index(name='n_ion')
arm_order = ['Intact','Castration','Castration + Enzalutamide']
arm_short = {'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}
arm_color = {'Intact':'#4a7fbf','Castration':'#e0955f','Castration + Enzalutamide':'#c44e52'}

# Panel b: normalized leaf-depth ratio Ion/Other + PATH z for the 10 samples with meaningful Ion counts
samples_pos = list(ion_counts[ion_counts['n_ion']>=4].sort_values(['treatment','n_ion'], ascending=[True, False])['sample'])
depth_rows = []
for s in samples_pos:
    treef = M/f'cassiopeia_trees_alt/{s}_tree_maxcut.nwk'
    if not treef.exists(): continue
    tr = Tree(str(treef), format=1)
    max_depth = 0
    leaf_depth = {}
    for l in tr.get_leaves():
        d = 0; cur = l
        while cur.up is not None:
            d += 1; cur = cur.up
        leaf_depth[l.name] = d
        max_depth = max(max_depth, d)
    sub = sc_fp[sc_fp['sample']==s].drop_duplicates('cellBC').set_index('cellBC')
    sub_matched = sub.reindex([l for l in leaf_depth if l in sub.index])
    sub_matched['norm_depth'] = [leaf_depth[c]/max_depth for c in sub_matched.index]
    ion_mask = sub_matched['putative_state']=='Ionocyte'
    ion_d = sub_matched.loc[ion_mask,'norm_depth']
    oth_d = sub_matched.loc[~ion_mask,'norm_depth']
    depth_rows.append(dict(sample=s, treatment=sub_matched['treatment'].iloc[0],
                            n_ion=ion_mask.sum(),
                            ion_depth=ion_d.mean(), other_depth=oth_d.mean()))
depth_df = pd.DataFrame(depth_rows)

# PATH heritability for Ionocyte state
path_h = pd.read_csv(OUT/'PATH_out/PATH_heritability_long.csv')
ion_h = path_h[path_h['state']=='Ionocyte'].sort_values('z', ascending=False)

# Panel d: sister-leaf state distribution (JZ204 example)
s = 'JZ204'
treef = M/f'cassiopeia_trees_alt/{s}_tree_maxcut.nwk'
tr = Tree(str(treef), format=1)
sub = sc_fp[sc_fp['sample']==s].drop_duplicates('cellBC').set_index('cellBC')
state_map = sub['putative_state'].to_dict()
sister_counts = Counter()
for l in tr.get_leaves():
    if state_map.get(l.name) != 'Ionocyte': continue
    parent = l.up
    if parent is None: continue
    for sib in [c for c in parent.get_leaves() if c.name != l.name]:
        st = state_map.get(sib.name)
        if st: sister_counts[st] += 1

# --- FIGURE ---
fig, axes = plt.subplots(2, 2, figsize=(15, 9))

# Panel a — Ionocyte counts per sample
ax = axes[0,0]
ion_sorted = ion_counts.copy()
ion_sorted['arm_order'] = ion_sorted['treatment'].apply(arm_order.index)
ion_sorted = ion_sorted.sort_values(['arm_order','n_ion'], ascending=[True, False])
for i, (_, row) in enumerate(ion_sorted.iterrows()):
    ax.bar(i, row['n_ion'], color=arm_color[row['treatment']], edgecolor='k', linewidth=0.3)
ax.set_xticks(range(len(ion_sorted)))
ax.set_xticklabels(ion_sorted['sample'], rotation=45, fontsize=8)
ax.set_ylabel('Ionocyte cell count')
ax.set_title(f'a. Ionocyte counts per sample (n={ion_sorted["n_ion"].sum()}; ~76% in treated arms)', fontsize=11, loc='left')
# legend
for arm in arm_order:
    ax.bar(0, 0, color=arm_color[arm], label=arm_short[arm])
ax.legend(loc='upper right', fontsize=9)
ax.grid(alpha=0.3, axis='y')

# Panel b — normalized leaf depth: Ion vs Other per sample
ax = axes[0,1]
x = np.arange(len(depth_df))
w = 0.35
ax.bar(x-w/2, depth_df['ion_depth'], w, label='Ionocyte',    color='#c44e52', edgecolor='k', linewidth=0.3)
ax.bar(x+w/2, depth_df['other_depth'], w, label='All other', color='#888', edgecolor='k', linewidth=0.3)
ax.set_xticks(x); ax.set_xticklabels(depth_df['sample'], rotation=45, fontsize=8)
ax.set_ylabel('Normalized leaf depth (0=root, 1=deepest)')
ax.set_title(f'b. Ionocytes sit ~5% shallower than average across samples (weak ancestor-like)', fontsize=11, loc='left')
ax.legend(fontsize=9); ax.grid(alpha=0.3, axis='y')
# annotate mean ratio
mean_ratio = (depth_df['ion_depth'] / depth_df['other_depth']).mean()
ax.text(0.98, 0.95, f'mean(Ion/Other) = {mean_ratio:.3f}', transform=ax.transAxes,
        ha='right', va='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

# Panel c — PATH heritability z per sample
ax = axes[1,0]
ion_h_plot = ion_h.copy()
# ion_h already has 'treatment' from PATH_heritability_long.csv
ion_h_plot['arm_order'] = ion_h_plot['treatment'].apply(arm_order.index)
ion_h_plot = ion_h_plot.sort_values(['arm_order','z'], ascending=[True, False])
colors = [arm_color[a] for a in ion_h_plot['treatment']]
ax.bar(range(len(ion_h_plot)), ion_h_plot['z'], color=colors, edgecolor='k', linewidth=0.3)
ax.axhline(2, color='red', ls='--', lw=1, label='z > 2 = detectable heritability (PATH)')
ax.axhline(0, color='k', lw=0.5)
ax.set_xticks(range(len(ion_h_plot))); ax.set_xticklabels(ion_h_plot['sample'], rotation=45, fontsize=8)
ax.set_ylabel('PATH z-score (Ionocyte heritability)')
ax.set_title(f'c. Ionocyte is lineage-heritable in Cas/CasEnz + JZ214 (z>2); scattered elsewhere', fontsize=11, loc='left')
ax.legend(loc='upper right', fontsize=9); ax.grid(alpha=0.3, axis='y')

# Panel d — sister state distribution (JZ204)
ax = axes[1,1]
sister_df = pd.DataFrame(sister_counts.most_common(10), columns=['state','count'])
colors = ['#c44e52' if s=='Ionocyte' else '#4a7fbf' if s=='Basal-like' else '#888' for s in sister_df['state']]
ax.barh(range(len(sister_df)), sister_df['count'], color=colors, edgecolor='k', linewidth=0.3)
ax.set_yticks(range(len(sister_df))); ax.set_yticklabels(sister_df['state'], fontsize=9)
ax.invert_yaxis()
ax.set_xlabel('# sister-leaf pairs')
ax.set_title(f'd. Ionocyte closest sister states (JZ204, 172 pairs): Basal-like dominant', fontsize=11, loc='left')
ax.grid(alpha=0.3, axis='x')

fig.suptitle('§7.2 · Ionocyte trajectory validation (n=297 Ionocytes across 25 mice)', fontsize=14, y=1.00)
fig.tight_layout()
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'Fig_7_2_Ionocyte.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'saved: {OUT}/Fig_7_2_Ionocyte.{{pdf,png}}')

# also save the tabular data for report inclusion
ion_counts.to_csv(OUT/'ionocyte_counts_per_sample.csv', index=False)
depth_df.to_csv(OUT/'ionocyte_depth_per_sample.csv', index=False)
ion_h.to_csv(OUT/'ionocyte_PATH_heritability.csv', index=False)
sister_df.to_csv(OUT/'ionocyte_JZ204_sister_states.csv', index=False)
print('saved: 4 companion CSV files')

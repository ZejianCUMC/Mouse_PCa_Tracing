#!/usr/bin/env python
"""Test whether the observed clonal dominance is greedy-specific (algorithm)
or shows up across all tree-building methods (biology).

Compare dominant subclade fraction at depth 1/2/3/5 across 5 solvers
(greedy, nj, spectral, maxcut, percolation) on the SAME per-sample character
matrices. If greedy is uniquely lopsided → algorithmic artifact. If all
solvers agree → biological signal.

Also runs a random-permutation null: shuffle leaves-to-tree-positions on the
greedy topology; the resulting dominant fraction distribution is a stringent
null that tests just the tree shape.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ete3 import Tree

M = '/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
OUT = Path(f'{M}/qc')

TREE_DIRS = {
    'greedy':      Path(f'{M}/analysis/hpc_tree_rebuild/v8_trees/greedy'),
    'nj':          Path(f'{M}/analysis/hpc_tree_rebuild/v8_trees/nj'),
    'spectral':    Path(f'{M}/analysis/hpc_tree_rebuild/v8_trees/spectral'),
    'maxcut':      Path(f'{M}/lineage_tracing/cassiopeia_trees_v8_rebuilt/maxcut'),
    'percolation': Path(f'{M}/lineage_tracing/cassiopeia_trees_v8_rebuilt/percolation'),
}
SOLVER_COLORS = {'greedy':'#1f77b4','nj':'#e0955f','spectral':'#2ca02c',
                 'maxcut':'#9467bd','percolation':'#8c564b'}

scfp = pd.read_csv(f'{M}/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv')
arm_of = scfp.groupby('sample')['treatment'].first().replace(
    {'Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}).to_dict()
ARM_COLORS = {'Intact':'#4a7fbf', 'Cas':'#e0955f', 'Cas+Enz':'#c44e52'}
ARM_ORDER  = ['Intact','Cas','Cas+Enz']

def frontier_at_depth(root, d):
    if d == 0: return [root]
    fr = [root]
    for _ in range(d):
        nf = []
        for n in fr:
            if n.is_leaf(): nf.append(n)
            else: nf.extend(n.children)
        fr = nf
    return fr

def dom_frac_at_depth(tree, d):
    sizes = [len(n.get_leaves()) for n in frontier_at_depth(tree, d)]
    return max(sizes) / sum(sizes)

# --- read all trees, compute dom_frac at multiple depths ---
DEPTHS = [1, 2, 3, 5]
rows = []
for solver, tdir in TREE_DIRS.items():
    print(f'Reading {solver}...')
    for f in sorted(tdir.glob('*_tree_*.nwk')):
        s = f.stem.split('_tree_')[0]
        # sanity: solver in filename
        try:
            t = Tree(str(f), format=1)
        except Exception:
            try: t = Tree(str(f), format=0)
            except Exception: continue
        row = {'sample': s, 'solver': solver, 'arm': arm_of.get(s,'?'), 'n_leaves': len(t.get_leaves())}
        for d in DEPTHS:
            row[f'dom_frac_d{d}'] = dom_frac_at_depth(t, d)
        rows.append(row)
df = pd.DataFrame(rows)

# --- keep only samples where we have all 5 solvers ---
sample_counts = df.groupby('sample')['solver'].nunique()
full = set(sample_counts[sample_counts == 5].index)
print(f'Samples with all 5 solvers: {len(full)}/{df["sample"].nunique()}')
df_full = df[df['sample'].isin(full)].copy()

df_full.to_csv(OUT/'dominance_cross_solver.csv', index=False)

# --- Random-permutation null on greedy topology ---
# For each greedy tree, shuffle the leaf order (relabel leaves), recompute dom_frac.
# The tree TOPOLOGY is fixed, only which cell sits where changes.
rng = np.random.default_rng(0)
null_rows = []
for f in sorted(TREE_DIRS['greedy'].glob('*_tree_greedy.nwk')):
    s = f.stem.replace('_tree_greedy','')
    try:
        t = Tree(str(f), format=1)
    except Exception:
        t = Tree(str(f), format=0)
    if s not in full: continue
    # baseline (real assignment)
    for d in DEPTHS:
        obs = dom_frac_at_depth(t, d)
        # null: permute leaf identities (topology unchanged) — dom_frac invariant under this!
        # Since dom_frac only depends on subtree sizes (not identities), this null gives THE SAME value.
        # So instead we need a DIFFERENT null — shuffle the topology.
    # Real null: build a random binary tree with same n leaves, compute dom_frac_d
    n = len(t.get_leaves())
    for _ in range(50):
        # random binary tree by iterative pair-joining
        active = [Tree(name=str(i)) for i in range(n)]
        while len(active) > 1:
            rng.shuffle(active)
            a, b = active.pop(), active.pop()
            parent = Tree(); parent.add_child(a); parent.add_child(b)
            active.append(parent)
        rt = active[0]
        for d in DEPTHS:
            null_rows.append({'sample': s, 'depth': d, 'dom_frac': dom_frac_at_depth(rt, d)})
null_df = pd.DataFrame(null_rows)

# --- Print summary ---
print('\n=== Median dom_frac at each depth per solver ===')
piv = df_full.groupby('solver')[[f'dom_frac_d{d}' for d in DEPTHS]].median() * 100
print(piv.round(1))
print('\n=== Random-binary-tree null (median dom_frac) ===')
print((null_df.groupby('depth')['dom_frac'].median()*100).round(1))

# --- FIGURE ---
fig = plt.figure(figsize=(15, 9))
gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.28, height_ratios=[1, 1])

# Panel A: box per solver per depth
ax = fig.add_subplot(gs[0, 0])
SOLVER_ORDER = ['greedy','nj','spectral','maxcut','percolation']
x = 0; xticks = []; xlabels = []
for d in DEPTHS:
    for solver in SOLVER_ORDER:
        vals = df_full.loc[df_full['solver']==solver, f'dom_frac_d{d}'].values * 100
        bp = ax.boxplot([vals], positions=[x], widths=0.65, patch_artist=True,
                        boxprops=dict(facecolor=SOLVER_COLORS[solver], alpha=0.55, color='#333'),
                        medianprops=dict(color='#111', linewidth=2),
                        whiskerprops=dict(color='#666'), capprops=dict(color='#666'), showfliers=False)
        j = np.random.default_rng(x).uniform(-0.15,0.15,size=len(vals))
        ax.scatter(np.full(len(vals),x)+j, vals, s=14, c=SOLVER_COLORS[solver], edgecolors='#333', lw=0.3, alpha=0.85, zorder=3)
        xticks.append(x); xlabels.append(solver[:3])
        x += 1
    # depth separator
    ax.axvline(x-0.5, color='#ccc', ls=':', lw=0.7); x += 0.5
    # depth label above
    ax.text(x - 3.5, 106, f'd={d}', ha='center', va='bottom', fontsize=11, fontweight='bold', color='#0f2747')
# random-binary-tree null median as a red dashed line per depth
xpos_offset = 0
for i, d in enumerate(DEPTHS):
    null_med = null_df[null_df['depth']==d]['dom_frac'].median() * 100
    ax.hlines(null_med, xmin=i*(len(SOLVER_ORDER)+0.5)-0.5, xmax=i*(len(SOLVER_ORDER)+0.5)+len(SOLVER_ORDER)-0.5,
              colors='#c44e52', linestyles='--', linewidth=1.4, zorder=2,
              label='random binary tree null' if i == 0 else None)
ax.set_xticks(xticks); ax.set_xticklabels(xlabels, rotation=0, fontsize=8)
ax.set_ylabel('Dominant subclade %', fontsize=11)
ax.set_title('A · Cross-solver dominance @ depths 1/2/3/5 — random-binary-null in dashed red',
             fontsize=11, fontweight='bold', color='#0f2747')
ax.set_ylim(0, 110); ax.axhline(50, color='#888', ls=':', lw=0.5)
ax.legend(loc='lower left', fontsize=8)
ax.grid(axis='y', alpha=0.25)

# Panel B: per-sample greedy vs NJ scatter at depth 3
ax = fig.add_subplot(gs[0, 1])
for arm in ARM_ORDER:
    sub = df_full[df_full['arm']==arm]
    g = sub[sub['solver']=='greedy'].set_index('sample')['dom_frac_d3']*100
    nj = sub[sub['solver']=='nj'].set_index('sample')['dom_frac_d3']*100
    common = g.index.intersection(nj.index)
    ax.scatter(g[common], nj[common], c=ARM_COLORS[arm], s=60, edgecolor='#111', lw=0.5,
               label=arm, alpha=0.85, zorder=3)
ax.plot([0,100],[0,100], color='#888', ls='--', lw=0.7, zorder=1)
ax.set_xlabel('Greedy — dom_frac @ d=3 (%)', fontsize=11)
ax.set_ylabel('NJ — dom_frac @ d=3 (%)', fontsize=11)
ax.set_title('B · Per-sample greedy vs NJ (does the algorithm choice matter?)',
             fontsize=11, fontweight='bold', color='#0f2747')
ax.set_xlim(30, 105); ax.set_ylim(30, 105)
ax.legend(loc='lower right', fontsize=9); ax.grid(alpha=0.3)

# Panel C: per-sample delta greedy - NJ at each depth
ax = fig.add_subplot(gs[1, 0])
delta_rows = []
for s in full:
    for d in DEPTHS:
        g = df_full[(df_full['sample']==s)&(df_full['solver']=='greedy')][f'dom_frac_d{d}'].values
        nj = df_full[(df_full['sample']==s)&(df_full['solver']=='nj')][f'dom_frac_d{d}'].values
        if len(g)==1 and len(nj)==1:
            delta_rows.append({'sample':s, 'arm':arm_of[s], 'depth':d, 'delta_pct':(g[0]-nj[0])*100})
dlt = pd.DataFrame(delta_rows)
for i, d in enumerate(DEPTHS):
    vals = dlt.loc[dlt['depth']==d, 'delta_pct'].values
    j = np.random.default_rng(i).uniform(-0.15, 0.15, size=len(vals))
    ax.scatter(np.full(len(vals),i)+j, vals, s=32, c=['#4a7fbf' if arm_of[s]=='Intact' else '#e0955f' if arm_of[s]=='Cas' else '#c44e52' for s in dlt.loc[dlt['depth']==d,'sample']],
               edgecolors='#111', lw=0.3, alpha=0.85, zorder=3)
    ax.hlines(np.median(vals), i-0.3, i+0.3, colors='#111', lw=2, zorder=4)
ax.axhline(0, color='#c44e52', ls='--', lw=1)
ax.set_xticks(range(len(DEPTHS))); ax.set_xticklabels([f'd={d}' for d in DEPTHS], fontsize=11)
ax.set_ylabel('Greedy − NJ  (percentage-point diff)', fontsize=11)
ax.set_title('C · Per-sample greedy-vs-NJ dominance gap (0 = algorithm-agnostic)',
             fontsize=11, fontweight='bold', color='#0f2747')
ax.grid(alpha=0.3)

# Panel D: same for greedy vs spectral (another non-greedy solver)
ax = fig.add_subplot(gs[1, 1])
delta_rows2 = []
for s in full:
    for d in DEPTHS:
        g = df_full[(df_full['sample']==s)&(df_full['solver']=='greedy')][f'dom_frac_d{d}'].values
        sp = df_full[(df_full['sample']==s)&(df_full['solver']=='spectral')][f'dom_frac_d{d}'].values
        if len(g)==1 and len(sp)==1:
            delta_rows2.append({'sample':s, 'arm':arm_of[s], 'depth':d, 'delta_pct':(g[0]-sp[0])*100})
dlt2 = pd.DataFrame(delta_rows2)
for i, d in enumerate(DEPTHS):
    vals = dlt2.loc[dlt2['depth']==d, 'delta_pct'].values
    j = np.random.default_rng(i).uniform(-0.15, 0.15, size=len(vals))
    ax.scatter(np.full(len(vals),i)+j, vals, s=32, c=['#4a7fbf' if arm_of[s]=='Intact' else '#e0955f' if arm_of[s]=='Cas' else '#c44e52' for s in dlt2.loc[dlt2['depth']==d,'sample']],
               edgecolors='#111', lw=0.3, alpha=0.85, zorder=3)
    ax.hlines(np.median(vals), i-0.3, i+0.3, colors='#111', lw=2, zorder=4)
ax.axhline(0, color='#c44e52', ls='--', lw=1)
ax.set_xticks(range(len(DEPTHS))); ax.set_xticklabels([f'd={d}' for d in DEPTHS], fontsize=11)
ax.set_ylabel('Greedy − Spectral  (pp diff)', fontsize=11)
ax.set_title('D · Same test vs Spectral solver',
             fontsize=11, fontweight='bold', color='#0f2747')
ax.grid(alpha=0.3)

fig.suptitle('Is dominance a greedy-solver artifact? — cross-solver comparison',
             fontsize=13, fontweight='bold', color='#0f2747', y=0.995)
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'dominance_cross_solver.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'\nsaved: {OUT}/dominance_cross_solver.png')

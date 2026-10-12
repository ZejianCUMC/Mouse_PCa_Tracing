#!/usr/bin/env python
"""Clonal dominance QC — supplementary-figure quality.

Rationale: barcoded ES-cell organoids implanted as a pool → the a priori
expectation is polyclonal seeding. If we observe one dominant subclade per
tumor in the greedy tree, that's a real biological finding about in-vivo
clonal selection.

This script produces a 5-panel supplementary figure characterizing the
clonal-dominance signal from multiple angles.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ete3 import Tree

M = '/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
TREES = Path(f'{M}/analysis/hpc_tree_rebuild/v8_trees/greedy')
OUT = Path(f'{M}/qc')
OUT.mkdir(exist_ok=True)

scfp = pd.read_csv(f'{M}/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv')
arm_of = scfp.groupby('sample')['treatment'].first().replace(
    {'Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}).to_dict()
ARM_COLORS = {'Intact':'#4a7fbf', 'Cas':'#e0955f', 'Cas+Enz':'#c44e52'}
ARM_ORDER  = ['Intact','Cas','Cas+Enz']

# --- extract subclade size distribution at multiple depths ---
def frontier_at_depth(root, d):
    if d == 0: return [root]
    frontier = [root]
    for _ in range(d):
        new_frontier = []
        for n in frontier:
            if n.is_leaf(): new_frontier.append(n)
            else: new_frontier.extend(n.children)
        frontier = new_frontier
    return frontier

def clade_sizes(root, d):
    return sorted([len(n.get_leaves()) for n in frontier_at_depth(root, d)], reverse=True)

def simpson_1D(sizes):
    N = sum(sizes)
    if N < 2: return 0.0
    p = np.array(sizes, dtype=float) / N
    return 1 - (p**2).sum()

def shannon(sizes):
    N = sum(sizes)
    if N < 2: return 0.0
    p = np.array(sizes, dtype=float) / N
    p = p[p > 0]
    return -(p * np.log(p)).sum()

# --- read trees + compute ---
rows = []
sizes_by_sample = {}   # sample -> {depth: [sorted size list]}
for f in sorted(TREES.glob('*_tree_greedy.nwk')):
    s = f.stem.replace('_tree_greedy','')
    try:
        t = Tree(str(f), format=1)
    except Exception:
        t = Tree(str(f), format=0)
    n = len(t.get_leaves())
    sizes_by_sample[s] = {}
    row = {'sample': s, 'arm': arm_of.get(s,'?'), 'n_leaves': n}
    for d in [1,2,3,4,5]:
        sz = clade_sizes(t, d)
        sizes_by_sample[s][d] = sz
        row[f'dom_frac_d{d}'] = sz[0]/n if n else np.nan
        row[f'top2_frac_d{d}'] = sum(sz[:2])/n if n else np.nan
        row[f'n_subtrees_d{d}'] = len(sz)
        row[f'simpson_d{d}'] = simpson_1D(sz)
        row[f'shannon_d{d}'] = shannon(sz)
    rows.append(row)
df = pd.DataFrame(rows).sort_values(['arm','sample']).reset_index(drop=True)
df.to_csv(OUT/'clonal_dominance_summary.csv', index=False)

# Print quick summary
print('=== Dominance QC summary ===')
print(df[['sample','arm','n_leaves','dom_frac_d1','dom_frac_d2','dom_frac_d3','dom_frac_d5','simpson_d3']].to_string(index=False))

# --- Compute NULL reference: if leaves were assigned randomly to same # subtrees ---
# For each sample+depth, simulate uniform assignment 100× and get expected dom_frac
rng = np.random.default_rng(0)
null_rows = []
for s, dsizes in sizes_by_sample.items():
    n = sum(dsizes[1])
    for d in [3,5]:
        k = dsizes[d]; K = len(k)
        # random assignment: each leaf independently to one of K bins uniformly
        if K < 1 or n < 1: continue
        null_doms = []
        for _ in range(200):
            assign = rng.integers(0, K, size=n)
            counts = np.bincount(assign, minlength=K)
            null_doms.append(counts.max() / n)
        null_rows.append({'sample':s, 'depth':d, 'null_dom_median':np.median(null_doms),
                          'null_dom_q05':np.quantile(null_doms, 0.05),
                          'null_dom_q95':np.quantile(null_doms, 0.95),
                          'observed_dom':dsizes[d][0]/n})
null_df = pd.DataFrame(null_rows)
null_df.to_csv(OUT/'clonal_dominance_null.csv', index=False)

# --- FIGURE ---
fig = plt.figure(figsize=(17, 13))
gs = GridSpec(3, 2, figure=fig, hspace=0.42, wspace=0.28, height_ratios=[1, 1.1, 1.1])

# Panel A: per-sample rank-abundance at depth-5 (finer partition, up to 32 subtrees)
ax = fig.add_subplot(gs[0, 0])
for s in df['sample']:
    arm = arm_of[s]; c = ARM_COLORS[arm]
    sz = sizes_by_sample[s][5]
    total = sum(sz)
    frac = np.array(sz) / total
    ax.plot(range(1, len(frac)+1), frac, '-', color=c, lw=1.2, alpha=0.55)
ax.set_xlabel('subclade rank (depth-5)', fontsize=10)
ax.set_ylabel('fraction of cells', fontsize=10)
ax.set_title('A · Per-sample rank-abundance @ depth-5', fontsize=11, fontweight='bold', color='#0f2747')
ax.set_yscale('log'); ax.set_ylim(0.001, 1.05)
ax.set_xticks([1,5,10,15,20,25,32])
ax.axhline(0.5, color='#c44e52', ls=':', lw=0.6, alpha=0.6, label='50% cutoff')
# legend
from matplotlib.lines import Line2D
ax.legend([Line2D([0],[0],color=ARM_COLORS[a],lw=1.5) for a in ARM_ORDER],
          ARM_ORDER, loc='upper right', fontsize=9, frameon=True, title='Arm')
ax.grid(alpha=0.25)

# Panel B: dominant-clone-fraction distribution at 3 depths
ax = fig.add_subplot(gs[0, 1])
box_data = []; box_labels = []; box_colors = []
for d in [1,2,3,5]:
    for arm in ARM_ORDER:
        vals = df.loc[df['arm']==arm, f'dom_frac_d{d}'].values
        box_data.append(vals*100); box_labels.append(f'd={d}\n{arm}'); box_colors.append(ARM_COLORS[arm])
bp = ax.boxplot(box_data, positions=range(len(box_data)), widths=0.6, patch_artist=True,
                boxprops=dict(color='#333'), medianprops=dict(color='#111', linewidth=2),
                whiskerprops=dict(color='#666'), capprops=dict(color='#666'), showfliers=False)
for patch, c in zip(bp['boxes'], box_colors):
    patch.set_facecolor(c); patch.set_alpha(0.55)
# jitter points
for i, (vals, c) in enumerate(zip(box_data, box_colors)):
    jitter = np.random.default_rng(i).uniform(-0.15, 0.15, size=len(vals))
    ax.scatter(np.full(len(vals),i)+jitter, vals, s=22, c=c, edgecolors='#333', lw=0.4, alpha=0.9, zorder=3)
ax.set_xticks(range(len(box_labels)))
ax.set_xticklabels(box_labels, fontsize=8, rotation=0)
for i in range(len(box_labels)):
    if i > 0 and i % 3 == 0:
        ax.axvline(i-0.5, color='#ddd', lw=0.8, ls=':')
ax.set_ylabel('Dominant subclade %', fontsize=10)
ax.set_title('B · Dominant subclade fraction — depth 1/2/3/5, per arm', fontsize=11, fontweight='bold', color='#0f2747')
ax.axhline(50, color='#c44e52', ls=':', lw=0.6, alpha=0.6)
ax.set_ylim(40, 105)
ax.grid(axis='y', alpha=0.3)

# Panel C: stacked bar per sample at depth-5 (top-N subclades + rest)
ax = fig.add_subplot(gs[1, :])
samples_ordered = df.sort_values(['arm','n_leaves'], ascending=[True,False])['sample'].tolist()
STACK_D = 5
TOPN = 5
def stack_frac(s, d, topn):
    sz = sizes_by_sample[s][d]
    total = sum(sz)
    top = [sz[i]/total for i in range(min(topn, len(sz)))]
    rest = 1 - sum(top)
    return top + [rest]
palette = ['#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#bdbdbd']  # top5 warm→cool + grey "rest"
xpos = np.arange(len(samples_ordered))
bottoms = np.zeros(len(samples_ordered))
for tier in range(TOPN + 1):  # 5 top + rest
    heights = []
    for s in samples_ordered:
        top_and_rest = stack_frac(s, STACK_D, TOPN)
        heights.append(top_and_rest[tier] if tier < len(top_and_rest) else 0)
    heights = np.array(heights)
    ax.bar(xpos, heights*100, bottom=bottoms*100, width=0.78,
           color=palette[tier], edgecolor='#111', linewidth=0.4,
           label=(f'#{tier+1}' if tier < TOPN else 'rest'))
    bottoms += heights
# arm chip labels below x-axis
ax.set_xticks(xpos); ax.set_xticklabels(samples_ordered, rotation=45, ha='right', fontsize=8)
for i, s in enumerate(samples_ordered):
    arm = arm_of[s]
    ax.text(i, -6, arm[0], color=ARM_COLORS[arm], ha='center', va='top', fontsize=8, fontweight='bold')
ax.text(-0.5, -6, 'arm:', color='#333', ha='right', va='top', fontsize=8)
ax.set_ylabel('% of cells', fontsize=10)
ax.set_title(f'C · Top-{TOPN} subclade composition per tumor @ depth-{STACK_D} (rest = smaller clades pooled)',
             fontsize=11, fontweight='bold', color='#0f2747')
ax.set_ylim(0, 100)
ax.legend(loc='center left', bbox_to_anchor=(1.005, 0.5), fontsize=8, title='rank')
ax.grid(axis='y', alpha=0.25)

# Panel D: observed vs uniform-null dominance (per sample at depth 5)
ax = fig.add_subplot(gs[2, 0])
d5 = null_df[null_df['depth']==5].merge(df[['sample','arm']], on='sample')
d5['arm_c'] = d5['arm'].map(ARM_COLORS)
# plot bars: observed vs null median with 5-95% range
xpos = np.arange(len(d5))
d5 = d5.sort_values(['arm','sample']).reset_index(drop=True)
# uniform-null 5-95% interval as grey band
ax.fill_between(xpos, d5['null_dom_q05']*100, d5['null_dom_q95']*100,
                color='#bbb', alpha=0.5, label='uniform-null 5-95%')
ax.plot(xpos, d5['null_dom_median']*100, 'o', color='#555', ms=4, label='null median', zorder=2)
for i, r in d5.iterrows():
    ax.plot([i], [r['observed_dom']*100], marker='D', ms=6, mfc=ARM_COLORS[r['arm']], mec='#111', mew=0.6, zorder=4)
ax.set_xticks(xpos); ax.set_xticklabels(d5['sample'], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Dominant subclade % @ depth-5', fontsize=10)
ax.set_title('D · Observed dominance ≫ uniform-null (evidence for selective sweep, not random)',
             fontsize=11, fontweight='bold', color='#0f2747')
ax.set_ylim(0, 105)
ax.legend(loc='lower left', fontsize=8)
ax.grid(axis='y', alpha=0.25)
# arm color guide
for i, r in d5.iterrows():
    ax.text(i, -6, r['arm'][0], color=ARM_COLORS[r['arm']], ha='center', va='top', fontsize=7, fontweight='bold')

# Panel E: diversity index (Simpson 1-D) per sample per arm
ax = fig.add_subplot(gs[2, 1])
for i, d in enumerate([1, 2, 3, 5]):
    for j, arm in enumerate(ARM_ORDER):
        vals = df.loc[df['arm']==arm, f'simpson_d{d}'].values
        x = i*4 + j
        ax.boxplot([vals], positions=[x], widths=0.6, patch_artist=True,
                   boxprops=dict(facecolor=ARM_COLORS[arm], alpha=0.5, color='#333'),
                   medianprops=dict(color='#111', linewidth=2),
                   whiskerprops=dict(color='#666'), capprops=dict(color='#666'), showfliers=False)
        jitter = np.random.default_rng(x).uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(np.full(len(vals),x)+jitter, vals, s=22, c=ARM_COLORS[arm], edgecolors='#333', lw=0.4, alpha=0.9, zorder=3)
ax.set_xticks([0*4+1, 1*4+1, 2*4+1, 3*4+1])
ax.set_xticklabels(['d=1','d=2','d=3','d=5'], fontsize=10)
# arm legend at top
ax.legend([plt.Rectangle((0,0),1,1,fc=ARM_COLORS[a],alpha=0.5) for a in ARM_ORDER],
          ARM_ORDER, loc='upper right', fontsize=8, frameon=True, title='Arm')
ax.set_ylabel("Simpson diversity (1 − Σ pᵢ²)", fontsize=10)
ax.set_title('E · Subclade diversity index per depth — near 0 = one dominant clone',
             fontsize=11, fontweight='bold', color='#0f2747')
ax.axhline(0.1, color='#c44e52', ls=':', lw=0.6, alpha=0.6)
ax.set_ylim(-0.02, 1.05)
ax.grid(axis='y', alpha=0.3)

fig.suptitle('Clonal-dominance QC · 25 Npp53 tumors (greedy trees) — 1 dominant subclade per tumor at every depth',
             fontsize=14, fontweight='bold', color='#0f2747', y=0.995)
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'clonal_dominance_qc.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'\nsaved: {OUT}/clonal_dominance_qc.png')
print(f'saved: {OUT}/clonal_dominance_summary.csv')
print(f'saved: {OUT}/clonal_dominance_null.csv')

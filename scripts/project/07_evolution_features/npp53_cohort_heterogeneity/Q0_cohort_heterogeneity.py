#!/usr/bin/env python
"""§10 · Q0 · Cohort-heterogeneity sanity check (BEFORE Q2 doublecheck).

Question
--------
Do tumors within the SAME treatment arm share more phylogenetic evolutionary structure
with each other than they do with tumors from OTHER arms?

Interpretation
--------------
- within-arm > between-arm  → arm has a shared program → Q1/Q2 arm-level tests are LEGITIMATE
- within-arm ≈ between-arm → tumors are individually idiosyncratic → arm-level tests
  are effectively testing noise, and any Q1/Q2 arm-level p-values should be re-interpreted
  as descriptive-only

Method
------
Per-tumor phylogenetic-pattern representation = flattened Fitch 12×12 parent→child
edge-transition count matrix (n_gains + n_losses; row=src state, col=dst state).
25 tumors → 25 × 144-dim vectors. Compute pairwise distance (Pearson correlation
distance = 1 − Pearson r) and Bray-Curtis (compositional distance).

Test
----
PERMANOVA-style: F-statistic = variance BETWEEN arms / variance WITHIN arms on the
25×25 distance matrix. Permute arm labels 9999 times to build null → 2-sided p.

Bonus: MDS embedding of the 25 tumors in 2D + colored by arm — visual sanity check.

Outputs (analysis/evolution_features_v2/lineage_6Q/):
  Q0_pairwise_phylo_distance.csv              25×25
  Q0_arm_dispersion_stats.csv                 mean/CV within vs between per metric
  Q0_permanova_pvalues.csv                    F, p, effect size per metric
  Q0_within_vs_between.{png,pdf}              violin: within-arm vs between-arm distance
  Q0_pairwise_heatmap.{png,pdf}               25×25 heatmap ordered by arm
  Q0_MDS_embedding.{png,pdf}                  2D MDS of tumors colored by arm
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.spatial.distance import squareform, pdist
from scipy.stats import pearsonr
from sklearn.manifold import MDS

M   = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
Q1  = M/'analysis/evolution_features_v2/lineage_6Q/per_sample_transition_counts_Q1.csv'
SCFP = M/'analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT = M/'analysis/evolution_features_v2/lineage_6Q'

STATES = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
          'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
          'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
DISP = {s: f'c{i+1}_{s}' for i,s in enumerate(STATES)}
NS = len(STATES)
ARM_ORDER  = ['Intact','Castration','Castration + Enzalutamide']
ARM_SHORT  = {'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}
ARM_COLORS = {'Intact':'#4a7fbf','Cas':'#e0955f','Cas+Enz':'#c44e52'}

# ── build per-sample 12x12 Fitch transition matrices (from Q1 long-format CSV) ──
print('[build] per-sample 12×12 Fitch transition matrices ...')
q1 = pd.read_csv(Q1)   # sample, treatment, src_state, dst_state, n_transitions
samples = sorted(q1['sample'].unique())
scfp = pd.read_csv(SCFP)
sample_arm = scfp.groupby('sample')['treatment'].first().to_dict()

# state->index
DISP_REV = {v:k for k,v in DISP.items()}
def state_to_idx(disp_name):
    return STATES.index(DISP_REV[disp_name])

per_sample_mat = {}
for samp in samples:
    m = np.zeros((NS, NS), dtype=int)
    for _, r in q1[q1['sample']==samp].iterrows():
        i = state_to_idx(r['src_state'])
        j = state_to_idx(r['dst_state'])
        m[i,j] = int(r['n_transitions'])
    per_sample_mat[samp] = m

# flatten each 12×12 → 144-dim
X = np.array([per_sample_mat[s].flatten() for s in samples])   # 25 × 144
arms_short = np.array([ARM_SHORT[sample_arm[s]] for s in samples])
print(f'  built X: {X.shape}   arms: {dict(zip(*np.unique(arms_short, return_counts=True)))}')

# ── compute pairwise distances (2 metrics) ────────────────────────────────
def corr_distance(X):
    """1 − Pearson r between rows; ignores samples with zero variance."""
    N = X.shape[0]
    D = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if i == j: D[i,j] = 0; continue
            xi, xj = X[i], X[j]
            if xi.std() == 0 or xj.std() == 0: D[i,j] = np.nan
            else: D[i,j] = 1 - pearsonr(xi, xj)[0]
    return D

# 1. Pearson correlation distance on raw transition counts (0-inflated)
D_corr = corr_distance(X)
# 2. Bray-Curtis on per-sample-normalized transition matrices (compositional)
X_norm = X / (X.sum(axis=1, keepdims=True) + 1e-9)   # sum-to-1 per sample
D_bray = squareform(pdist(X_norm, metric='braycurtis'))

# ── PERMANOVA (arm labels vs distance) ────────────────────────────────────
def permanova_F(D, groups):
    """Anderson 2001 pseudo-F statistic. Returns (F, mean_within, mean_between)."""
    N = D.shape[0]
    # remove rows/cols with any NaN distance
    mask = ~np.isnan(D).any(axis=1)
    D = D[np.ix_(mask, mask)]; groups = groups[mask]
    N = D.shape[0]
    # within-group sum of squared distances
    unique_g = np.unique(groups)
    SS_within = 0
    for g in unique_g:
        idx = np.where(groups == g)[0]
        if len(idx) < 2: continue
        SS_within += np.sum(D[np.ix_(idx, idx)]**2) / (2 * len(idx))
    # total sum of squared distances
    SS_total = np.sum(D**2) / (2 * N)
    SS_between = SS_total - SS_within
    df_between = len(unique_g) - 1
    df_within  = N - len(unique_g)
    F = (SS_between / df_between) / (SS_within / df_within) if SS_within > 0 else np.nan
    # scalar diagnostics
    within_pairs, between_pairs = [], []
    for i in range(N):
        for j in range(i+1, N):
            if groups[i] == groups[j]: within_pairs.append(D[i,j])
            else: between_pairs.append(D[i,j])
    return F, np.mean(within_pairs), np.mean(between_pairs)

def permanova_p(D, groups, n_perm=9999, seed=0):
    F_obs, mw, mb = permanova_F(D, groups)
    rng = np.random.default_rng(seed)
    F_perm = np.empty(n_perm)
    for k in range(n_perm):
        g_p = rng.permutation(groups)
        F_perm[k], _, _ = permanova_F(D, g_p)
    p = (np.sum(F_perm >= F_obs) + 1) / (n_perm + 1)
    return F_obs, p, mw, mb

# run for both distance metrics
results = []
for name, D in [('pearson_corr_dist', D_corr), ('bray_curtis', D_bray)]:
    F, p, mw, mb = permanova_p(D, arms_short, n_perm=9999, seed=0)
    results.append({'metric':name, 'F':round(F,3), 'permanova_p':round(p,4),
                    'mean_within_arm':round(mw,3), 'mean_between_arm':round(mb,3),
                    'within_over_between':round(mw/mb,3) if mb>0 else np.nan})
    print(f'  {name:20s}   F={F:.3f}  perm p={p:.4f}   within={mw:.3f}  between={mb:.3f}  ratio={mw/mb:.3f}')

pd.DataFrame(results).to_csv(OUT/'Q0_permanova_pvalues.csv', index=False)

# save the 25×25 distance matrices
pd.DataFrame(D_corr, index=samples, columns=samples).to_csv(OUT/'Q0_pairwise_phylo_distance_corr.csv')
pd.DataFrame(D_bray, index=samples, columns=samples).to_csv(OUT/'Q0_pairwise_phylo_distance_bray.csv')

# ── figure 1: within vs between violin (2 metrics) ─────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
for ax, name, D in zip(axes, ['pearson_corr_dist','bray_curtis'], [D_corr, D_bray]):
    within, between = [], []
    for i in range(len(samples)):
        for j in range(i+1, len(samples)):
            if np.isnan(D[i,j]): continue
            if arms_short[i] == arms_short[j]:
                within.append(D[i,j])
            else:
                between.append(D[i,j])
    parts = ax.violinplot([within, between], positions=[0,1], widths=0.7,
                          showmeans=False, showmedians=False, showextrema=False)
    for pc, col in zip(parts['bodies'], ['#4a7fbf', '#c44e52']):
        pc.set_facecolor(col); pc.set_edgecolor('none'); pc.set_alpha(0.35)
    for i, (data, col) in enumerate(zip([within, between], ['#4a7fbf','#c44e52'])):
        xj = np.random.default_rng(i).uniform(-0.10, 0.10, size=len(data))
        ax.scatter(np.full(len(data), i) + xj, data, s=25, c='white',
                   edgecolors='#222', linewidths=0.5, zorder=3)
        ax.hlines(np.median(data), i-0.20, i+0.20, colors='#222', lw=1.4, zorder=4)
    ax.set_xticks([0,1])
    ax.set_xticklabels([f'within-arm\nn={len(within)} pairs', f'between-arm\nn={len(between)} pairs'], fontsize=10)
    ax.set_ylabel('phylogenetic distance', fontsize=10)
    r = [r for r in results if r['metric']==name][0]
    verdict = '⚠ arm-level tests INVALID' if r['within_over_between'] >= 1 else 'arm-level tests OK'
    ax.set_title(f'{name}\nF={r["F"]}  PERMANOVA p={r["permanova_p"]}\nwithin/between = {r["within_over_between"]}  ({verdict})',
                 fontsize=10, color='#0f2747')
    ax.grid(axis='y', ls=':', alpha=0.4)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.suptitle('Cohort-heterogeneity — do within-arm tumors share phylogenetic pattern more than between-arm?',
             fontsize=12, fontweight='bold', color='#0f2747', y=1.05)
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'Q0_within_vs_between.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'[saved] {OUT/"Q0_within_vs_between.png"}')

# ── figure 2: 25×25 pairwise heatmap ordered by arm ────────────────────────
order = np.argsort([{'Intact':0,'Cas':1,'Cas+Enz':2}[a] for a in arms_short])
samples_ord = np.array(samples)[order]
arms_ord    = arms_short[order]
S_ord = (1 - D_corr[np.ix_(order, order)])   # SIMILARITY (higher = more similar)
np.fill_diagonal(S_ord, np.nan)                # mask diagonal (self-similarity = 1, uninformative)
fig, ax = plt.subplots(figsize=(10, 9))
# palette: blue-white-red (RdBu_r); higher similarity = redder, lower = bluer.
# Diagonal (self-similarity=1, uninformative) is NaN → color it explicit grey (not white) so it's clearly "no data" not "low similarity"
sim_cmap = plt.get_cmap('RdBu_r').copy()
sim_cmap.set_bad(color='#bdbdbd')   # grey for NaN diagonal
im = ax.imshow(S_ord, cmap=sim_cmap,
               vmin=np.nanpercentile(S_ord, 5),
               vmax=np.nanpercentile(S_ord, 99))
ax.set_xticks(range(len(samples_ord))); ax.set_yticks(range(len(samples_ord)))
ax.set_xticklabels(samples_ord, rotation=90, fontsize=8)
ax.set_yticklabels(samples_ord, fontsize=8)
plt.colorbar(im, ax=ax, fraction=0.03, label='Pearson r (phylogenetic similarity)  ·  higher = more similar')
# arm boundary lines
prev = arms_ord[0]; boundaries = []
for k, a in enumerate(arms_ord):
    if a != prev: boundaries.append(k); prev = a
for b in boundaries:
    ax.axvline(b-0.5, color='k', lw=1.5)
    ax.axhline(b-0.5, color='k', lw=1.5)
# arm labels along diagonal
arm_starts = [0] + boundaries + [len(arms_ord)]
for i in range(len(arm_starts)-1):
    mid = (arm_starts[i] + arm_starts[i+1] - 1) / 2
    ax.text(-1.8, mid, arms_ord[arm_starts[i]], rotation=90, ha='center', va='center',
            fontweight='bold', color=ARM_COLORS[arms_ord[arm_starts[i]]], fontsize=10)
ax.set_title(f'Pairwise phylogenetic similarity (25×25, Pearson r) ordered by arm  ·  diagonal masked\n'
             f'PERMANOVA F={results[0]["F"]}  p={results[0]["permanova_p"]}  within/between distance ratio={results[0]["within_over_between"]}',
             fontsize=11, color='#0f2747', fontweight='bold')
plt.tight_layout()
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'Q0_pairwise_heatmap.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'[saved] {OUT/"Q0_pairwise_heatmap.png"}')

# ── figure 3: 2D MDS embedding of tumors, colored by arm ───────────────────
D_use = D_corr.copy()
D_use[np.isnan(D_use)] = np.nanmean(D_use)   # impute for MDS
mds = MDS(n_components=2, dissimilarity='precomputed', random_state=0)
xy = mds.fit_transform(D_use)
fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
for arm in ['Intact','Cas','Cas+Enz']:
    m = arms_short == arm
    ax.scatter(xy[m,0], xy[m,1], s=110, c=ARM_COLORS[arm], edgecolors='#222',
               linewidths=0.8, alpha=0.85, label=f'{arm} (n={m.sum()})')
    # per-arm centroid + hull
    if m.sum() >= 2:
        cx, cy = xy[m].mean(axis=0)
        ax.scatter(cx, cy, s=280, c=ARM_COLORS[arm], edgecolors='#222',
                   linewidths=1.8, marker='X', zorder=5)
for i, s in enumerate(samples):
    ax.annotate(s, xy[i], fontsize=7, ha='left', va='bottom', color='#333',
                xytext=(3,3), textcoords='offset points')
ax.set_xlabel('MDS-1', fontsize=11); ax.set_ylabel('MDS-2', fontsize=11)
ax.set_title(f'2D MDS of 25 tumors on Fitch-transition-matrix similarity  '
             f'(F={results[0]["F"]}, PERMANOVA p={results[0]["permanova_p"]})',
             fontsize=11, color='#0f2747', fontweight='bold')
ax.legend(fontsize=10, frameon=False, loc='best')
ax.grid(ls=':', alpha=0.4)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'Q0_MDS_embedding.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'[saved] {OUT/"Q0_MDS_embedding.png"}')

# ── console verdict ────────────────────────────────────────────────────────
print('\n=== Q0 VERDICT ===')
for r in results:
    ratio = r['within_over_between']
    verdict = ('WITHIN-ARM MORE DIVERSE than between-arm — arm-level tests INVALID due to heterogeneity'
               if ratio >= 1.0
               else 'within-arm tighter than between-arm — arm-level tests are legitimate')
    print(f'  {r["metric"]:20s}  F={r["F"]}  perm p={r["permanova_p"]}  within/between={ratio}  → {verdict}')

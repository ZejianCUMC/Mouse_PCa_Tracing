#!/usr/bin/env python
"""§10 · Q0 · Cohort-heterogeneity — PATH-Z-score featurization (Schiffman & Landau 2024).

Same question as Q0_cohort_heterogeneity.py (Fitch featurization) but re-run with the
peer-reviewed PATH phylogenetic-autocorrelation Z-score matrix as the phylogenetic signature
per tumor.

Featurization per tumor
-----------------------
- Read {sample}_PATH_zscore.csv from PATH_out/ (per-sample s×s Moran's I Z-score
  standardized against a phylo-permutation null: permute cell labels but keep the tree
  topology + branch lengths + marginal state frequencies fixed).
- Diagonal Z_αα = state heritability (phylo-autocorrelation within-state); off-diagonal
  Z_αβ = state cross-correlation (co-clustering positive / repulsion negative).
- Pad missing-state rows/cols with 0 (null = "no autocorrelation") to obtain a full
  12×12 matrix for every tumor → flatten → 144-dim signature per tumor.

Advantages over Fitch featurization
-----------------------------------
- Moran's I integrates over the WHOLE tree topology (weight = one_node_tree_dist),
  not just parent→child edges (Fitch is myopic).
- Z-standardization against phylo-permutation null makes cross-tree comparison
  unit-consistent: a Z=+2 in a small tumor means the same as Z=+2 in a large tumor.
- Branch-length aware.

Distances used
--------------
- Pearson correlation distance (1 - r)          — matches Fitch Q0 metric
- Cosine distance (1 - cos)                     — natural on signed Z-vectors

Test
----
PERMANOVA F (Anderson 2001) with 9999 arm-label permutations.

Outputs (analysis/evolution_features_v2/lineage_6Q/):
  PATH_zmatrices/{sample}_Zmatrix.csv          25× per-sample 12×12 (missing states padded 0)
  Q0_PATH_pairwise_phylo_distance_corr.csv     25×25
  Q0_PATH_pairwise_phylo_distance_cosine.csv   25×25
  Q0_PATH_permanova_pvalues.csv                F, p, within/between ratio per metric
  Q0_PATH_within_vs_between_pyviz.{png,pdf}    matplotlib backup panel
  Q0_PATH_vs_Fitch_comparison.csv              side-by-side stats
  Q0_PATH_pair_long.csv                        long-format pairs (for R ggpubr plotting)
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform, pdist
from scipy.stats import pearsonr

# ── paths ─────────────────────────────────────────────────────────────────────
M     = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
PATH_ = M/'analysis/evolution_features_v2/PATH_out'
SCFP  = M/'analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT   = M/'analysis/evolution_features_v2/lineage_6Q'
ZDIR  = OUT/'PATH_zmatrices'
ZDIR.mkdir(exist_ok=True)
FITCH_PERM_CSV = OUT/'Q0_permanova_pvalues.csv'   # reference from Fitch Q0

STATES = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
          'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
          'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
NS = len(STATES)
ARM_ORDER = ['Intact','Castration','Castration + Enzalutamide']
ARM_SHORT = {'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}

# ── sample -> arm map ─────────────────────────────────────────────────────────
scfp = pd.read_csv(SCFP)
sample_arm = scfp.groupby('sample')['treatment'].first().to_dict()
samples = sorted([s for s in sample_arm if (PATH_/f'{s}_PATH_zscore.csv').exists()])
print(f'[load] {len(samples)} samples with PATH Z-score CSVs')
arms_short = np.array([ARM_SHORT[sample_arm[s]] for s in samples])
print(f'  arms: {dict(zip(*np.unique(arms_short, return_counts=True)))}')

# ── build 12×12 padded PATH Z-matrix per sample ───────────────────────────────
print('[build] per-sample 12×12 PATH Z-score matrices (missing states padded 0)...')
per_sample_Z = {}
per_sample_covered = {}
for s in samples:
    Z_sub = pd.read_csv(PATH_/f'{s}_PATH_zscore.csv', index_col=0)
    # sanity: Z_sub should be square with states as row+col names
    assert list(Z_sub.index) == list(Z_sub.columns), f'{s} zscore CSV not square/aligned'
    covered = list(Z_sub.index)
    per_sample_covered[s] = covered
    # build 12×12, pad missing with 0
    Z12 = pd.DataFrame(0.0, index=STATES, columns=STATES)
    # only fill cells where BOTH row & col state are present in this tumor
    shared = [x for x in covered if x in STATES]
    Z12.loc[shared, shared] = Z_sub.loc[shared, shared].values
    # replace NaN (PATH can return NaN for singleton states) with 0
    Z12 = Z12.fillna(0.0)
    per_sample_Z[s] = Z12
    # write per-sample padded CSV
    Z12.to_csv(ZDIR/f'{s}_Zmatrix.csv')
print(f'  saved {len(samples)} 12×12 Z-matrices to {ZDIR}/')

# state coverage summary
coverage_df = pd.DataFrame({
    'sample': samples,
    'arm': arms_short,
    'n_states_covered': [len(per_sample_covered[s]) for s in samples],
})
coverage_df.to_csv(ZDIR/'PATH_zmatrix_state_coverage.csv', index=False)
print(f'  state-coverage distribution: '
      f'min={coverage_df["n_states_covered"].min()}, '
      f'median={int(coverage_df["n_states_covered"].median())}, '
      f'max={coverage_df["n_states_covered"].max()}')

# ── flatten each 12×12 → 144-dim ──────────────────────────────────────────────
X = np.array([per_sample_Z[s].values.flatten() for s in samples])   # 25 × 144
print(f'[flatten] X shape = {X.shape}')

# ── compute pairwise distances (2 metrics) ────────────────────────────────────
def corr_distance(X):
    """1 − Pearson r between rows; ignores samples with zero variance."""
    N = X.shape[0]
    D = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j: D[i, j] = 0; continue
            xi, xj = X[i], X[j]
            if xi.std() == 0 or xj.std() == 0: D[i, j] = np.nan
            else: D[i, j] = 1 - pearsonr(xi, xj)[0]
    return D

D_corr   = corr_distance(X)
D_cosine = squareform(pdist(X, metric='cosine'))
print(f'  D_corr range   [{np.nanmin(D_corr):.3f}, {np.nanmax(D_corr):.3f}]')
print(f'  D_cosine range [{np.nanmin(D_cosine):.3f}, {np.nanmax(D_cosine):.3f}]')

# ── PERMANOVA (Anderson 2001 pseudo-F) ────────────────────────────────────────
def permanova_F(D, groups):
    """Anderson 2001 pseudo-F. Returns (F, mean_within, mean_between)."""
    mask = ~np.isnan(D).any(axis=1)
    D = D[np.ix_(mask, mask)]; groups = groups[mask]
    N = D.shape[0]
    unique_g = np.unique(groups)
    SS_within = 0.0
    for g in unique_g:
        idx = np.where(groups == g)[0]
        if len(idx) < 2: continue
        SS_within += np.sum(D[np.ix_(idx, idx)]**2) / (2 * len(idx))
    SS_total = np.sum(D**2) / (2 * N)
    SS_between = SS_total - SS_within
    df_between = len(unique_g) - 1
    df_within  = N - len(unique_g)
    F = (SS_between / df_between) / (SS_within / df_within) if SS_within > 0 else np.nan
    within_pairs, between_pairs = [], []
    for i in range(N):
        for j in range(i+1, N):
            if groups[i] == groups[j]: within_pairs.append(D[i, j])
            else: between_pairs.append(D[i, j])
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
print('[PERMANOVA] running 9999 permutations for each metric...')
results = []
for name, D in [('pearson_corr_dist', D_corr), ('cosine_dist', D_cosine)]:
    F, p, mw, mb = permanova_p(D, arms_short, n_perm=9999, seed=0)
    ratio = mw/mb if mb > 0 else np.nan
    results.append({'metric': name, 'F': round(F, 3), 'permanova_p': round(p, 4),
                    'mean_within_arm': round(mw, 4), 'mean_between_arm': round(mb, 4),
                    'within_over_between': round(ratio, 3) if not np.isnan(ratio) else np.nan})
    print(f'  {name:20s}   F={F:.3f}  perm p={p:.4f}   within={mw:.3f}  '
          f'between={mb:.3f}  ratio={ratio:.3f}')

pd.DataFrame(results).to_csv(OUT/'Q0_PATH_permanova_pvalues.csv', index=False)
pd.DataFrame(D_corr,   index=samples, columns=samples).to_csv(OUT/'Q0_PATH_pairwise_phylo_distance_corr.csv')
pd.DataFrame(D_cosine, index=samples, columns=samples).to_csv(OUT/'Q0_PATH_pairwise_phylo_distance_cosine.csv')

# ── long-format pair table for R ggpubr plotting ──────────────────────────────
print('[export] pair table for R ggpubr...')
pair_rows = []
for name, D in [('pearson_corr_dist', D_corr), ('cosine_dist', D_cosine)]:
    for i in range(len(samples)):
        for j in range(i+1, len(samples)):
            if np.isnan(D[i, j]): continue
            kind = 'within-arm' if arms_short[i] == arms_short[j] else 'between-arm'
            pair_rows.append({
                'metric': name, 'kind': kind,
                'sample_i': samples[i], 'sample_j': samples[j],
                'arm_i': arms_short[i], 'arm_j': arms_short[j],
                'distance': D[i, j], 'similarity': 1 - D[i, j],
            })
pair_df = pd.DataFrame(pair_rows)
pair_df.to_csv(OUT/'Q0_PATH_pair_long.csv', index=False)

# ── side-by-side comparison vs Fitch ──────────────────────────────────────────
print('[compare] Fitch vs PATH side-by-side...')
fitch = pd.read_csv(FITCH_PERM_CSV)
# align on Pearson metric for the primary claim
fitch_row = fitch.iloc[0].to_dict()   # first row is pearson_corr_dist
path_row  = results[0]                # first is pearson_corr_dist

cmp = pd.DataFrame([
    {'method':'Fitch (parent→child transition counts)',
     'featurization':'12×12 count matrix (flattened 144-dim)',
     'distance_metric':'Pearson correlation distance (1-r)',
     'F': fitch_row['F'], 'permanova_p': fitch_row['permanova_p'],
     'mean_within_arm': fitch_row['mean_within_arm'],
     'mean_between_arm': fitch_row['mean_between_arm'],
     'within_over_between': fitch_row['within_over_between'],
     'arm_signal_significant_p<0.01': fitch_row['permanova_p'] < 0.01,
     'within_tighter_than_between': fitch_row['within_over_between'] < 1.0},
    {'method':'PATH (Moran-I phylo-autocorrelation Z-matrix)',
     'featurization':'12×12 Z-score matrix (flattened 144-dim, missing-state → 0 null)',
     'distance_metric':'Pearson correlation distance (1-r)',
     'F': path_row['F'], 'permanova_p': path_row['permanova_p'],
     'mean_within_arm': path_row['mean_within_arm'],
     'mean_between_arm': path_row['mean_between_arm'],
     'within_over_between': path_row['within_over_between'],
     'arm_signal_significant_p<0.01': path_row['permanova_p'] < 0.01,
     'within_tighter_than_between': path_row['within_over_between'] < 1.0},
    # extra second-metric rows for reference
    {'method':'Fitch (parent→child transition counts)',
     'featurization':'12×12 count matrix (flattened 144-dim)',
     'distance_metric':'Bray-Curtis on row-normalized counts',
     'F': fitch.iloc[1]['F'], 'permanova_p': fitch.iloc[1]['permanova_p'],
     'mean_within_arm': fitch.iloc[1]['mean_within_arm'],
     'mean_between_arm': fitch.iloc[1]['mean_between_arm'],
     'within_over_between': fitch.iloc[1]['within_over_between'],
     'arm_signal_significant_p<0.01': fitch.iloc[1]['permanova_p'] < 0.01,
     'within_tighter_than_between': fitch.iloc[1]['within_over_between'] < 1.0},
    {'method':'PATH (Moran-I phylo-autocorrelation Z-matrix)',
     'featurization':'12×12 Z-score matrix (flattened 144-dim, missing-state → 0 null)',
     'distance_metric':'Cosine distance (natural on signed Z)',
     'F': results[1]['F'], 'permanova_p': results[1]['permanova_p'],
     'mean_within_arm': results[1]['mean_within_arm'],
     'mean_between_arm': results[1]['mean_between_arm'],
     'within_over_between': results[1]['within_over_between'],
     'arm_signal_significant_p<0.01': results[1]['permanova_p'] < 0.01,
     'within_tighter_than_between': results[1]['within_over_between'] < 1.0},
])
cmp.to_csv(OUT/'Q0_PATH_vs_Fitch_comparison.csv', index=False)
print(cmp[['method','distance_metric','F','permanova_p','within_over_between',
           'arm_signal_significant_p<0.01','within_tighter_than_between']].to_string(index=False))

# ── concordance verdict ──────────────────────────────────────────────────────
fitch_sig  = fitch_row['permanova_p'] < 0.01 and fitch_row['within_over_between'] < 1.0
path_sig   = path_row['permanova_p']  < 0.01 and path_row['within_over_between']  < 1.0
if fitch_sig and path_sig:
    verdict = 'AGREE — both methods show significant within-arm > between-arm similarity (p<0.01, ratio<1). Arm signal is method-independent.'
elif (not fitch_sig) and (not path_sig):
    verdict = 'AGREE — both methods show NO arm signal (p>=0.01 or ratio>=1). Arm-level tests are effectively noise under both featurizations.'
else:
    verdict = f'DISAGREE — Fitch(sig={fitch_sig}) vs PATH(sig={path_sig}). Investigate WHY: Fitch is myopic (parent-child only), PATH integrates whole tree; heterogeneity may be branch-length- or long-range-topology-dependent.'
print(f'\n[verdict] {verdict}')

# ── matplotlib backup 2-panel figure (Fitch vs PATH) ─────────────────────────
print('[plot] matplotlib backup 2-panel figure...')
# reload Fitch pairs from CSV
fitch_dcorr = pd.read_csv(OUT/'Q0_pairwise_phylo_distance_corr.csv', index_col=0)
fitch_sim = 1 - fitch_dcorr.values     # Pearson r as similarity
path_sim  = 1 - D_corr                 # PATH Pearson r similarity

def pairs_by_kind(sim, samples_list, arms_arr):
    within, between = [], []
    n = len(samples_list)
    for i in range(n):
        for j in range(i+1, n):
            v = sim[i, j]
            if np.isnan(v): continue
            (within if arms_arr[i] == arms_arr[j] else between).append(v)
    return within, between

fitch_samples = list(fitch_dcorr.index)
fitch_arms = np.array([ARM_SHORT[sample_arm[s]] for s in fitch_samples])
fitch_within, fitch_between = pairs_by_kind(fitch_sim, fitch_samples, fitch_arms)
path_within,  path_between  = pairs_by_kind(path_sim,  samples, arms_short)

fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), constrained_layout=True)
for ax, (data_w, data_b), title, F, p in [
    (axes[0], (fitch_within, fitch_between),
     f'Fitch parent→child\nF={fitch_row["F"]}, p={fitch_row["permanova_p"]}, '
     f'ratio(sim W/B)={np.median(fitch_within)/np.median(fitch_between):.2f}',
     fitch_row['F'], fitch_row['permanova_p']),
    (axes[1], (path_within, path_between),
     f'PATH Moran-I Z-matrix\nF={path_row["F"]}, p={path_row["permanova_p"]}, '
     f'ratio(sim W/B)={np.median(path_within)/np.median(path_between):.2f}',
     path_row['F'], path_row['permanova_p']),
]:
    parts = ax.violinplot([data_w, data_b], positions=[0, 1], widths=0.75,
                          showmeans=False, showmedians=False, showextrema=False)
    for pc, col in zip(parts['bodies'], ['#4a7fbf','#c44e52']):
        pc.set_facecolor(col); pc.set_edgecolor('none'); pc.set_alpha(0.45)
    # inner white boxplot (width=0.1), no dots per user preference
    bp = ax.boxplot([data_w, data_b], positions=[0, 1], widths=0.10,
                    patch_artist=True, showfliers=False)
    for box in bp['boxes']:
        box.set(facecolor='white', edgecolor='#222', linewidth=0.8)
    for element in ['whiskers','medians','caps']:
        for item in bp[element]: item.set(color='#222', linewidth=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'within-arm\nn={len(data_w)} pairs',
                        f'between-arm\nn={len(data_b)} pairs'], fontsize=10)
    ax.set_ylabel('phylogenetic similarity (higher = more similar)', fontsize=10)
    ax.set_title(title, fontsize=10.5, color='#0f2747', fontweight='bold')
    ax.grid(axis='y', ls=':', alpha=0.4)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.suptitle('Q0 cohort heterogeneity — Fitch vs PATH featurization  '
             f'(concordance: {"AGREE" if fitch_sig==path_sig else "DISAGREE"})',
             fontsize=12, fontweight='bold', color='#0f2747')
for fmt in ['png','pdf']:
    fig.savefig(OUT/f'Q0_PATH_within_vs_between_pyviz.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'[saved] {OUT/"Q0_PATH_within_vs_between_pyviz.png"}')

# ── console final summary ────────────────────────────────────────────────────
print('\n=== Q0-PATH FINAL ===')
print(f'  Fitch  (Pearson dist): F={fitch_row["F"]}, p={fitch_row["permanova_p"]}, ratio={fitch_row["within_over_between"]}')
print(f'  PATH   (Pearson dist): F={path_row["F"]}, p={path_row["permanova_p"]}, ratio={path_row["within_over_between"]}')
print(f'  PATH   (Cosine dist ): F={results[1]["F"]}, p={results[1]["permanova_p"]}, ratio={results[1]["within_over_between"]}')
print(f'  Verdict: {verdict}')

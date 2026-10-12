#!/usr/bin/env python
"""§10 evolution — per-sample state-specificity, three scalars, three arms.

Question: does treatment shift where a mouse's tumor sits on the
"idiosyncratic vs canonical state repertoire" axis?

Three per-sample scalars (all computed on the same state-proportion vector p_i
of shape 12 = the 12 putative_states in sc_fitness_plasticity_v10.csv):

  (1) f_rare_state_cells
      = fraction of the mouse's cells in cohort-rare states
      rare state = one whose universality tier is *worse* than "shared-majority"
      i.e. present at >=10 cells in FEWER than 12/25 mice
      (matches the §10.1 tier gate — c4, c8, c9, c10, c11, c12 are "rare")

  (2) KL_from_cohort
      = KL(p_mouse || p_cohort_mean)
      where p_cohort_mean = mean of p_i across all 25 mice
      Laplace-smoothed with eps = 1e-6 so log(0) does not fire

  (3) shannon_states  (base 2) + effective_n_states = 2**shannon_states
      = H_2(p_mouse) over 12 states — repertoire breadth
      low = mono-lineaged (few active states); high = broad repertoire

Then per-arm distribution + 3-arm Kruskal-Wallis + pairwise Welch t-tests.

Input : sc_fitness_plasticity_v10.csv (columns: sample, treatment, putative_state)
Output: sample_specificity_per_arm.{pdf,png,csv}  in the same folder.
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import kruskal, ttest_ind

# ---- config ---------------------------------------------------------------
CELL_FLOOR      = 10      # state is "present in a sample" if >=10 cells (matches §10.1)
UBIQUITY_GATE   = 12      # state is "rare" if present in <12/25 samples (i.e. worse than shared-majority)
EPS             = 1e-6    # laplace smoothing for KL

STATES = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
          'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
          'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
STATE_DISP = {s: f'c{i+1}_{s}' for i,s in enumerate(STATES)}

ARM_ORDER  = ['Intact', 'Castration', 'Castration + Enzalutamide']
ARM_SHORT  = {'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}
ARM_COLORS = {'Intact':'#4a7fbf','Castration':'#e0955f','Castration + Enzalutamide':'#c44e52'}

M   = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
IN  = M/'analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT = M/'analysis/evolution_features_v2'
OUT.mkdir(exist_ok=True)

# ---- 1. build per-mouse state-count matrix -------------------------------
scfp = pd.read_csv(IN)
scfp = scfp[scfp.putative_state.isin(STATES)].copy()
print(f'loaded {len(scfp):,} cells across {scfp["sample"].nunique()} mice / {scfp["putative_state"].nunique()} states')

# sample -> arm
sample_arm = scfp.groupby('sample')['treatment'].first()
sample_arm = sample_arm.reindex(sample_arm.index.tolist())

count_mat  = (scfp.groupby(['sample','putative_state']).size()
                .unstack(fill_value=0).reindex(columns=STATES, fill_value=0))
count_mat  = count_mat.reindex(sample_arm.index)
print(f'count matrix: {count_mat.shape}')

# ---- 2. cohort-wide state universality (matches §10.1) --------------------
n_samples_with_state = (count_mat >= CELL_FLOOR).sum(axis=0)
rare_states = n_samples_with_state[n_samples_with_state < UBIQUITY_GATE].index.tolist()
print(f'\ncohort universality (n samples with >={CELL_FLOOR} cells):')
for s in STATES:
    tier = 'RARE' if s in rare_states else 'shared'
    print(f'  {STATE_DISP[s]:24s}  n_samples={int(n_samples_with_state[s]):2d} / 25  [{tier}]')
print(f'\n{len(rare_states)} rare states -> {rare_states}')

# ---- 3. per-sample scalars -----------------------------------------------
tot            = count_mat.sum(axis=1).clip(lower=1)
n_rare_cells   = count_mat[rare_states].sum(axis=1)
f_rare         = (n_rare_cells / tot).rename('f_rare_state_cells')

# proportion vector p_i
prop           = count_mat.div(tot, axis=0)
prop_smooth    = (prop + EPS)
prop_smooth    = prop_smooth.div(prop_smooth.sum(axis=1), axis=0)
p_cohort_mean  = prop_smooth.mean(axis=0)

kl_vals = (prop_smooth * np.log2(prop_smooth.div(p_cohort_mean, axis=1))).sum(axis=1)
kl_vals.name = 'KL_from_cohort'

# Shannon (base 2)
def shannon(p):
    p = p[p > 0]
    return -(p * np.log2(p)).sum()
shannon_vals = prop.apply(shannon, axis=1).rename('shannon_states')   # bits
eff_n_states = (2**shannon_vals).rename('effective_n_states')

per_sample = pd.DataFrame({
    'sample'          : count_mat.index,
    'treatment'       : sample_arm.values,
    'n_cells_total'   : tot.values.astype(int),
    'n_cells_rare'    : n_rare_cells.values.astype(int),
    'f_rare_state_cells' : f_rare.values,
    'KL_from_cohort'  : kl_vals.values,
    'shannon_states'  : shannon_vals.values,
    'effective_n_states' : eff_n_states.values,
})
per_sample.to_csv(OUT/'sample_specificity_per_arm.csv', index=False)
print(f'\n[saved] {OUT/"sample_specificity_per_arm.csv"}')

# ---- 4. per-arm stats -----------------------------------------------------
def summarise(col):
    rows = []
    grouped = {a: per_sample[per_sample.treatment == a][col].dropna().values for a in ARM_ORDER}
    for a in ARM_ORDER:
        v = grouped[a]
        rows.append({'arm': ARM_SHORT[a], 'n': len(v),
                     'mean': np.mean(v), 'median': np.median(v),
                     'sd': np.std(v, ddof=1) if len(v)>1 else np.nan})
    stat, p_kw = kruskal(*[grouped[a] for a in ARM_ORDER])
    pairs = []
    for i,ai in enumerate(ARM_ORDER):
        for aj in ARM_ORDER[i+1:]:
            t, p = ttest_ind(grouped[ai], grouped[aj], equal_var=False)
            pairs.append({'pair': f'{ARM_SHORT[ai]} vs {ARM_SHORT[aj]}',
                          't': t, 'p_welch': p})
    return pd.DataFrame(rows), p_kw, pd.DataFrame(pairs)

print('\n' + '='*80)
scalar_cols = ['f_rare_state_cells','KL_from_cohort','shannon_states','effective_n_states']
stat_records = {}
for c in scalar_cols:
    summ, pkw, pairs = summarise(c)
    stat_records[c] = (summ, pkw, pairs)
    print(f'\n--- {c} ---')
    print(summ.to_string(index=False))
    print(f'  KW p = {pkw:.3g}')
    print(pairs.to_string(index=False))

# ---- 5. figure -----------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
def plot_panel(ax, col, ylab, title):
    data = [per_sample[per_sample.treatment == a][col].dropna().values for a in ARM_ORDER]
    parts = ax.violinplot(data, positions=range(len(ARM_ORDER)), widths=0.75,
                          showmeans=False, showmedians=False, showextrema=False)
    for pc, a in zip(parts['bodies'], ARM_ORDER):
        pc.set_facecolor(ARM_COLORS[a]); pc.set_edgecolor('none'); pc.set_alpha(0.35)
    for i,a in enumerate(ARM_ORDER):
        v = data[i]
        xj = np.random.default_rng(i).uniform(-0.10, 0.10, size=len(v))
        ax.scatter(np.full(len(v), i) + xj, v, s=32, c=ARM_COLORS[a],
                   edgecolors='#222', linewidths=0.5, zorder=3)
        # median bar
        ax.hlines(np.median(v), i-0.25, i+0.25, colors='#222', lw=1.3, zorder=4)
    ax.set_xticks(range(len(ARM_ORDER)))
    ax.set_xticklabels([ARM_SHORT[a] + f'\nn={len(d)}' for a,d in zip(ARM_ORDER,data)], fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)
    _, p_kw, pairs = stat_records[col]
    ax.set_title(f'{title}\nKW p={p_kw:.3g}', fontsize=11, color='#0f2747', loc='left')
    # pairwise p annotations at top
    ymax = max([np.max(d) for d in data if len(d)]); ymin = min([np.min(d) for d in data if len(d)])
    yr = ymax - ymin if ymax > ymin else max(abs(ymax), 1e-6)
    yoff = ymax + 0.05*yr
    for k, r in pairs.iterrows():
        p_str = f'p={r.p_welch:.2g}'
        idx_map = {'Intact vs Cas':(0,1), 'Intact vs Cas+Enz':(0,2), 'Cas vs Cas+Enz':(1,2)}
        i,j = idx_map[r['pair']]
        y = yoff + k*0.06*yr
        ax.plot([i,j], [y,y], color='#555', lw=0.8)
        ax.text((i+j)/2, y + 0.008*yr, p_str, ha='center', fontsize=8, color='#333')
    ax.set_ylim(ymin - 0.05*yr, yoff + 3*0.06*yr + 0.05*yr)
    ax.grid(axis='y', ls=':', alpha=0.4)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)

plot_panel(axes[0], 'f_rare_state_cells',
           f'Fraction of cells in rare states\n(rare = coverage < {UBIQUITY_GATE}/25)',
           'A · Rare-state occupancy')
plot_panel(axes[1], 'KL_from_cohort',
           'KL(p_mouse || p_cohort_mean)  [bits]',
           'B · State-composition atypicality')
plot_panel(axes[2], 'effective_n_states',
           'effective # states  (2^Shannon)',
           'C · Repertoire breadth')

fig.suptitle('§10 · Per-sample state-specificity across treatment arms  (n=25 mice; 14 Intact / 8 Cas / 3 Cas+Enz)',
             fontsize=12, fontweight='bold', color='#0f2747', y=1.05)
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'sample_specificity_per_arm.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'\n[saved] {OUT/"sample_specificity_per_arm.png"}')
print(f'[saved] {OUT/"sample_specificity_per_arm.pdf"}')

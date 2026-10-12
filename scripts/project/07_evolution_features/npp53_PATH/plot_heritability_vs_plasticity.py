#!/usr/bin/env python
"""Per-state scatter: PATH heritability (x) vs Fitch plasticity (y).
Heritability = mean diagonal of PATH xcor phylocorr matrix across all
25 samples where the state is present.
Plasticity   = mean per-cell Fitch-parsimony plasticity across all cells
of that state in the cohort."""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

BASE = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
PATH_DIR = BASE / 'analysis/evolution_features_v2/PATH_out'
SCFP     = BASE / 'analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT      = BASE / 'analysis/evolution_features_v2/carta_fine12/biology_review'

STATE_ORDER = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
               'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
               'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
SHORT = {s: f'c{i+1}' for i, s in enumerate(STATE_ORDER)}
STATE_COLORS = {
    'Luminal-Distal':'#1f5fa8','Luminal-Proximal':'#9467bd','Luminal-ARdep':'#d23030',
    'Luminal-CR':'#f1903d','Luminal-Persistent':'#117c70','Luminal-Squamous':'#f78ec7',
    'Basal-like':'#8c564b','Partial-EMT':'#bcbd22','Mes-Quiescent':'#c49c94',
    'Mes-Collagen':'#5b5b5b','Mes-Stemness':'#8c1a1a','Mes-myCAF':'#4daf4a',
}

# ── (1) Aggregate PATH heritability across samples ───────────────────────
her_vals = {s: [] for s in STATE_ORDER}
samples = sorted(set(p.stem.replace('_PATH_phylocorr','') for p in PATH_DIR.glob('*_PATH_phylocorr.csv')))
print(f'{len(samples)} sample phylocorr CSVs')
for samp in samples:
    csv = PATH_DIR / f'{samp}_PATH_phylocorr.csv'
    df = pd.read_csv(csv, index_col=0)
    for s in df.index:
        if s in her_vals:
            her_vals[s].append(df.loc[s, s])

her_summary = []
for s in STATE_ORDER:
    v = her_vals[s]
    her_summary.append({
        'state': s, 'short': SHORT[s], 'n_samples': len(v),
        'her_mean': np.mean(v) if v else np.nan,
        'her_std':  np.std(v)  if len(v) > 1 else 0.0,
    })
her_df = pd.DataFrame(her_summary)
print('\n=== Heritability per state ===')
print(her_df.to_string(index=False))

# ── (2) Aggregate plasticity per state ───────────────────────────────────
scfp = pd.read_csv(SCFP)
scfp = scfp[scfp['putative_state'].isin(STATE_ORDER)].copy()
plast_summary = scfp.groupby('putative_state').agg(
    plast_mean=('plasticity','mean'),
    plast_std=('plasticity','std'),
    n_cells=('plasticity','count'),
).reset_index().rename(columns={'putative_state':'state'})
print('\n=== Plasticity per state ===')
print(plast_summary.to_string(index=False))

# ── (3) Merge + rank correlation ─────────────────────────────────────────
df = her_df.merge(plast_summary, on='state')
df = df.dropna(subset=['her_mean','plast_mean']).reset_index(drop=True)
rp_stat, rp_p = stats.pearsonr(df['her_mean'], df['plast_mean'])
rs_stat, rs_p = stats.spearmanr(df['her_mean'], df['plast_mean'])
print(f'\nn states: {len(df)}')
print(f'Pearson r  = {rp_stat:+.3f}  (p={rp_p:.3f})')
print(f'Spearman ρ = {rs_stat:+.3f}  (p={rs_p:.3f})')

df.to_csv(OUT/'npp53_heritability_vs_plasticity.tsv', sep='\t', index=False)

# ── (4) Scatter ──────────────────────────────────────────────────────────
# Per-state label offset directions curated to avoid neighbor + quadrant-label collisions.
# Values: (dx, dy, ha, va) where dx/dy are in data units.
LABEL_OFFSETS = {
    'Luminal-Distal':     (+0.055, +0.008, 'left',  'center'), # c1 — push far RIGHT (dot at 0.262, label starts at 0.317)
    'Luminal-Proximal':   (+0.020, +0.020, 'left',  'bottom'), # c2 — up-right, +0.008 gap from c6 dot
    'Luminal-ARdep':      (-0.008, -0.022, 'right', 'top'),    # c3 — small left drop; keeps text on-canvas + off y-title
    'Luminal-CR':         (+0.030, -0.020, 'left',  'top'),    # c4 — DOWN-right so it sits below dot, clears c12
    'Luminal-Persistent': (+0.020, -0.030, 'left',  'top'),    # c5 — DOWN-right, clears c3 dot
    'Luminal-Squamous':   (+0.014, +0.035, 'left',  'bottom'), # c6 — up-RIGHT (above c2 by ~0.008) — stays inside xlim=[0, ...]
    'Basal-like':         (-0.030, -0.028, 'right', 'top'),    # c7 — push DOWN-LEFT (dot at 0.219, label ends ~0.189) clears c1
    'Partial-EMT':        (+0.020, -0.028, 'left',  'top'),    # c8 — down-right, Δy vs c12 = 0.038 (safe)
    'Mes-Quiescent':      (+0.020, +0.010, 'left',  'center'), # c9 — push RIGHT into empty margin (dot at 0.366, label at 0.386-0.469) — clears c4 label
    'Mes-Collagen':       (-0.014, +0.022, 'right', 'bottom'), # c10 — up-left
    'Mes-Stemness':       (+0.014, -0.024, 'left',  'top'),    # c11 — down-right
    'Mes-myCAF':          (+0.032, +0.028, 'left',  'bottom'), # c12 — up-right, ~0.06 above c4 to clear label collision
}

fig, ax = plt.subplots(figsize=(12.5, 8.5))
for _, r in df.iterrows():
    x, y = r['her_mean'], r['plast_mean']
    c = STATE_COLORS[r['state']]
    ax.errorbar(x, y, xerr=r['her_std']/np.sqrt(max(1, r['n_samples'])),
                yerr=r['plast_std']/np.sqrt(max(1, r['n_cells'])),
                ecolor=c, elinewidth=1.2, capsize=3, alpha=0.55, zorder=2)
    ax.scatter([x], [y], s=320, c=c, edgecolor='#111', linewidth=1.0, zorder=3)
    ax.annotate(r['short'], (x, y), fontsize=10, fontweight='bold',
                color='white', ha='center', va='center', zorder=4)
    dx, dy, ha, va = LABEL_OFFSETS[r['state']]
    label = f"{r['short']} {r['state']}\nn={int(r['n_cells']):,} cells · {int(r['n_samples'])} samples"
    ax.annotate(label, (x, y), xytext=(x + dx, y + dy),
                fontsize=8.5, color=c, fontweight='bold', ha=ha, va=va, zorder=5,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, pad=1.5))

# Quadrant reference lines at cohort means
ax.axvline(df['her_mean'].mean(), color='#888', ls='--', lw=0.8, alpha=0.5)
ax.axhline(df['plast_mean'].mean(), color='#888', ls='--', lw=0.8, alpha=0.5)

# FIXED axis limits per user request: x=[0, 0.6], y=[0, 0.4]. This leaves the
# corner quadrant labels in the empty margins well away from the state clusters
ax.set_xlim(0.0, 0.6)
ax.set_ylim(0.0, 0.4)
xlo, xhi = 0.0, 0.6; ylo, yhi = 0.0, 0.4

# Quadrant corner labels in the four TRUE corners of the fixed [0,0.6]×[0,0.4]
# box — the state cluster occupies only ~0.0-0.4 × 0.05-0.30, so the extreme
# corners are safe empty-margin real estate.
xr = 0.012; yr = 0.015
for (x_pos, y_pos, ha, va, text) in [
    (xhi - xr, yhi - yr, 'right', 'top',    'HIGH heritability\nHIGH plasticity'),
    (xlo + xr, yhi - yr, 'left',  'top',    'LOW heritability\nHIGH plasticity'),
    (xhi - xr, ylo + yr, 'right', 'bottom', 'HIGH heritability\nLOW plasticity'),
    (xlo + xr, ylo + yr, 'left',  'bottom', 'LOW heritability\nLOW plasticity'),
]:
    ax.text(x_pos, y_pos, text, ha=ha, va=va, fontsize=9, style='italic',
             color='#555', alpha=0.95, fontweight='bold',
             bbox=dict(facecolor='white', edgecolor='#ccc', alpha=0.9,
                        boxstyle='round,pad=0.4'))

ax.set_xlabel('PATH heritability  (mean sibling-pair phylo-corr across samples where state is present)',
              fontsize=11, fontweight='bold')
ax.set_ylabel('Fitch plasticity  (mean per-cell subclade state-change rate)',
              fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.grid(True, alpha=0.25)

fig.suptitle('Per-state PATH heritability × Fitch plasticity  ·  Npp53 cohort (25 samples)',
             fontsize=13.5, fontweight='bold', color='#0f2747', y=0.995)
fig.text(0.5, -0.01,
         f'12 CARTA states  ·  x = PATH xcor diagonal (state autocorrelation on sibling-leaf pairs)  ·  y = Fitch parsimony subclade transition rate.  '
         f'Cross-state rank correlation: Spearman ρ = {rs_stat:+.3f} (p = {rs_p:.3f}).  Error bars = SEM.',
         ha='center', fontsize=9, color='#555', style='italic')
plt.subplots_adjust(top=0.90, bottom=0.12, left=0.10, right=0.96)
for fmt in ['png','pdf']:
    fig.savefig(OUT/f'npp53_heritability_vs_plasticity_scatter.{fmt}', dpi=200,
                 bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'\nsaved: {OUT}/npp53_heritability_vs_plasticity_scatter.png')

#!/usr/bin/env python
"""Per-sample CARTA best-k selection — v2 (final: search range k=1..8).

Selection method: **maximum perpendicular distance from chord** on the
(k, normalized-discrepancy) Pareto curve. This is the formal version of what
CARTA authors do visually — "identify the elbow of the Pareto front" — and
is robust to the trivial k=1 baseline that dominates kneedle-default.

Algorithm (per sample):
  1. Read objective for k=1..K where K = min(n_states, 8)
     (k > n_states is a duplicate plateau — do not use for selection)
  2. Normalize discrepancy y_norm = obj / n_leaves  (CARTA paper's convention)
  3. Normalize both axes to [0,1] range
  4. Chord = line from (k_min, y_norm_min) to (k_max, y_norm_max)
  5. For each intermediate k, compute perpendicular distance from (k, y_norm(k)) to chord
  6. Best k* = arg max perpendicular distance

The chord-distance method has been the recommended way to detect the "knee"
in Pareto-optimal fronts since Satopää et al. 2011 (kneedle uses it too, with
extra smoothing) but without kneedle's convexity constraint and S-sensitivity
tunable, both of which bias k selection at low state counts.

Compares against a battery of alternate methods for reporting only.
Also shows kneedle for legacy comparison.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    from kneed import KneeLocator
except ImportError:
    KneeLocator = None

M = '/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
R = Path(f'{M}/analysis/evolution_features_v2/carta_fine12/results_fine12_greedy')
OUT = Path(f'{M}/analysis/evolution_features_v2/carta_fine12')

K_MIN, K_MAX = 1, 8

arm_scfp = pd.read_csv(f'{M}/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv')
arm_of = arm_scfp.groupby('sample')['treatment'].first().replace(
    {'Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}).to_dict()
ARM_COLORS = {'Intact':'#4a7fbf','Cas':'#e0955f','Cas+Enz':'#c44e52'}
ARM_ORDER  = ['Intact','Cas','Cas+Enz']

summary = pd.read_csv(OUT/'fine12_greedy_summary.tsv', sep='\t')
n_states = summary.set_index('sample')['n_states'].to_dict()
n_leaves = summary.set_index('sample')['n_leaves'].to_dict()

def read_curve(sample):
    rows = []
    for k in range(K_MIN, K_MAX + 1):
        f = R/sample/f'{sample}_results_{k}.txt'
        if not f.exists(): continue
        df = pd.read_csv(f, sep='\t')
        rows.append({'k': int(df.iloc[0]['k']), 'obj': float(df.iloc[0]['objective_score'])})
    d = pd.DataFrame(rows).drop_duplicates('k').sort_values('k').reset_index(drop=True)
    d['obj_norm'] = d['obj'] / n_leaves[sample]
    # plateau: obj not decreasing anymore (within tolerance) → k > n_states typically
    d['is_plateau'] = False
    for i in range(1, len(d)):
        if d['obj'].iloc[i] >= d['obj'].iloc[i-1] - 1e-9:
            d.loc[i:, 'is_plateau'] = True
            break
    return d

def chord_distance_elbow(xs, ys):
    """Perpendicular distance from each point to the chord (first, last).
    Returns the x that maximizes the distance (the mathematical elbow).
    Both axes normalized to [0,1] before distance computation."""
    if len(xs) < 3: return xs[0] if xs else np.nan
    xs = np.array(xs, dtype=float); ys = np.array(ys, dtype=float)
    # normalize to [0,1]
    xn = (xs - xs.min()) / max(1e-12, xs.max() - xs.min())
    yn = (ys - ys.min()) / max(1e-12, ys.max() - ys.min())
    # chord from (0,1) to (1,0)  (curve is decreasing, so left-top → right-bottom)
    # perpendicular distance from (xn[i], yn[i]) to that chord
    # chord vector = (1, -1); line equation: x + y = 1  →  distance = |x+y-1| / sqrt(2)
    d = np.abs(xn + yn - 1) / np.sqrt(2)
    # elbow = point with maximum distance ABOVE the chord (i.e., x+y < 1 for a convex-decreasing curve)
    # but for a monotonically-decreasing curve after normalization, the maximum d corresponds
    # to the "most convex" point, which IS the elbow
    return int(xs[int(np.argmax(d))])

def kneedle(xs, ys, S=1.0):
    if KneeLocator is None or len(xs) < 3: return xs[0]
    try:
        kl = KneeLocator(xs, ys, S=S, curve='convex', direction='decreasing')
        return int(kl.elbow) if kl.elbow is not None else xs[-1]
    except Exception:
        return xs[-1]

# --- Per-sample selection ---
rows = []
curves = {}
for samp in sorted(n_states.keys()):
    d = read_curve(samp)
    curves[samp] = d
    # Use only NON-plateau points for elbow (plateau tail distorts geometry)
    d_use = d[~d['is_plateau']].reset_index(drop=True)
    xs = d_use['k'].tolist()
    ys = d_use['obj_norm'].tolist()
    k_chord   = chord_distance_elbow(xs, ys)
    k_kneedle = kneedle(xs, ys, S=1.0)
    n_st = n_states[samp]
    rows.append({
        'sample': samp, 'arm': arm_of[samp], 'n_states': n_st, 'n_leaves': n_leaves[samp],
        'k_max_informative': int(d_use['k'].max()),
        'k_star_chord':   k_chord,     # <-- FINAL — chord-distance elbow
        'k_kneedle_S1':   k_kneedle,   # legacy comparison
    })
tab = pd.DataFrame(rows).sort_values(['arm','sample']).reset_index(drop=True)
tab.to_csv(OUT/'best_k_per_sample_v2_chord.csv', index=False)

print('=' * 90)
print('Best-k per sample — CHORD-DISTANCE method (final, search k=1..8)')
print('=' * 90)
print(tab.to_string(index=False))

print('\n=== Distribution of k_star_chord ===')
for k in sorted(tab['k_star_chord'].unique()):
    ss = sorted(tab.loc[tab['k_star_chord']==k, 'sample'].tolist())
    print(f'  k*={k} ({len(ss)} samples): {ss}')

# --- 25-panel plot with chord-elbow highlighted ---
samples_ordered = sorted(n_states.keys(), key=lambda s: (ARM_ORDER.index(arm_of[s]), s))
fig, axes = plt.subplots(5, 5, figsize=(16, 13), constrained_layout=True)
for ax, samp in zip(axes.flat, samples_ordered):
    d = curves[samp]; arm = arm_of[samp]; color = ARM_COLORS[arm]
    d_inf = d[~d['is_plateau']]; d_pl = d[d['is_plateau']]
    ax.plot(d_inf['k'], d_inf['obj_norm'], 'o-', color=color, lw=1.7, ms=5, zorder=3, label='informative')
    if len(d_pl):
        bridge_k = list(d_inf['k'].tail(1)) + list(d_pl['k'])
        bridge_y = list(d_inf['obj_norm'].tail(1)) + list(d_pl['obj_norm'])
        ax.plot(bridge_k, bridge_y, ':', color='#888', lw=1.2, zorder=2)
        ax.plot(d_pl['k'], d_pl['obj_norm'], 'o', color='#888', ms=4, mfc='white', zorder=2, label='plateau')
    # Chord line (from first to last informative point)
    if len(d_inf) >= 2:
        ax.plot([d_inf['k'].iloc[0], d_inf['k'].iloc[-1]],
                [d_inf['obj_norm'].iloc[0], d_inf['obj_norm'].iloc[-1]],
                '-', color='#bbb', lw=0.7, zorder=1)
    row = tab[tab['sample']==samp].iloc[0]
    kc = row['k_star_chord']
    y_at = d.loc[d['k']==kc, 'obj_norm'].values
    if len(y_at):
        ax.scatter([kc], [y_at[0]], s=140, edgecolor='#111', facecolor='#ffde3a',
                   lw=1.5, marker='*', zorder=6, label=f'chord k*={kc}')
    kn = row['k_kneedle_S1']
    if kn != kc:
        y_kn = d.loc[d['k']==kn, 'obj_norm'].values
        if len(y_kn):
            ax.scatter([kn], [y_kn[0]], s=60, edgecolor='#c44e52', facecolor='none',
                       lw=1.2, marker='o', zorder=5, label=f'kneedle k={kn}')
    ax.set_title(f'{samp} ({arm})  n_st={n_states[samp]}  k*={kc}',
                 fontsize=9, color='#0f2747', fontweight='bold', pad=3)
    ax.set_xticks([1,2,3,4,5,6,7,8])
    ax.set_xlabel('k', fontsize=8); ax.set_ylabel('obj / n_leaves', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.suptitle("CARTA fine12 · normalized objective vs k (k=1..8)",
             fontsize=13, fontweight="bold", color="#0f2747")
fig.text(0.5, 0.98, "★ = chord-distance elbow (k*)    ○ = kneedle S=1    — = chord line", ha="center", fontsize=11, color="#444", style="italic")
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'carta_fine12_obj_curves_v2.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'\nsaved: {OUT}/carta_fine12_obj_curves_v2.png')

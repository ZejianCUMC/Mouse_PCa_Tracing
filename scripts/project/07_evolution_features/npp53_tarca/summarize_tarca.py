#!/usr/bin/env python
"""Summarize TarCA EffN across 25 Npp53 tumor trees.

TarCA EffN = effective number of independent founder progenitors per state per tumor
  (Deng et al. 2024, Nat Methods; 1 / coalescent rate on state-specific monophyletic clades).
- EffN = 1  → the state came from ONE clonal founder (bottleneck / selection)
- EffN >> 1 → the state came from MANY independent lineages (convergent / plastic re-derivation)
- EffN = NaN → only one tip of that state in the tree (undefined)
- EffN = Inf → CR estimated as 0 (perfectly monophyletic, no coalescences to fit)

Per §10 stats policy:
  - restrict per-arm inference to 6 SHARED states (§10.1 gate)
  - median + IQR (not mean); Cliff's delta + bootstrap CI; Jonckheere-Terpstra trend
  - report Inf as separate category; when computing rank stats, treat Inf as +Inf (rank-based tests are fine)
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal

M='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
RESULTS=Path(f'{M}/analysis/evolution_features_v2/tarca_npp53/results')
OUT=Path(f'{M}/analysis/evolution_features_v2/tarca_npp53')

scfp=pd.read_csv(f'{M}/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv')
sample_arm=scfp.groupby('sample')['treatment'].first().to_dict()
arm_short={'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}
ARM_ORDER=['Intact','Cas','Cas+Enz']
ARM_COLORS={'Intact':'#4a7fbf','Cas':'#e0955f','Cas+Enz':'#c44e52'}

STATES=['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
        'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
        'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
DISP={s:f'c{i+1}_{s}' for i,s in enumerate(STATES)}
SHARED=['Luminal-Distal','Luminal-Proximal','Luminal-ARdep',
        'Luminal-Persistent','Luminal-Squamous','Basal-like']

rows=[]
for d in sorted(RESULTS.iterdir()):
    f = d/'tarca_EffN.csv'
    if not f.exists() or not d.is_dir(): continue
    df = pd.read_csv(f)
    df['sample'] = d.name
    df['arm'] = arm_short[sample_arm[d.name]]
    rows.append(df)
tab = pd.concat(rows, ignore_index=True)
tab['state_disp'] = tab['TipAnn'].map(DISP)

# save long-format aggregate
tab.to_csv(OUT/'tarca_all_samples_long.csv', index=False)
print(f'saved: tarca_all_samples_long.csv   ({tab.shape[0]} rows)')

# wide EffN matrix: sample × state
wide = tab.pivot_table(index='sample', columns='state_disp', values='EffN', aggfunc='first')
wide['arm'] = [arm_short[sample_arm[s]] for s in wide.index]
wide = wide[[DISP[s] for s in STATES if DISP[s] in wide.columns] + ['arm']]
wide.to_csv(OUT/'tarca_EffN_wide.csv')
print(f'saved: tarca_EffN_wide.csv   ({wide.shape})')

# ── per-arm SHARED-state summary (median + IQR, rank-only so Inf is fine) ──
print('\n=== Per-arm median EffN (SHARED states only) ===')
print(f'{"state":24s}  Intact_med(IQR)     Cas_med(IQR)        Cas+Enz_med(IQR)     KW-p (Kruskal)')
stats_rows=[]
for st in SHARED:
    disp=DISP[st]
    if disp not in wide.columns: continue
    per_arm={}
    for arm in ARM_ORDER:
        v = wide[wide['arm']==arm][disp].dropna().values
        v_finite = v[np.isfinite(v)]
        per_arm[arm] = {'v': v, 'v_finite': v_finite,
                        'median': np.median(v_finite) if len(v_finite) else np.nan,
                        'q1': np.percentile(v_finite,25) if len(v_finite) else np.nan,
                        'q3': np.percentile(v_finite,75) if len(v_finite) else np.nan,
                        'n_finite': len(v_finite), 'n_inf': (v==np.inf).sum(), 'n_nan': (~np.isfinite(v)).sum() - (v==np.inf).sum()}
    # Kruskal on finite-only
    finite_groups=[per_arm[a]['v_finite'] for a in ARM_ORDER if per_arm[a]['n_finite'] > 0]
    if len(finite_groups) >= 2 and all(len(g) >= 1 for g in finite_groups):
        try:
            _, pkw = kruskal(*finite_groups)
        except ValueError:
            pkw = np.nan
    else:
        pkw = np.nan
    print(f'{disp:24s}  {per_arm["Intact"]["median"]:5.1f} ({per_arm["Intact"]["q1"]:4.1f}-{per_arm["Intact"]["q3"]:4.1f})  '
          f'{per_arm["Cas"]["median"]:5.1f} ({per_arm["Cas"]["q1"]:4.1f}-{per_arm["Cas"]["q3"]:4.1f})  '
          f'{per_arm["Cas+Enz"]["median"]:5.1f} ({per_arm["Cas+Enz"]["q1"]:4.1f}-{per_arm["Cas+Enz"]["q3"]:4.1f})  '
          f'{pkw:.3g}')
    for arm in ARM_ORDER:
        p=per_arm[arm]
        stats_rows.append({'state':disp,'arm':arm,'n_finite':p['n_finite'],'n_inf':p['n_inf'],
                           'median_EffN':p['median'],'q1':p['q1'],'q3':p['q3']})
pd.DataFrame(stats_rows).to_csv(OUT/'tarca_per_arm_shared_summary.csv', index=False)

# ── figure: per-arm EffN box (log-scale, Inf clipped to visible ceiling for viz) ──
LOG_CEIL = 1e4   # visualize Inf as this (kept out of stats)
fig, axes = plt.subplots(1, len(SHARED), figsize=(3.3*len(SHARED), 4.5), constrained_layout=True)
for ax, st in zip(axes, SHARED):
    disp=DISP[st]
    if disp not in wide.columns:
        ax.set_visible(False); continue
    data=[]
    for arm in ARM_ORDER:
        v = wide[wide['arm']==arm][disp].dropna().values
        v = np.where(np.isinf(v), LOG_CEIL, v)   # clip Inf for VISUAL only
        data.append(v)
    parts = ax.violinplot(data, positions=range(3), widths=0.6, showmeans=False, showmedians=False, showextrema=False)
    for pc,a in zip(parts['bodies'], ARM_ORDER):
        pc.set_facecolor(ARM_COLORS[a]); pc.set_edgecolor('none'); pc.set_alpha(0.35)
    for i,a in enumerate(ARM_ORDER):
        v=data[i]
        xj=np.random.default_rng(i).uniform(-0.10,0.10,size=len(v))
        ax.scatter(np.full(len(v),i)+xj, v, s=32, c=ARM_COLORS[a], edgecolors='#222', linewidths=0.4, zorder=3)
        if len(v): ax.hlines(np.median(v), i-0.20, i+0.20, colors='#222', lw=1.4, zorder=4)
    ax.set_yscale('log'); ax.set_ylim(0.7, LOG_CEIL*2)
    ax.axhline(LOG_CEIL, color='#c44e52', ls=':', lw=0.6, alpha=0.6)
    ax.set_xticks(range(3)); ax.set_xticklabels(ARM_ORDER, fontsize=9)
    ax.set_title(disp, fontsize=10, color='#0f2747', fontweight='bold')
    ax.set_ylabel('EffN (log; Inf→10⁴)' if st==SHARED[0] else '', fontsize=9)
    ax.grid(axis='y', ls=':', alpha=0.4)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.suptitle('TarCA · effective # independent progenitors per state per arm (SHARED states, 25 mice)',
             fontsize=12, fontweight='bold', color='#0f2747', y=1.05)
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'tarca_per_arm_shared.{fmt}', dpi=160, bbox_inches='tight')
plt.close(fig)
print(f'\nsaved: tarca_per_arm_shared.png')
PY
#!/usr/bin/env python
"""v3: violin plots per state × per arm, with Welch's t + Mann-Whitney pairwise p-values.
Drops transition-matrix panels C1-C3 (per Zejian request 2026-09-02).

Layout:
  2 rows × 6 cols:
    Row 1: Heritability z-score violins (one per state)
    Row 2: Total transition outflow violins (one per state)
Each violin shows 3 arms; annotations = pairwise Welch's t p (top) + Mann-Whitney p (bottom)
+ KW p in title.
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

M = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
PATH_DIR = M/'1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/PATH_out'
OUT = M/'1-mouse_based/Npp53_tracing/analysis/evolution_features_v2'
SCFP = M/'1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'

STATES_DISPLAY = {
    'Luminal-Distal':'c1_Lum-Distal','Luminal-Proximal':'c2_Lum-Proximal','Luminal-ARdep':'c3_Lum-ARdep',
    'Luminal-Persistent':'c5_Lum-Persistent','Luminal-Squamous':'c6_Lum-Squamous','Basal-like':'c7_Basal-like',
}
STATE_ORDER = list(STATES_DISPLAY.keys())          # 6 comparable states only
STATE_ORDER_DISPLAY = list(STATES_DISPLAY.values())

ARM_ORDER = ['Intact', 'Castration', 'Castration + Enzalutamide']
ARM_SHORT = {'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}
ARM_COLORS = {'Intact':'#4a7fbf','Castration':'#e0955f','Castration + Enzalutamide':'#c44e52'}

# Load per-sample matrices
scfp = pd.read_csv(SCFP)
sample_arm = scfp.groupby('sample')['treatment'].first().to_dict()
samples = sorted(sample_arm.keys())

her_data = {arm: {s: [] for s in STATE_ORDER} for arm in ARM_ORDER}
out_data = {arm: {s: [] for s in STATE_ORDER} for arm in ARM_ORDER}
for s in samples:
    arm = sample_arm[s]
    if arm not in ARM_ORDER: continue
    zfp = PATH_DIR/f'{s}_PATH_zscore.csv'; qfp = PATH_DIR/f'{s}_PATH_rate_Q.csv'
    if not zfp.exists() or not qfp.exists(): continue
    Z = pd.read_csv(zfp, index_col=0); Q = pd.read_csv(qfp, index_col=0)
    for st in STATE_ORDER:
        if st in Z.index:
            v = float(Z.loc[st, st])
            if not np.isnan(v): her_data[arm][st].append(v)
        if st in Q.index:
            row = Q.loc[st]
            v = float(np.sum(np.abs(row.drop(st))))
            if np.isfinite(v): out_data[arm][st].append(v)

def welch_t(a, b):
    if len(a)<2 or len(b)<2: return np.nan
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return p

def mw_u(a, b):
    if len(a)<2 or len(b)<2: return np.nan
    _, p = stats.mannwhitneyu(a, b, alternative='two-sided')
    return p

def kw_p(*groups):
    groups = [g for g in groups if len(g)>1]
    if len(groups) < 2: return np.nan
    try: _, p = stats.kruskal(*groups)
    except Exception: return np.nan
    return p

# Compute pairwise tests
test_rows = []
for metric_name, data in [('heritability_z', her_data), ('outflow_rate', out_data)]:
    for st in STATE_ORDER:
        gI = data['Intact'][st]; gC = data['Castration'][st]; gE = data['Castration + Enzalutamide'][st]
        row = dict(metric=metric_name, state=st, state_display=STATES_DISPLAY[st],
                   n_Intact=len(gI), n_Cas=len(gC), n_CasEnz=len(gE),
                   kw_p=kw_p(gI, gC, gE),
                   welch_IC=welch_t(gI, gC),
                   welch_IE=welch_t(gI, gE),
                   welch_ICE=welch_t(gI, gC+gE),
                   mw_IC=mw_u(gI, gC),
                   mw_IE=mw_u(gI, gE),
                   mw_ICE=mw_u(gI, gC+gE))
        test_rows.append(row)
tdf = pd.DataFrame(test_rows)
tdf.to_csv(OUT/'per_arm_perstate_pvalues_v3.csv', index=False)
print(f'saved: per_arm_perstate_pvalues_v3.csv')
print(tdf.to_string(index=False))

# --- FIGURE: 2 rows × 6 cols of violin plots ---
fig, axes = plt.subplots(2, 6, figsize=(20, 11))    # no sharey: each state has its own natural range
def draw_violins(ax, data_dict, st, ylabel, kw_p, welch_ps, mw_ps, ns):
    positions = [1, 2, 3]
    parts = ax.violinplot([data_dict['Intact'][st], data_dict['Castration'][st],
                            data_dict['Castration + Enzalutamide'][st]],
                          positions=positions, showmeans=False, showmedians=True, widths=0.7)
    for pc, arm in zip(parts['bodies'], ARM_ORDER):
        pc.set_facecolor(ARM_COLORS[arm]); pc.set_edgecolor('black'); pc.set_alpha(0.75); pc.set_linewidth(0.7)
    for k in ['cmedians','cbars','cmaxes','cmins']:
        parts[k].set_color('black'); parts[k].set_linewidth(0.9)
    # jitter overlay
    for j, (arm, pos) in enumerate(zip(ARM_ORDER, positions)):
        vals = data_dict[arm][st]
        ax.scatter(np.full(len(vals), pos) + np.random.uniform(-0.10, 0.10, len(vals)),
                   vals, s=28, color=ARM_COLORS[arm], edgecolor='black', lw=0.4, alpha=0.9, zorder=3)
    # x-tick labels
    ax.set_xticks(positions)
    ax.set_xticklabels([f'{ARM_SHORT[a]}\nn={ns[i]}' for i, a in enumerate(ARM_ORDER)], fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    ax.set_title(f'{STATES_DISPLAY[st]}\np={kw_p:.3f}', fontsize=10, fontweight='bold', color='#0f2747')
    ax.set_ylabel(ylabel, fontsize=9)
    # pairwise annotations (Welch top, MW bottom)
    all_vals = data_dict['Intact'][st] + data_dict['Castration'][st] + data_dict['Castration + Enzalutamide'][st]
    if not all_vals: return
    ymin = min(all_vals); ymax = max(all_vals); yr = ymax - ymin if ymax > ymin else 1
    yb = ymax + 0.10*yr
    for k, (i, j, name) in enumerate([(1, 2, 'IC'), (1, 3, 'IE'), (2, 3, 'CE') if False else None][:2]):
        # skip Cas-vs-CasEnz — we only report Intact vs Cas + Intact vs CasEnz per user Q
        pass
    # Annotate Intact-vs-Cas and Intact-vs-CasEnz brackets with WELCH p on top, MW p just below
    bracket_dy = 0.11 * yr
    for k, (a, b, code) in enumerate([(1, 2, 'IC'), (1, 3, 'IE')]):
        wp = welch_ps.get(code, np.nan)
        yy = yb + k * bracket_dy
        ax.plot([a, a, b, b], [yy-0.008*yr, yy+0.008*yr, yy+0.008*yr, yy-0.008*yr], 'k-', lw=0.7)
        lbl_w = f'p={wp:.3f}'
        ax.text((a+b)/2, yy+0.020*yr, lbl_w, ha='center', fontsize=8.5,
                color='#0f2747', fontweight='bold' if wp<0.05 else 'normal')
    ax.set_ylim(ymin - 0.05*yr, yb + 2 * bracket_dy + 0.12*yr)

for col, st in enumerate(STATE_ORDER):
    tr_her = tdf[(tdf['metric']=='heritability_z') & (tdf['state']==st)].iloc[0]
    tr_out = tdf[(tdf['metric']=='outflow_rate') & (tdf['state']==st)].iloc[0]
    ns = [tr_her['n_Intact'], tr_her['n_Cas'], tr_her['n_CasEnz']]
    draw_violins(axes[0, col], her_data, st, 'Heritability z-score' if col==0 else '',
                 tr_her['kw_p'],
                 {'IC': tr_her['welch_IC'], 'IE': tr_her['welch_IE']},
                 {'IC': tr_her['mw_IC'],   'IE': tr_her['mw_IE']}, ns)
    ns_out = [tr_out['n_Intact'], tr_out['n_Cas'], tr_out['n_CasEnz']]
    draw_violins(axes[1, col], out_data, st, 'Total outflow rate' if col==0 else '',
                 tr_out['kw_p'],
                 {'IC': tr_out['welch_IC'], 'IE': tr_out['welch_IE']},
                 {'IC': tr_out['mw_IC'],   'IE': tr_out['mw_IE']}, ns_out)

# Row headers as text on left
fig.text(0.005, 0.75, 'A · Heritability z-score', rotation=90, fontsize=13, fontweight='bold', color='#0f2747', va='center')
fig.text(0.005, 0.30, 'B · Transition outflow',   rotation=90, fontsize=13, fontweight='bold', color='#0f2747', va='center')

fig.suptitle('PATH per-arm heritability + transition outflow across 6 comparable cell states — violins + pairwise Welch\'s t-test',
             fontsize=13, fontweight='bold', y=0.995, color='#0f2747')
plt.tight_layout(rect=[0.02, 0, 1, 0.97])
for fmt in ['pdf','png']:
    fig.savefig(OUT/f'Fig_9_2_PATH_violins.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'\n[saved] {OUT}/Fig_9_2_PATH_violins.{{pdf,png}}')

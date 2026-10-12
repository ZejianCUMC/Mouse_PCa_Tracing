#!/usr/bin/env python
"""pySCENIC downstream — per-state top regulons on the 18,162 tumor cells.

Per-state top-regulon workflow (mirrors karthaus_v3clean_pyscenic pattern):
  1. Load auc_mtx.csv (cells × regulons)
  2. Join to sc_fitness_plasticity_v10.csv (cellBC → putative_state, treatment)
  3. Per state: mean AUCell per regulon → z-score across states
  4. Top-10 regulons per state (by state-z) → table + heatmap
  5. Bonus: top-15 arm-differential regulons (Intact vs Cas vs Cas+Enz)
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

M = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
AUC   = M/'pysenic/downstream/auc_mtx.csv'
SCFP  = M/'analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT   = M/'pysenic/downstream'

STATES = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
          'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
          'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
DISP  = {s: f'c{i+1}_{s}' for i,s in enumerate(STATES)}

print('[load]')
auc = pd.read_csv(AUC, index_col=0)
scfp = pd.read_csv(SCFP)
print(f'  auc_mtx: {auc.shape} (cells × regulons)')
print(f'  scfp   : {len(scfp):,} cells')

# Regulon names in pySCENIC output are like "Ar(+)" / "Foxa1(+)" — strip the "(+)" for readability
auc.columns = [c.replace('(+)','').strip() for c in auc.columns]
print(f'  n regulons after cleanup: {auc.shape[1]}')

# ---- align by cellBC ----
# auc_mtx.csv rows are keyed by "<orig.ident>_<barcode>-1" (h5ad-style).
# scfp uses (sample=JZ*, cellBC=16nt). Bridge via bare 16-nt barcode + sample-id via h5ad obs.
import anndata as ad
H5 = M/'data/npp53_merge/npp53_edited.h5ad'
a  = ad.read_h5ad(H5, backed='r')
bridge = a.obs[['sample_id']].copy()
bridge['bc16']  = bridge.index.str.split('_').str[1].str.replace('-1','', regex=False)
bridge['scfp_key'] = bridge['sample_id'].astype(str) + '::' + bridge['bc16']
bridge = bridge.reset_index().rename(columns={'index':'auc_key'})
# scfp_key on scfp side
scfp['scfp_key'] = scfp['sample'].astype(str) + '::' + scfp['cellBC'].astype(str)

# Merge auc.index → scfp_key → putative_state/treatment
auc = auc.reset_index().rename(columns={'index':'auc_key'})
merged = auc.merge(bridge[['auc_key','scfp_key']], on='auc_key', how='inner') \
            .merge(scfp[['scfp_key','putative_state','treatment']], on='scfp_key', how='inner')
print(f'  merged shape: {merged.shape}   ({merged.shape[0]:,} cells with regulon + state)')
if merged.shape[0] < 15000:
    raise SystemExit(f'  ⚠ merge too small — abort')

# reset auc to the aligned matrix + carry state/treatment as sidecar Series
regulon_cols = [c for c in merged.columns if c not in ('auc_key','scfp_key','putative_state','treatment')]
auc  = merged[regulon_cols].copy()
auc.index = merged['auc_key'].values
scfp = merged[['putative_state','treatment']].copy()
scfp.index = merged['auc_key'].values
print(f'  aligned auc: {auc.shape}')

# ---- per-state mean AUCell ----
scfp['state'] = scfp['putative_state'].map(DISP)
per_state = auc.groupby(scfp['state']).mean().T   # regulon × state
per_state = per_state[[DISP[s] for s in STATES if DISP[s] in per_state.columns]]
per_state.to_csv(OUT/'per_state_mean_AUCell.csv')

# ---- z-score across states (per regulon) — 'specificity' ----
per_state_z = per_state.sub(per_state.mean(axis=1), axis=0).div(per_state.std(axis=1) + 1e-9, axis=0)
per_state_z.to_csv(OUT/'per_state_z_AUCell.csv')

# ---- top-10 regulons per state (by state-z) ----
top_per_state = {}
lines = []
for state in per_state_z.columns:
    top = per_state_z[state].sort_values(ascending=False).head(10)
    top_per_state[state] = top
    lines.append(f'\n{state:24s}   top-10 regulons (z-across-states):')
    for reg, z in top.items():
        mean_auc = per_state.loc[reg, state]
        lines.append(f'   {reg:20s}   z={z:+5.2f}   mean_AUCell={mean_auc:.4f}')

top_report = '\n'.join(lines)
(OUT/'top10_regulons_per_state.md').write_text('# Top-10 regulons per tumor state (mean AUCell z-scored across states)\n' + top_report)
print(top_report[:3200])
print('...')

# ---- heatmap: union of top-5 per state × 12 states ----
top5_set = set()
for state in per_state_z.columns:
    top5_set |= set(per_state_z[state].sort_values(ascending=False).head(5).index)
htm = per_state_z.loc[sorted(top5_set)]
htm.to_csv(OUT/'top5_regulons_per_state_z_heatmap.csv')
print(f'\ntop-5-union heatmap: {htm.shape}')

# render heatmap — clip z at [-3, 3], hierarchical clustering on BOTH rows + cols
import seaborn as sns
vmax = 3
htm_clip = htm.clip(lower=-vmax, upper=vmax)
g = sns.clustermap(
    htm_clip,
    cmap='RdBu_r', vmin=-vmax, vmax=vmax, center=0,
    figsize=(9.5, max(10, len(htm)*0.22)),
    method='average', metric='correlation',      # correlation dist = co-varying regulons cluster
    row_cluster=True, col_cluster=True,
    dendrogram_ratio=(0.09, 0.08),
    cbar_kws={'label':'AUCell z-score across states (clipped ±3)'},
    cbar_pos=(0.02, 0.85, 0.03, 0.10),
    xticklabels=True, yticklabels=True,
    linewidths=0, linecolor='none',
)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10)
g.ax_heatmap.tick_params(axis='y', labelsize=8)
g.ax_heatmap.set_xlabel(''); g.ax_heatmap.set_ylabel('')
g.fig.suptitle(f'pySCENIC · top-5 regulons per state × 12 tumor states  ({htm.shape[0]} regulons, union) — clustered both axes',
               fontsize=12, fontweight='bold', color='#0f2747', y=1.00)
for fmt in ['pdf','png']:
    g.fig.savefig(OUT/f'top5_regulons_per_state_heatmap.{fmt}', dpi=160, bbox_inches='tight')
plt.close(g.fig)
print(f'[saved] {OUT/"top5_regulons_per_state_heatmap.png"}  (row+col hclust, method=average, metric=correlation)')

# ---- per-arm mean (bonus) ----
arm_short = {'Intact':'Intact','Castration':'Cas','Castration + Enzalutamide':'Cas+Enz'}
scfp['arm'] = scfp['treatment'].map(arm_short)
per_arm = auc.groupby(scfp['arm']).mean().T
per_arm = per_arm[['Intact','Cas','Cas+Enz']]
per_arm_z = per_arm.sub(per_arm.mean(axis=1), axis=0).div(per_arm.std(axis=1)+1e-9, axis=0)
top_by_arm_diff = (per_arm.max(axis=1) - per_arm.min(axis=1)).sort_values(ascending=False).head(15)
print(f'\n=== top 15 regulons by arm-range (Intact / Cas / Cas+Enz) ===')
for reg, diff in top_by_arm_diff.items():
    row = per_arm.loc[reg]
    print(f'   {reg:22s}  I={row["Intact"]:.3f}  Cas={row["Cas"]:.3f}  Cas+Enz={row["Cas+Enz"]:.3f}   range={diff:.3f}')
per_arm.to_csv(OUT/'per_arm_mean_AUCell.csv')
per_arm_z.to_csv(OUT/'per_arm_z_AUCell.csv')
print(f'\n[saved] {OUT/"per_arm_mean_AUCell.csv"}')
print(f'[saved] {OUT/"top10_regulons_per_state.md"}')

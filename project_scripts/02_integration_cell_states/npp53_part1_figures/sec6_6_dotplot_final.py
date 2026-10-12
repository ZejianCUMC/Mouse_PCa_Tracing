#!/usr/bin/env python
# =============================================================================
# Fig 6.6 · Final tumor marker dotplot
# 12 states (Cycling + Ionocyte removed); top-5 rank_genes_groups markers per state.
# No GFP/mCherry reporter channels; no Ionocyte markers (Foxi1, Ascl3, Cftr).
# =============================================================================
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import scanpy as sc, pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

M   = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
OUT = M/'1-mouse_based/Npp53_tracing/figures/npp53_plots'
STATE_ORDER = [
    'Luminal-Distal', 'Luminal-Proximal', 'Luminal-ARdep', 'Luminal-CR',
    'Luminal-Persistent', 'Luminal-Squamous',
    'Basal-like', 'Partial-EMT',
    'Mes-Quiescent', 'Mes-Collagen', 'Mes-Stemness', 'Mes-myCAF'
]  # 12 states — Cycling and Ionocyte excluded
EXCLUDE_STATES = {'Cycling', 'Ionocyte'}
EXCLUDE_MARKERS = {'Gfp','GFP','mCherry','Mcherry','Foxi1','Ascl3','Cftr','Atp6v1b1','Atp6v0d2'}

print('[1] Load h5ad + filter cells')
a = sc.read_h5ad(M/'1-mouse_based/Npp53_tracing/data/npp53_merge/npp53_edited.h5ad')
print(f'  full: {a.n_obs} cells')
keep = ~a.obs['putative_state'].isin(EXCLUDE_STATES)
a_f = a[keep].copy()
print(f'  after excluding {EXCLUDE_STATES}: {a_f.n_obs} cells')
a_f.obs['putative_state'] = pd.Categorical(a_f.obs['putative_state'], categories=STATE_ORDER, ordered=True)

print('[2] Ensure log-normalized layer + rank_genes_groups')
if 'log1p' not in a_f.uns:
    # canonical: normalize then log1p
    sc.pp.normalize_total(a_f, target_sum=1e4)
    sc.pp.log1p(a_f)

print('[3] rank_genes_groups (Wilcoxon) per state')
sc.tl.rank_genes_groups(a_f, groupby='putative_state', method='wilcoxon', n_genes=50)
top_by_state = {}
for st in STATE_ORDER:
    names = pd.DataFrame(a_f.uns['rank_genes_groups']['names'])[st].tolist()
    # keep first N passing exclusion + not-noise + exists in var
    keep = []
    for g in names:
        if g in EXCLUDE_MARKERS: continue
        if any(g.startswith(p) for p in ['Rps','Rpl','Mrps','Mrpl','mt-','Gm','MT-']): continue
        if g.endswith('Rik'): continue
        if g not in a_f.var_names: continue
        keep.append(g)
        if len(keep) >= 5: break
    top_by_state[st] = keep
    print(f'  {st}: {keep}')

# de-dup across states while preserving order per group
seen = set(); marker_dict = {}
for st in STATE_ORDER:
    ug = [g for g in top_by_state[st] if g not in seen]
    seen.update(ug)
    marker_dict[st] = ug
n_markers = sum(len(v) for v in marker_dict.values())
print(f'[4] {n_markers} unique markers total')

print('[5] Draw dotplot')
dp = sc.pl.dotplot(a_f, marker_dict, groupby='putative_state',
                   dendrogram=False, swap_axes=False,
                   var_group_positions=None, var_group_labels=None,
                   standard_scale='var', return_fig=True, show=False,
                   cmap='Reds', figsize=(max(14, n_markers*0.22), 4.5))
dp.savefig(OUT/'sec6_6_dotplot_final_12states.png', dpi=170, bbox_inches='tight')
dp.savefig(OUT/'sec6_6_dotplot_final_12states.pdf', bbox_inches='tight')
plt.close('all')
print(f'saved: {OUT}/sec6_6_dotplot_final_12states.{{png,pdf}}')
print(f'\nstate counts: {a_f.obs["putative_state"].value_counts().to_dict()}')

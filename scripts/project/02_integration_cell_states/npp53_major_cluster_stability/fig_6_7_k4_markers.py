#!/usr/bin/env python
"""Update Fig 6.7 — enriched marker-gene evidence for the k=4 major cell-type grouping.

Panels:
  A · Top-5 rank_genes_groups markers × 4 major types (dot plot)
  B · Same markers across the 12 fine states (proves the 4 majors emerge from markers)
  C · Canonical prostate lineage markers × 4 major types (independent QC)

Uses the exact k=4 Method-4 grouping (chooseR stability 0.996) that has already
been validated in the major_cluster_stability analysis.
"""
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

M = '/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
OUT = Path(f'{M}/analysis/evolution_features_v2/major_cluster_stability')

FINE_TO_K4 = {
    'Luminal-Distal':'DistalAR_Lum','Luminal-ARdep':'DistalAR_Lum','Luminal-CR':'DistalAR_Lum',
    'Luminal-Proximal':'ProxBasalSq','Luminal-Persistent':'ProxBasalSq',
    'Luminal-Squamous':'ProxBasalSq','Basal-like':'ProxBasalSq',
    'Mes-Collagen':'CoreMes','Mes-Stemness':'CoreMes','Mes-myCAF':'CoreMes',
    'Partial-EMT':'EMT_Bridge','Mes-Quiescent':'EMT_Bridge',
}
K4_ORDER   = ['DistalAR_Lum','ProxBasalSq','EMT_Bridge','CoreMes']
# canonical Npp53 c-numbering (Fig 6.6 convention — FIXED, position = cluster ID)
CANONICAL_ORDER = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
                    'Luminal-Persistent','Luminal-Squamous','Basal-like',
                    'Partial-EMT','Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
CID_PREFIX  = {s: f'c{i+1}_{s}' for i, s in enumerate(CANONICAL_ORDER)}
# Display order — group same-k4 states contiguously (per Zejian): DistalAR (c1,c3,c4) → ProxBasalSq (c2,c5,c6,c7) → EMT_Bridge (c8,c9) → CoreMes (c10,c11,c12)
FINE_ORDER = ['Luminal-Distal','Luminal-ARdep','Luminal-CR',
              'Luminal-Proximal','Luminal-Persistent','Luminal-Squamous','Basal-like',
              'Partial-EMT','Mes-Quiescent',
              'Mes-Collagen','Mes-Stemness','Mes-myCAF']

# Expanded canonical prostate-lineage markers — added crucial per-state markers Zejian flagged as missing
# Sources: Karthaus 2020, Crowley 2020, Xin 2013, Kirk 2020, Yu 2026 (bladder Xenium), MJZ pipeline memory
CANONICAL = {
    'Distal AR-lum (c1/c3/c4)':  ['Krt8','Krt18','Cldn3','Cldn7','Foxa1','Ar','Nkx3-1','Fkbp5','Klk3',
                                   'Pbsn','Tgm4','Nectin4','Elf3','Klk4','Piezo2','Wfdc2','Clu'],
    'Prox/Basal/Squamous (c2/c5/c6/c7)':
                                  ['Ly6d','Perp','Fxyd3','Sfn','Krt5','Krt14','Trp63','Snai2',
                                   'Krt6a','Krt15','Ppp1r1b','Msmb','Sbp','Aqp3','Krt19'],
    'EMT bridge (c8/c9)':         ['Igfbp2','Vim','Spp1','Lgals1','Cd63','Zeb1','Zeb2','Snai1',
                                   'Fn1','Cdh2','Nes','Cd44','Tagln','Cnn1','Meg3'],
    'Core mesenchymal (c10/c11/c12)':
                                  ['Col1a1','Col1a2','Col3a1','Sparc','Serpinh1','S100a4',
                                   'Acta2','Fap','Pdgfra','Postn','Mmp2','Loxl2','Timp1','Dcn','Sox2'],
}

print('Loading h5ad ...')
a = sc.read_h5ad(f'{M}/data/npp53_merge/npp53_edited.h5ad')
print(f'  n_obs={a.n_obs}, n_vars={a.n_vars}')

# Exclude Cycling + Ionocyte (per §10 policy)
a = a[a.obs['putative_state'].isin(FINE_ORDER)].copy()
print(f'  after excluding Cycling+Ionocyte: n_obs={a.n_obs}')

# Assign k4 major cell type
a.obs['k4_major'] = a.obs['putative_state'].map(FINE_TO_K4).astype('category')
a.obs['k4_major'] = a.obs['k4_major'].cat.set_categories(K4_ORDER, ordered=True)
a.obs['putative_state'] = a.obs['putative_state'].astype('category').cat.set_categories(FINE_ORDER, ordered=True)

# Log-normalize if needed (X should already be log1p from CellRanger + scanpy workflow — check)
# Use raw layer if available for markers
if 'counts' in a.layers:
    a_marker = sc.AnnData(X=a.layers['counts'].copy(), obs=a.obs.copy(), var=a.var.copy())
    sc.pp.normalize_total(a_marker, target_sum=1e4)
    sc.pp.log1p(a_marker)
else:
    a_marker = a.copy()

print('\nRunning rank_genes_groups per k4_major (Wilcoxon) ...')
sc.tl.rank_genes_groups(a_marker, groupby='k4_major', method='wilcoxon', n_genes=200)
# Extract top markers per group (excluding transgenes/mito)
TOP_N = 5
EXCLUDE_PREF = ('mt-','Rp[sl]','Malat1','Neat1','Gm[0-9]','Xist','Cas9','mCherry','GFP','Cre','Tomato')
import re
excl = re.compile('|'.join(EXCLUDE_PREF))
top_by_grp = {}
for grp in K4_ORDER:
    names = a_marker.uns['rank_genes_groups']['names'][grp]
    picked = []
    seen = set()
    for g in names:
        if excl.match(g): continue
        if g in seen: continue
        picked.append(g); seen.add(g)
        if len(picked) >= TOP_N: break
    top_by_grp[grp] = picked
    print(f'  {grp}: {picked}')

# Assemble Panel-A marker set
panel_a_markers = []
panel_a_group   = []
for grp in K4_ORDER:
    for g in top_by_grp[grp]:
        if g in a_marker.var_names and g not in panel_a_markers:
            panel_a_markers.append(g); panel_a_group.append(grp)

# Compute mean-expression + fraction-expressed per (group, gene) — for dot plot
def dotplot_data(adata, genes, groupby, cats):
    X = adata[:, genes].X
    if hasattr(X, 'toarray'): X = X.toarray()
    df = pd.DataFrame(X, columns=genes, index=adata.obs.index)
    df[groupby] = adata.obs[groupby].values
    means = df.groupby(groupby, observed=False)[genes].mean().reindex(cats)
    fracs = df.groupby(groupby, observed=False)[genes].apply(lambda g: (g > 0).mean()).reindex(cats)
    # scale means per gene (z-score across groups)
    scaled = (means - means.mean(axis=0)) / (means.std(axis=0) + 1e-9)
    return means, scaled, fracs

means_A, scaled_A, fracs_A = dotplot_data(a_marker, panel_a_markers, 'k4_major', K4_ORDER)
means_B, scaled_B, fracs_B = dotplot_data(a_marker, panel_a_markers, 'putative_state', FINE_ORDER)
# Apply c[id]_ prefix to Panel B y-axis labels
scaled_B.index = [CID_PREFIX[s] for s in scaled_B.index]
fracs_B.index  = [CID_PREFIX[s] for s in fracs_B.index]

# Genes to drop from canonical panel (Zejian: explicit + poorly-expressed markers)
DROP_GENES = {'Ar','Nkx3-1','Pbsn','Tgm4','Cnn1','Meg3','Postn','Krt19','Krt15','Msmb','Sbp','Snai1','Klk4'}
# Canonical markers — keep only those in var_names (avoid inventing genes)
canonical_flat, canonical_group, missing = [], [], []
for cat, gs in CANONICAL.items():
    for g in gs:
        if g in DROP_GENES: continue                        # explicitly excluded by user
        if g in a_marker.var_names and g not in canonical_flat:
            canonical_flat.append(g); canonical_group.append(cat)
        elif g not in a_marker.var_names:
            missing.append(g)
if missing: print(f'[canonical missing from var_names, dropped]: {missing}')
print(f'[user-dropped from canonical panel]: {sorted(DROP_GENES)}')
# Panel C: canonical × 12 fine states (was previously × 4 major types — Zejian wants finer resolution)
means_C, scaled_C, fracs_C = dotplot_data(a_marker, canonical_flat, 'putative_state', FINE_ORDER)
scaled_C.index = [CID_PREFIX[s] for s in scaled_C.index]
fracs_C.index  = [CID_PREFIX[s] for s in fracs_C.index]

# --- FIGURE ---
K4_COLORS = {'DistalAR_Lum':'#1f5fa8', 'ProxBasalSq':'#00a19b',
             'EMT_Bridge':'#e69f00',   'CoreMes':'#d13845'}
CAT_COLORS = {'Distal AR-lum (c1/c3/c4)':'#1f5fa8',
              'Prox/Basal/Squamous (c2/c5/c6/c7)':'#00a19b',
              'EMT bridge (c8/c9)':'#e69f00',
              'Core mesenchymal (c10/c11/c12)':'#d13845'}
# k4-major color for each fine state (for left color bar on Panels B & C)
STATE_TO_K4COLOR = {CID_PREFIX[s]: K4_COLORS[FINE_TO_K4[s]] for s in FINE_ORDER}

def draw_dotplot(ax, scaled, fracs, group_of_gene, group_colors, x_labels, y_labels,
                 max_size=190, left_bar_colors=None, is_fine_state=False):
    """left_bar_colors: dict row_label -> hex color (for k4 assignment strip); if None, skip.
    When present, the color strip sits to the LEFT of the y-tick text (Zejian's layout)."""
    n_g = len(scaled.columns); n_r = len(scaled.index)
    for i, row in enumerate(scaled.index):
        for j, gene in enumerate(scaled.columns):
            sz = fracs.loc[row, gene] * max_size + 3
            val = np.clip(scaled.loc[row, gene], -2, 2)
            ax.scatter(j, i, s=sz, c=[val], cmap='RdBu_r', vmin=-2, vmax=2, edgecolors='#333', lw=0.3)
    ax.set_xticks(range(n_g))
    ax.set_xticklabels(scaled.columns, rotation=90, fontsize=8)
    # colored bar under x-tick — one per gene (gene-block colors)
    for j, gene in enumerate(scaled.columns):
        c = group_colors.get(group_of_gene[j], '#333')
        ax.plot([j-0.3, j+0.3], [n_r+0.2, n_r+0.2], color=c, lw=3, solid_capstyle='butt')
    ax.set_yticks(range(n_r))
    ax.set_yticklabels(scaled.index if y_labels else [], fontsize=9)
    ax.set_ylim(-0.7, n_r + 0.7)
    ax.set_xlim(-0.7, n_g - 0.3)
    ax.invert_yaxis()
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis='both', color='#e8e8e8', lw=0.5, zorder=0)
    # NOTE: left color bar is now drawn in a dedicated side axis via draw_left_strip()
    # (so it sits to the LEFT of the y-tick text, not between text and plot).

def draw_left_strip(strip_ax, row_labels, colors_map):
    """[text (right-aligned) | color-rect] strip left of a dot-plot panel.
    Layout ends flush against plot's y-axis, so color rect sits between text and plot."""
    n = len(row_labels)
    # Color rect on the LEFT (annotation bar), text on the RIGHT — layout order [color | text | plot]
    for i, r in enumerate(row_labels):
        strip_ax.add_patch(plt.Rectangle((0.02, i-0.42), 0.14, 0.84,
                                          color=colors_map.get(r, '#333'), ec='#333', lw=0.5))
        strip_ax.text(0.20, i, r, ha='left', va='center', fontsize=9, color='#111',
                       clip_on=False)
    strip_ax.set_xlim(0, 1); strip_ax.set_ylim(-0.7, n + 0.7)
    strip_ax.invert_yaxis()
    strip_ax.axis('off')

# Compute figure size — 2 rows now (Panel B removed per Zejian)
n_A = len(panel_a_markers); n_C = len(canonical_flat)
fig_w = max(15, max(n_A, n_C)*0.42 + 4.5)
fig_h = 9   # reduced from 13 since Panel B gone
fig = plt.figure(figsize=(fig_w, fig_h))

# 2-row × 3-col: [color_strip (left) | main plot | legend (right)]
from matplotlib.gridspec import GridSpec
gs = GridSpec(2, 3, figure=fig, hspace=0.75, wspace=0.03,
              width_ratios=[5.0, max(n_A, n_C), 3.2],
              height_ratios=[1.1, 2.4])

axA       = fig.add_subplot(gs[0, 1])
legA_ax   = fig.add_subplot(gs[0, 2]); legA_ax.axis('off')
axC_strip = fig.add_subplot(gs[1, 0])
axC       = fig.add_subplot(gs[1, 1])

# Panel A — data-driven top-N × 4 major types
draw_dotplot(axA, scaled_A, fracs_A, panel_a_group, K4_COLORS,
             x_labels=True, y_labels=True)
# Panel B (formerly bottom) — canonical markers × 12 fine states with c[id]_ prefix + left k4 strip
draw_dotplot(axC, scaled_C, fracs_C, canonical_group, CAT_COLORS,
             x_labels=True, y_labels=False)
draw_left_strip(axC_strip, list(scaled_C.index), STATE_TO_K4COLOR)

# ---- Legends (in dedicated right column, NO overlap with plot area) ----
# Colorbar — anchored to Panel B/C
cax = fig.add_axes([0.905, 0.28, 0.010, 0.32])
sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=plt.Normalize(vmin=-2, vmax=2))
cb = fig.colorbar(sm, cax=cax); cb.set_label('z-scored mean expression', fontsize=9)
cb.ax.tick_params(labelsize=8)
# Size legend — right of colorbar
size_ax = fig.add_axes([0.945, 0.28, 0.045, 0.32]); size_ax.axis('off')
size_ax.text(0.5, 1.05, '% expressing', fontsize=9, ha='center', fontweight='bold')
for i, f in enumerate([0.25, 0.5, 0.75, 1.0]):
    size_ax.scatter(0.35, 0.92 - i*0.24, s=f*190 + 3, c='#666', edgecolor='#333')
    size_ax.text(0.65, 0.92 - i*0.24, f'{int(f*100)}%', va='center', fontsize=8)
size_ax.set_xlim(0, 1); size_ax.set_ylim(0, 1.1)
# k4-major color legend — inside Panel-A legend column
legA_ax.text(0.5, 1.02, 'k=4 major cell type', fontsize=10, ha='center',
             fontweight='bold', color='#0f2747')
for i, k in enumerate(K4_ORDER):
    legA_ax.add_patch(plt.Rectangle((0.05, 0.85 - i*0.15), 0.14, 0.10,
                                      color=K4_COLORS[k], ec='#333'))
    legA_ax.text(0.24, 0.90 - i*0.15, k, fontsize=9, va='center')
legA_ax.set_xlim(0,1); legA_ax.set_ylim(0, 1.1)

# NO suptitle (per user)

for fmt in ['pdf','png']:
    fig.savefig(OUT/f'fig_6_7_k4_major_markers_v2.{fmt}', dpi=180, bbox_inches='tight')
plt.close(fig)
print(f'\nsaved: {OUT}/fig_6_7_k4_major_markers_v2.png (NEW filename — original v1 preserved)')

# Also save the top-markers table
tab = pd.DataFrame({'k4_major': panel_a_group, 'marker_gene': panel_a_markers})
tab.to_csv(OUT/'fig_6_7_top_markers_per_k4.csv', index=False)
print(f'saved: {OUT}/fig_6_7_top_markers_per_k4.csv')

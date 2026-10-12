#!/usr/bin/env python
"""Dotplot of leading-edge genes from the 2 most promising canonical Hallmark
pathways enriched in Fate 2 (c3+c4+c5) via clusterProfiler Wilcoxon-DE ranking:
  - HALLMARK_INTERFERON_ALPHA_RESPONSE
  - HALLMARK_CHOLESTEROL_HOMEOSTASIS
(Peroxisome dropped — shares Fdps/Idi1/Ech1 with Cholesterol and is redundant.)
Rows = 12 CARTA states in fate-map trajectory order (Fate 2 target rows highlighted).
c3 relabeled from Luminal-ARdep → Luminal-ARdep for display only (h5ad key unchanged).
"""
import numpy as np, pandas as pd
import anndata as ad
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BASE = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
H5AD = BASE/'data/npp53_merge/npp53_edited.h5ad'
GSEA = BASE/'analysis/evolution_features_v2/carta_fine12/biology_review/gsea_fate2_Hallmark_result.tsv'
OUT  = BASE/'analysis/evolution_features_v2/carta_fine12/biology_review'

STATE_ORDER_FATE = [
    'Luminal-Distal','Luminal-Proximal','Luminal-Squamous','Basal-like',
    'Luminal-Persistent','Luminal-ARdep','Luminal-CR',
    'Partial-EMT','Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF',
]
SHORT = {s: f'c{i+1}' for i, s in enumerate([
    'Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR','Luminal-Persistent',
    'Luminal-Squamous','Basal-like','Partial-EMT','Mes-Quiescent','Mes-Collagen',
    'Mes-Stemness','Mes-myCAF'])}
STATE_COLORS = {
    'Luminal-Distal':'#1f5fa8','Luminal-Proximal':'#9467bd','Luminal-ARdep':'#d23030',
    'Luminal-CR':'#f1903d','Luminal-Persistent':'#117c70','Luminal-Squamous':'#f78ec7',
    'Basal-like':'#8c564b','Partial-EMT':'#bcbd22','Mes-Quiescent':'#c49c94',
    'Mes-Collagen':'#5b5b5b','Mes-Stemness':'#8c1a1a','Mes-myCAF':'#4daf4a',
}
FATE2 = {'Luminal-ARdep','Luminal-CR','Luminal-Persistent'}
# Display re-label (h5ad putative_state values unchanged; only y-tick text changes)
DISPLAY = {s: s for s in STATE_ORDER_FATE}
DISPLAY['Luminal-ARdep'] = 'Luminal-ARdep'

# ── Load leading-edge from GSEA result and preserve LE-ranked order ─────
gsea = pd.read_csv(GSEA, sep='\t').set_index('ID')
TARGET_PATHS = [
    ('IFNα Response',            'HALLMARK_INTERFERON_ALPHA_RESPONSE',   '#7a3ea5'),
    ('Cholesterol Biosynthesis', 'HALLMARK_CHOLESTEROL_HOMEOSTASIS',     '#c44e52'),
]
GENE_GROUPS = []
for label, path, color in TARGET_PATHS:
    le = gsea.loc[path, 'core_enrichment'].split('/')
    print(f'{label} ({path}): {len(le)} LE genes: {le}')
    GENE_GROUPS.append((label, le, color))

# ── Compute mean expr + frac expressing per (state × gene) ──────────────
print('loading h5ad ...')
adata = ad.read_h5ad(H5AD)
adata = adata[adata.obs['putative_state'].isin(STATE_ORDER_FATE)].copy()
print(f'  {adata.shape[0]} cells across 12 states')

all_genes = []
group_boundaries = []
group_labels = []
group_colors  = []
for gname, genes, gcolor in GENE_GROUPS:
    genes_present = [g for g in genes if g in adata.var_names]
    all_genes.extend(genes_present)
    group_boundaries.append(len(all_genes))
    group_labels.append(gname)
    group_colors.append(gcolor)

n_states = len(STATE_ORDER_FATE)
n_genes  = len(all_genes)
print(f'  n_genes total: {n_genes}')
mean_expr = np.zeros((n_states, n_genes))
frac_expr = np.zeros((n_states, n_genes))
for i, s in enumerate(STATE_ORDER_FATE):
    idx = np.where(adata.obs['putative_state'].values == s)[0]
    if len(idx) == 0: continue
    X = adata[idx, all_genes].X
    if hasattr(X,'toarray'): X = X.toarray()
    frac_expr[i, :] = (X > 0).mean(axis=0)
    mean_expr[i, :] = X.mean(axis=0)

mu = mean_expr.mean(axis=0, keepdims=True)
sd = mean_expr.std(axis=0, keepdims=True); sd[sd==0] = 1
z_expr = (mean_expr - mu) / sd

# ── Draw dotplot ────────────────────────────────────────────────────────
fig_w = 0.42 * n_genes + 6
fig_h = 0.55 * n_states + 3
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
cmap = LinearSegmentedColormap.from_list('bwr_soft', ['#3355aa','#f5f5f5','#c22f2f'], N=256)
max_frac = 0.9; max_dot = 420

# Highlight Fate 2 rows (c5, c3, c4) with a pink background band
for i, s in enumerate(STATE_ORDER_FATE):
    if s in FATE2:
        ax.axhspan(i-0.4, i+0.4, facecolor='#fce4ec', alpha=0.5, zorder=0)

for i in range(n_states):
    for j in range(n_genes):
        f = min(frac_expr[i, j], max_frac)
        if f < 0.01: continue
        color = cmap((np.clip(z_expr[i, j], -2, 2) + 2) / 4)
        size = 40 + (max_dot - 40) * (f / max_frac)
        ax.scatter(j, i, s=size, c=[color], edgecolor='#333', linewidth=0.4, zorder=3)

ax.set_yticks(range(n_states))
ax.set_yticklabels([f'{SHORT[s]} {DISPLAY[s]}' for s in STATE_ORDER_FATE], fontsize=10)
for ytick, s in zip(ax.get_yticklabels(), STATE_ORDER_FATE):
    ytick.set_color(STATE_COLORS[s]); ytick.set_fontweight('bold')
ax.set_xticks(range(n_genes))
ax.set_xticklabels(all_genes, rotation=60, ha='right', fontsize=9)
ax.set_xlim(-0.5, n_genes-0.5); ax.set_ylim(n_states-0.5, -0.5)
ax.grid(True, alpha=0.25, zorder=1)
ax.tick_params(top=False, right=False, left=False, bottom=False)
for spine in ('top','right','left','bottom'): ax.spines[spine].set_alpha(0.3)

# Group brackets + labels above the plot
prev = 0
for gname, gb, gcolor in zip(group_labels, group_boundaries, group_colors):
    start = prev; end = gb - 1
    mid = (start + end) / 2
    ax.plot([start-0.4, end+0.4], [-1.2, -1.2],
             color=gcolor, lw=2.8, solid_capstyle='round', clip_on=False)
    ax.text(mid, -1.7, gname, ha='center', va='bottom',
             fontsize=11, fontweight='bold', color=gcolor, clip_on=False)
    prev = gb

# Colorbar + size legend
cbar_ax = fig.add_axes([0.93, 0.35, 0.010, 0.30])
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-2, vmax=2))
cbar = fig.colorbar(sm, cax=cbar_ax); cbar.set_label('z-scored mean expr.\n(across 12 states)',
                                                       fontsize=8, rotation=270, labelpad=18)
cbar.ax.tick_params(labelsize=7)
sz_ax = fig.add_axes([0.93, 0.10, 0.05, 0.15]); sz_ax.axis('off')
for k, f in enumerate([0.10, 0.30, 0.60, 0.90]):
    y = 0.85 - k * 0.20
    sz_ax.scatter([0.15], [y], s=40 + (max_dot - 40) * (f / max_frac),
                   c='#999', edgecolor='#333', linewidth=0.4)
    sz_ax.text(0.42, y, f'{int(f*100)}%', va='center', fontsize=8, color='#333')
sz_ax.text(0.30, 1.02, 'frac expr.', ha='center', fontsize=8, color='#333', fontweight='bold')
sz_ax.set_xlim(0, 1); sz_ax.set_ylim(0, 1)

fig.suptitle('Fate 2 = ISG–Sterol biosynthesis program (c5 Persistent + c3 ARdep + c4 CR, pink band) — leading-edge genes from clusterProfiler GSEA',
             fontsize=12, fontweight='bold', color='#0f2747', y=1.005)
fig.text(0.5, -0.015,
         'Leading-edge = genes contributing to the running enrichment score up to peak (from clusterProfiler::GSEA).   '
         'Column order preserves LE ranking; dot color = per-column z-score; dot size = fraction of cells expressing.',
         ha='center', fontsize=9, color='#555', style='italic')

plt.subplots_adjust(top=0.88, bottom=0.22, left=0.10, right=0.90)
for fmt in ['png','pdf']:
    fig.savefig(OUT/f'npp53_fate2_LE_canonical_dotplot.{fmt}', dpi=200,
                 bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'saved: {OUT}/npp53_fate2_LE_canonical_dotplot.png')

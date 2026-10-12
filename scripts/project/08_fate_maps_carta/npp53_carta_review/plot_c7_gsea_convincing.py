#!/usr/bin/env python
"""c7 Basal-like — final GSEA bar chart with ONLY convincing pathways
(pathway-specific canonical LE genes present, not stress-generic / misnomer /
redundant). Single-panel: NES on x, color = -log10(FDR), y = pathway."""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

BASE = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
GSEA = BASE/'analysis/evolution_features_v2/carta_fine12/biology_review/gsea_c7_combined_result.tsv'
OUT  = BASE/'analysis/evolution_features_v2/carta_fine12/biology_review'

# CONVINCING pathways only (from LE canonical-marker review — see prior turn):
#   REMOVED: TNFA_NFKB (stress-generic AP-1 IEGs), WNT_BETA (misnomer = Notch),
#            KRAS_DN (basal-identity not KRAS), ESTROGEN_LATE (generic epithelial),
#            ADIPOGENESIS (mitochondrial-redundant with OxPhos),
#            MTORC1_SIGNALING (generic proteostasis).
CONVINCING = [
    'HALLMARK_TGF_BETA_SIGNALING',        # UP — Skil, Id2/3, Smad1, Pmepa1 canonical
    'HALLMARK_P53_PATHWAY',                # UP — Sfn, Cdkn1a, Btg2, Perp canonical
    'Bluemn_FGF_pathway',                  # UP — Fgfr2, Fgfr3, Dusp6 canonical FGF
    'HALLMARK_APICAL_JUNCTION',            # UP — Cdh1/3, laminins, integrins canonical
    'HALLMARK_PROTEIN_SECRETION',          # DOWN — Sec24d, Copb2, Sec22b, Arf1 canonical COP
    'HALLMARK_FATTY_ACID_METABOLISM',      # DOWN — 11/30 canonical β-ox enzymes
    'HALLMARK_OXIDATIVE_PHOSPHORYLATION',  # DOWN — 22/47 canonical ETC subunits
]

# Pretty labels (with pathway family tags for legibility)
PRETTY = {
    'HALLMARK_TGF_BETA_SIGNALING':        'TGF-β signaling (Hallmark)',
    'HALLMARK_P53_PATHWAY':                'p53 stress response (Hallmark)',
    'Bluemn_FGF_pathway':                  'FGF signaling (Bluemn — CRPC)',
    'HALLMARK_APICAL_JUNCTION':            'Apical junction / epithelial (Hallmark)',
    'HALLMARK_PROTEIN_SECRETION':          'COP protein secretion (Hallmark)',
    'HALLMARK_FATTY_ACID_METABOLISM':      'Fatty-acid β-oxidation (Hallmark)',
    'HALLMARK_OXIDATIVE_PHOSPHORYLATION':  'Oxidative phosphorylation (Hallmark)',
}

g = pd.read_csv(GSEA, sep='\t')
g = g[g['ID'].isin(CONVINCING)].copy()
g['neg_log10_padj'] = -np.log10(np.maximum(g['p.adjust'], 1e-300))
g['label'] = g['ID'].map(PRETTY)
g = g.sort_values('NES', ascending=False)   # top = strongest UP; bottom = strongest DOWN
print(g[['ID','NES','p.adjust','neg_log10_padj','setSize']].to_string(index=False))

# ── Draw ────────────────────────────────────────────────────────────────
vmax = g['neg_log10_padj'].max() * 1.05
norm = Normalize(vmin=0, vmax=vmax)
cmap = cm.get_cmap('viridis')

fig, ax = plt.subplots(figsize=(12.5, 6.5))
y = np.arange(len(g))
colors = [cmap(norm(v)) for v in g['neg_log10_padj']]
bars = ax.barh(y, g['NES'].values, color=colors,
                edgecolor='#333', linewidth=0.5, height=0.72)
ax.set_yticks(y)
ax.set_yticklabels(g['label'].values, fontsize=11)
ax.invert_yaxis()   # top of chart = strongest UP
ax.axvline(0, color='#333', lw=1.0)
ax.set_xlabel('NES (Normalized Enrichment Score)', fontsize=12, fontweight='bold')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.grid(True, axis='x', alpha=0.3)

# NES + padj + n annotation: always sit on the OUTBOARD side of the bar tail,
# i.e. RIGHT of a positive-NES bar tip AND RIGHT of a negative-NES bar's
# near-zero end (i.e. INSIDE the negative bar, next to x=0 on the right side
# of the bar tail). This keeps ALL annotations to the RIGHT of x=0, never
# spilling into the y-axis pathway-label region on the left.
for bar, nes, padj, size in zip(bars, g['NES'], g['p.adjust'], g['setSize']):
    if nes >= 0:
        # positive bar: annotation to the RIGHT of the bar tip
        ax.text(nes + 0.06, bar.get_y() + bar.get_height()/2,
                f'NES={nes:+.2f}   padj={padj:.1e}   n={int(size)}',
                va='center', ha='left', fontsize=9.5, color='#111', fontweight='bold')
    else:
        # negative bar: annotation to the RIGHT of x=0 (i.e. inside the bar
        # graphics on the near-zero side), ha=left so it reads leftward
        ax.text(0.06, bar.get_y() + bar.get_height()/2,
                f'NES={nes:+.2f}   padj={padj:.1e}   n={int(size)}',
                va='center', ha='left', fontsize=9.5, color='#111', fontweight='bold')

# Wider symmetric x-limit so annotations always fit outside the bar tail
xmax = max(abs(g['NES'].max()), abs(g['NES'].min())) * 1.95
ax.set_xlim(-xmax, xmax)

# UP / DOWN axis banners
ax.text(xmax*0.98, -0.55, 'ENRICHED IN c7 (UP)', ha='right', va='center',
         fontsize=9, color='#c22f2f', fontweight='bold', alpha=0.7)
ax.text(-xmax*0.98, -0.55, 'DEPLETED IN c7 (DOWN)', ha='left', va='center',
         fontsize=9, color='#3355aa', fontweight='bold', alpha=0.7)

# Colorbar OUTSIDE the plot on the right — bumped further out from panel
cbar_ax = fig.add_axes([0.98, 0.20, 0.014, 0.60])
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('−log$_{10}$(FDR)', fontsize=11, rotation=270, labelpad=22, fontweight='bold')
cbar.ax.tick_params(labelsize=9)
# Significance guide-lines (labels moved BELOW the line to avoid overlapping ticks)
for thr, label, offset_y in [(-np.log10(0.05), 'FDR=0.05', -0.35),
                              (-np.log10(0.01), 'FDR=0.01', -0.35)]:
    if thr <= vmax:
        cbar.ax.axhline(thr, color='#e74c3c', lw=1)
        cbar.ax.text(3.5, thr + offset_y, label, va='center', fontsize=7.5,
                      color='#e74c3c', fontweight='bold')

fig.suptitle('c7 Basal-like — Convincing GSEA hits (canonical LE-gene-supported, 4 UP + 3 DOWN)',
             fontsize=13.5, fontweight='bold', color='#0f2747', y=1.00)
fig.text(0.5, -0.02,
         'clusterProfiler GSEA on Wilcoxon-DE-ranked c7 vs 11 other CARTA states, 75-pathway library (Hallmark + PCa/CRPC).  '
         'Only pathways with pathway-specific canonical LE genes shown (see LE-canonical review); stress-generic AP-1 / redundant / misnomer paths dropped.',
         ha='center', fontsize=9.5, color='#555', style='italic')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.32, right=0.94)
for fmt in ['png','pdf']:
    fig.savefig(OUT/f'npp53_c7_gsea_convincing.{fmt}', dpi=200,
                 bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'saved: {OUT}/npp53_c7_gsea_convincing.png')

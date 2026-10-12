#!/usr/bin/env python
"""Fig 6.4a + 6.4b · 4-row x 3-col UMAP/tSNE/Diffmap grids.
   Row 1: cell state · Row 2: Intact vs Cas · Row 3: 3-level treatment · Row 4: Intact vs Cas+Enz
   FIX (2026-09-02): rows 2 + 4 now FILTER OUT cells outside the pairwise arms (no gray NA cells).
   FIX (2026-09-02b): dedicated legend-slot column (col 1) so legend text no longer clips into tSNE;
                      bigger legend/title fonts; row labels added on far left;
                      no em-dashes (mid-dot separator throughout)."""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec

M = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
OUT = M/'1-mouse_based/Npp53_tracing/figures/npp53_plots'

# c[N]_ display map (14 -> also handles Ionocyte + Cycling for Fig 6.4a)
STATE_MAP_14 = {
    'Luminal-Distal':'c1_Lum-Distal','Luminal-Proximal':'c2_Lum-Proximal','Luminal-ARdep':'c3_Lum-ARdep',
    'Luminal-CR':'c4_Lum-CR','Luminal-Persistent':'c5_Lum-Persistent','Luminal-Squamous':'c6_Lum-Squamous',
    'Basal-like':'c7_Basal-like','Partial-EMT':'c8_Partial-EMT','Mes-Quiescent':'c9_Mes-Quiescent',
    'Mes-Collagen':'c10_Mes-Collagen','Mes-Stemness':'c11_Mes-Stemness','Mes-myCAF':'c12_Mes-myCAF',
    'Ionocyte':'c13_Ionocyte','Cycling':'c14_Cycling',
}
STATE_COLORS_14 = {
    'c1_Lum-Distal':'#26a69a','c2_Lum-Proximal':'#ec407a','c3_Lum-ARdep':'#7cb342','c4_Lum-CR':'#ffb300',
    'c5_Lum-Persistent':'#5c6bc0','c6_Lum-Squamous':'#8d6e63','c7_Basal-like':'#4a7fbf','c8_Partial-EMT':'#ffca28',
    'c9_Mes-Quiescent':'#bcaaa4','c10_Mes-Collagen':'#78909c','c11_Mes-Stemness':'#ce93d8','c12_Mes-myCAF':'#a1887f',
    'c13_Ionocyte':'#c44e52','c14_Cycling':'#9e9e9e',
}
TREAT_COLORS_2 = {'Intact':'#4a7fbf','Castration':'#e0955f'}
TREAT_COLORS_3 = {'Intact':'#4a7fbf','Castration':'#e0955f','Castration + Enzalutamide':'#c44e52'}
TREAT_COLORS_IE = {'Intact':'#4a7fbf','Cas+Enz':'#c44e52'}

def build_grid(cell_set_tag, umap_csv, tsne_csv, dm_csv, obs_csv, outfile, states_present):
    print(f'\n=== {cell_set_tag} -> {outfile} ===')
    u = pd.read_csv(M/f'1-mouse_based/Npp53_tracing/data/npp53_merge/{umap_csv}', index_col=0)
    t = pd.read_csv(M/f'1-mouse_based/Npp53_tracing/data/npp53_merge/{tsne_csv}', index_col=0)
    d = pd.read_csv(M/f'1-mouse_based/Npp53_tracing/data/npp53_merge/{dm_csv}',   index_col=0)
    obs = pd.read_csv(M/f'1-mouse_based/Npp53_tracing/data/npp53_merge/{obs_csv}',index_col=0)
    idx = u.index
    state_display = obs.loc[idx, 'putative_state'].map(STATE_MAP_14)
    treat = obs.loc[idx, 'treatment']
    n = len(idx); print(f'  cells: {n}, NA state: {state_display.isna().sum()}')
    U = u.values; T = t.values; D = d.values

    row_labels = [
        f'Cell state · {len(states_present)}',
        'Intact vs Castration',
        '3-level treatment',
        'Intact vs Cas+Enz',
    ]

    # Layout: 4 rows x 4 cols.
    #   col 0 = UMAP, col 1 = legend slot (no axes), col 2 = tSNE, col 3 = Diffmap (widened).
    fig = plt.figure(figsize=(26, 28))
    gs = gridspec.GridSpec(
        4, 4, figure=fig,
        width_ratios=[1, 0.75, 1, 1.4],
        hspace=0.30, wspace=0.16,
        left=0.055, right=0.985, top=0.955, bottom=0.02,
    )

    def panel(ax, XY, mask, colors, title, legend_title=None, legend_ncol=1, dot_size=1.6):
        # plot only cells in mask, colored by category (colors_key set outside)
        for cat, col in colors.items():
            m2 = mask & (colors_key == cat)
            if m2.sum() == 0: continue
            ax.scatter(XY[m2.values, 0], XY[m2.values, 1],
                       c=col, s=dot_size, alpha=0.85, edgecolors='none', label=cat)
        ax.set_title(title, fontsize=17, fontweight='bold', color='#0f2747', pad=8)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_edgecolor('#ccc'); sp.set_linewidth(0.5)
        if legend_title:
            ax.legend(
                fontsize=13.5, title=legend_title, title_fontsize=15, frameon=False,
                loc='center left', bbox_to_anchor=(1.04, 0.5),
                markerscale=4.0, borderpad=0.6, labelspacing=0.55,
                handletextpad=0.5, ncol=legend_ncol,
            )
        ax.set_xlabel(''); ax.set_ylabel('')

    row_col0_axes = []   # for placing row labels

    def build_row(row_i, mask, colors, key_series, titles, legend_title):
        nonlocal colors_key
        colors_key = key_series
        axes = []
        for j, gs_col in enumerate([0, 2, 3]):
            ax = fig.add_subplot(gs[row_i, gs_col])
            axes.append(ax)
        # Per-panel dot sizes: UMAP + tSNE = 1.6 (dense clouds); Diffmap = 7 (sparse fan, needs bigger dots)
        DOT_SIZES = [1.6, 1.6, 7.0]
        for j, (ax, XY, title) in enumerate(zip(axes, [U, T, D], titles)):
            panel(ax, XY, mask, colors, title,
                  legend_title=legend_title if j == 0 else None,
                  dot_size=DOT_SIZES[j])
        row_col0_axes.append(axes[0])

    colors_key = None   # placeholder overwritten inside build_row

    # === Row 1: cell state (all cells, colored by state) ===
    build_row(
        0, state_display.notna(),
        {s: STATE_COLORS_14[s] for s in states_present},
        state_display,
        ['UMAP · cell state', 'tSNE · cell state', 'Diffmap (DC2/DC3) · cell state'],
        legend_title='Cell state',
    )

    # === Row 2: Intact vs Castration (filter out Cas+Enz) ===
    build_row(
        1, treat.isin(['Intact','Castration']),
        TREAT_COLORS_2, treat,
        ['UMAP · Intact vs Castration',
         'tSNE · Intact vs Castration',
         'Diffmap (DC2/DC3) · Intact vs Castration'],
        legend_title='Treatment',
    )

    # === Row 3: 3-level treatment (all cells) ===
    build_row(
        2, treat.notna(),
        TREAT_COLORS_3, treat,
        ['UMAP · 3-level treatment',
         'tSNE · 3-level treatment',
         'Diffmap (DC2/DC3) · 3-level treatment'],
        legend_title='Treatment',
    )

    # === Row 4: Intact vs Cas+Enz (filter out Castration-only) ===
    treat_ie = treat.replace({'Castration + Enzalutamide': 'Cas+Enz'})
    build_row(
        3, treat.isin(['Intact','Castration + Enzalutamide']),
        TREAT_COLORS_IE, treat_ie,
        ['UMAP · Intact vs Cas+Enz',
         'tSNE · Intact vs Cas+Enz',
         'Diffmap (DC2/DC3) · Intact vs Cas+Enz'],
        legend_title='Treatment',
    )

    # Row labels on far-left (vertical). Needs positions -> draw first.
    fig.canvas.draw()
    for i, ax0 in enumerate(row_col0_axes):
        pos = ax0.get_position()
        y_center = (pos.y0 + pos.y1) / 2
        fig.text(
            0.010, y_center, row_labels[i],
            fontsize=16, fontweight='bold', color='#0f2747',
            ha='left', va='center', rotation=90,
        )

    fig.suptitle(cell_set_tag, fontsize=19, fontweight='bold', color='#0f2747', y=0.985)
    for fmt in ['pdf','png']:
        fig.savefig(OUT/outfile.replace('.png', f'.{fmt}'), dpi=160, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved: {outfile}')

# Fig 6.4a · no-Cycling (13 states, including Ionocyte)
build_grid(
    'Fig 6.4a · Cycling-removed view (13 states, 18,541 cells)',
    'umap_no_cycling.csv','tsne_no_cycling.csv','diffmap_no_cycling.csv','obs_no_cycling.csv',
    'sec6_4a_no_cycling_grid.png',
    states_present=[STATE_MAP_14[s] for s in ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
                                                'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
                                                'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF','Ionocyte']]
)

# Fig 6.4b · no-Cycling + no-Ionocyte (12 states)
build_grid(
    'Fig 6.4b · Final tumor cell state · No-Cycling + No-Ionocyte (12 states, 18,241 cells)',
    'umap_no_cycling_no_ionocyte.csv','tsne_no_cycling_no_ionocyte.csv','diffmap_no_cycling_no_ionocyte.csv','obs_no_cycling_no_ionocyte.csv',
    'sec6_4b_no_cycling_no_ionocyte_grid.png',
    states_present=[STATE_MAP_14[s] for s in ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
                                                'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
                                                'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']]
)
print('\ndone')

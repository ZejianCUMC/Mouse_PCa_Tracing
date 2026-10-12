#!/usr/bin/env python
"""Yang 2022 KP-Tracer style Meta-vote plots for FitnessSignature construction.

Panels:
  A. Consensus heatmap - top FSU + top FSD genes × 25 tumors, colored by per-tumor log2FC
  B. Distribution of n_tumors DE per gene (histogram of freq_total by direction)
  C. Waterfall of top-40 FSU + top-40 FSD genes (n tumors "up" or "dn", horizontal bars)
  D. Meta-consensus volcano - mean_log2FC vs consistency; color = direction
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
mpl.rcParams['pdf.fonttype']=42; mpl.rcParams['ps.fonttype']=42
mpl.rcParams['font.family']='sans-serif'

FS = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/FitnessSignature')
OUT = FS/'meta_vote_yang2022_plots'
OUT.mkdir(parents=True, exist_ok=True)

# ---- Load data ----
meta = pd.read_csv(FS/'FitnessSignature_meta.csv')
pergene = pd.read_csv(FS/'FitnessSignature_per_gene_per_sample.csv.gz')
print(f'meta rows: {len(meta):,}  ·  per_gene rows: {len(pergene):,}  ·  n_tumors: {pergene["sample"].nunique()}')

# ---- Load FSU / FSD gene lists (final meta-vote picks: freq_maj >= 2 & ratio > 0.5) ----
FSU = [g.strip() for g in (FS/'FitnessPosSignature_All_mouse.txt').read_text().splitlines() if g.strip()]
FSD = [g.strip() for g in (FS/'FitnessNegSignature_All_mouse.txt').read_text().splitlines() if g.strip()]
print(f'FSU (mouse): {len(FSU):,}  ·  FSD (mouse): {len(FSD):,}')

# ---- Panel A: consensus heatmap of top FSU + FSD genes × 25 tumors ----
N_SHOW = 40  # top N of each direction
FDR_CUT, LFC_CUT = 0.10, np.log2(1.25)

meta_up = meta[meta.majority_dir=='up'].sort_values(['freq_maj','mean_log2FC'], ascending=[False,False])
meta_dn = meta[meta.majority_dir=='dn'].sort_values(['freq_maj','mean_log2FC'], ascending=[False,True])
# restrict to those that actually passed final filter (in FSU/FSD)
meta_up_in = meta_up[meta_up.gene.isin(FSU)].head(N_SHOW)
meta_dn_in = meta_dn[meta_dn.gene.isin(FSD)].head(N_SHOW)
show_genes = meta_up_in.gene.tolist() + meta_dn_in.gene.tolist()

# Build gene × sample log2FC matrix
tumors = sorted(pergene['sample'].unique())
mat = pergene[pergene.gene.isin(show_genes)].pivot_table(
        index='gene', columns='sample', values='log2FC', aggfunc='first').reindex(show_genes)[tumors]
# significance mask
sig = pergene[pergene.gene.isin(show_genes)].pivot_table(
        index='gene', columns='sample', values='padj', aggfunc='first').reindex(show_genes)[tumors]

fig, ax = plt.subplots(figsize=(12, max(6, len(show_genes)*0.14)), constrained_layout=True)
lim = 3
im = ax.imshow(mat.values, aspect='auto', cmap='RdBu_r', vmin=-lim, vmax=lim, interpolation='nearest')
# overlay significant cells with black dot
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        p = sig.iat[i,j] if not pd.isna(sig.iat[i,j]) else 1
        lfc = mat.iat[i,j] if not pd.isna(mat.iat[i,j]) else 0
        if p < FDR_CUT and abs(lfc) > LFC_CUT:
            ax.plot(j, i, marker='.', color='black', markersize=2)
ax.set_yticks(range(len(show_genes)))
ax.set_yticklabels(show_genes, fontsize=6)
ax.set_xticks(range(len(tumors)))
ax.set_xticklabels(tumors, rotation=90, fontsize=7)
# horizontal separator between FSU and FSD blocks
ax.axhline(y=N_SHOW-0.5, color='black', linewidth=1)
ax.text(-1.5, N_SHOW/2 - 0.5, 'FSU\n(fitness-up)', ha='right', va='center', fontsize=9,
        fontweight='bold', color='#c44e52', rotation=90)
ax.text(-1.5, N_SHOW + N_SHOW/2 - 0.5, 'FSD\n(fitness-dn)', ha='right', va='center', fontsize=9,
        fontweight='bold', color='#4a7fbf', rotation=90)
plt.colorbar(im, ax=ax, fraction=0.02, label='per-tumor log2FC (top vs bottom fitness decile)')
ax.set_title(f'A · Meta-consensus heatmap - top-{N_SHOW} FSU + top-{N_SHOW} FSD genes across 25 tumors\n(black dot = tumor passes FDR<{FDR_CUT} & |LFC|>{LFC_CUT:.2f})',
             fontsize=11)
fig.savefig(OUT/'A_metavote_heatmap.pdf', dpi=180, bbox_inches='tight')
fig.savefig(OUT/'A_metavote_heatmap.png', dpi=140, bbox_inches='tight')
plt.close(fig)
print(f'  saved A_metavote_heatmap.png')

# ---- Panel B: distribution of freq_total (n_tumors DE) per direction ----
fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
bins = np.arange(0, meta.freq_total.max()+2) - 0.5
ax.hist(meta[meta.majority_dir=='up'].freq_total.values, bins=bins,
        color='#c44e52', alpha=0.7, label=f'up-direction (n={(meta.majority_dir=="up").sum():,})', edgecolor='#7a2b2c')
ax.hist(meta[meta.majority_dir=='dn'].freq_total.values, bins=bins,
        color='#4a7fbf', alpha=0.7, label=f'dn-direction (n={(meta.majority_dir=="dn").sum():,})', edgecolor='#2a4870')
ax.axvline(x=2, color='black', linestyle='--', linewidth=1, label='meta-vote threshold (>=2 tumors)')
ax.set_xlabel('n_tumors gene is DE (FDR<0.10 & |LFC|>log2(1.25))')
ax.set_ylabel('n_genes')
ax.set_title('B · Distribution of per-gene DE count across 25 tumors')
ax.legend(fontsize=9, loc='upper right')
fig.savefig(OUT/'B_metavote_ntumor_hist.pdf', dpi=180, bbox_inches='tight')
fig.savefig(OUT/'B_metavote_ntumor_hist.png', dpi=140, bbox_inches='tight')
plt.close(fig)
print(f'  saved B_metavote_ntumor_hist.png')

# ---- Panel C: Waterfall of top-N FSU / FSD genes by DE-count ----
N_WATER = 30
top_up = meta[meta.gene.isin(FSU)].sort_values(['freq_maj','mean_log2FC'], ascending=[True, True]).tail(N_WATER)
top_dn = meta[meta.gene.isin(FSD)].sort_values(['freq_maj','mean_log2FC'], ascending=[True, False]).tail(N_WATER)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, max(6, N_WATER*0.24)), constrained_layout=True)
axL.barh(range(N_WATER), top_up.freq_maj.values, color='#c44e52')
axL.set_yticks(range(N_WATER)); axL.set_yticklabels(top_up.gene.values, fontsize=8)
axL.set_xlabel(f'n_tumors DE up  (threshold FDR<{FDR_CUT}, |LFC|>{LFC_CUT:.2f})')
axL.set_title(f'C-left · Top {N_WATER} FSU genes by meta-vote count')
axL.axvline(x=2, color='black', linestyle='--', linewidth=1)
axL.set_xlim(0, meta.freq_total.max()+1)

axR.barh(range(N_WATER), top_dn.freq_maj.values, color='#4a7fbf')
axR.set_yticks(range(N_WATER)); axR.set_yticklabels(top_dn.gene.values, fontsize=8)
axR.set_xlabel(f'n_tumors DE dn  (threshold FDR<{FDR_CUT}, |LFC|>{LFC_CUT:.2f})')
axR.set_title(f'C-right · Top {N_WATER} FSD genes by meta-vote count')
axR.axvline(x=2, color='black', linestyle='--', linewidth=1)
axR.set_xlim(0, meta.freq_total.max()+1)
fig.savefig(OUT/'C_metavote_waterfall.pdf', dpi=180, bbox_inches='tight')
fig.savefig(OUT/'C_metavote_waterfall.png', dpi=140, bbox_inches='tight')
plt.close(fig)
print(f'  saved C_metavote_waterfall.png')

# ---- Panel D: Meta-consensus scatter (volcano-like) ----
fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
m = meta.copy()
m['in_FSU'] = m.gene.isin(FSU)
m['in_FSD'] = m.gene.isin(FSD)
m['color'] = '#dcdcdc'
m.loc[m.in_FSU, 'color'] = '#c44e52'
m.loc[m.in_FSD, 'color'] = '#4a7fbf'
m['size'] = 3 + m.freq_total * 1.5
# jitter freq_maj slightly for readability
np.random.seed(0)
xjit = m.freq_maj + np.random.uniform(-0.15, 0.15, size=len(m))
ax.scatter(xjit, m.mean_log2FC, c=m.color, s=m['size'], alpha=0.55, edgecolors='none', rasterized=True)
ax.axhline(0, color='#666', linewidth=0.6)
ax.axvline(2, color='black', linestyle='--', linewidth=0.8, alpha=0.6, label='meta-vote threshold ≥2')
ax.set_xlabel('n_tumors DE in majority direction (freq_maj)')
ax.set_ylabel('mean log2FC across tumors (meta-consensus)')
ax.set_title(f'D · Meta-consensus scatter  ·  FSU={m.in_FSU.sum():,}  ·  FSD={m.in_FSD.sum():,}  ·  neither={((~m.in_FSU)&(~m.in_FSD)).sum():,}')
# legend
from matplotlib.lines import Line2D
ax.legend([Line2D([0],[0],marker='o',color='w',markerfacecolor='#c44e52',markersize=9,label=f'FSU ({m.in_FSU.sum():,})'),
           Line2D([0],[0],marker='o',color='w',markerfacecolor='#4a7fbf',markersize=9,label=f'FSD ({m.in_FSD.sum():,})'),
           Line2D([0],[0],marker='o',color='w',markerfacecolor='#dcdcdc',markersize=9,label=f'not-selected ({((~m.in_FSU)&(~m.in_FSD)).sum():,})'),
           Line2D([0],[0],color='black',linestyle='--',label='≥2 tumors threshold')],
          [f'FSU ({m.in_FSU.sum():,})', f'FSD ({m.in_FSD.sum():,})',
           f'not-selected ({((~m.in_FSU)&(~m.in_FSD)).sum():,})', '≥2 tumors threshold'],
          loc='lower right', fontsize=9, frameon=True)
# annotate top-10 FSU + top-10 FSD
for _, r in m[m.in_FSU].nlargest(8, 'freq_maj').iterrows():
    ax.annotate(r.gene, (r.freq_maj, r.mean_log2FC), fontsize=7, color='#7a2b2c',
                xytext=(3, 3), textcoords='offset points')
for _, r in m[m.in_FSD].nlargest(8, 'freq_maj').iterrows():
    ax.annotate(r.gene, (r.freq_maj, r.mean_log2FC), fontsize=7, color='#2a4870',
                xytext=(3, -8), textcoords='offset points')
fig.savefig(OUT/'D_metavote_consensus_scatter.pdf', dpi=180, bbox_inches='tight')
fig.savefig(OUT/'D_metavote_consensus_scatter.png', dpi=140, bbox_inches='tight')
plt.close(fig)
print(f'  saved D_metavote_consensus_scatter.png')

# ---- combined 4-panel figure ----
from matplotlib.image import imread
import os
imgs = [OUT/'A_metavote_heatmap.png', OUT/'B_metavote_ntumor_hist.png',
        OUT/'C_metavote_waterfall.png', OUT/'D_metavote_consensus_scatter.png']
sizes = [imread(str(p)).shape for p in imgs]
# Just save individual + combined via subplot mosaic
fig = plt.figure(figsize=(16, 18), constrained_layout=True)
gs = fig.add_gridspec(3, 2)
for ax_pos, im_path, letter in [((slice(0,2), 0), imgs[0], 'A'),
                                  ((0, 1),  imgs[1], 'B'),
                                  ((1, 1),  imgs[3], 'D'),
                                  ((2, slice(0,2)), imgs[2], 'C')]:
    ax = fig.add_subplot(gs[ax_pos])
    ax.imshow(imread(str(im_path)))
    ax.set_axis_off()
fig.suptitle('Yang 2022 KP-Tracer style Meta-vote panels for Npp53 FitnessSignature construction', fontsize=13, fontweight='bold')
fig.savefig(OUT/'metavote_combined_ABCD.pdf', dpi=140, bbox_inches='tight')
fig.savefig(OUT/'metavote_combined_ABCD.png', dpi=110, bbox_inches='tight')
plt.close(fig)
print(f'  saved metavote_combined_ABCD.png')
print(f'\nDone.  Outputs → {OUT}')

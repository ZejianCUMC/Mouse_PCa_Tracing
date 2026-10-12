#!/usr/bin/env python
"""Yang 2022 Fig-2C-style meta-vote sign-consistency scatter for Npp53 FitnessSignature.

X = sign_consistency = up - dn  (right wing = consensus HIGH fitness)
Y = n_tumors_DE = freq_total
Colors:
    grey  = non-consensus / not in either signature
    red   = FSU (right wing, "Consensus high fitness")
    blue  = FSD (left  wing, "Consensus low  fitness")

Style follows REFERENCE_yang2022_fig2c_style_notes.md.
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['font.family']  = 'sans-serif'

# ── paths ────────────────────────────────────────────────────────────────────
FS  = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/'
           '1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/'
           'FitnessSignature')
OUT = FS / 'meta_vote_yang2022_plots'
OUT.mkdir(parents=True, exist_ok=True)

# ── data ─────────────────────────────────────────────────────────────────────
# Build meta from the FULL per-gene per-sample DE table so the "rejected but considered" grey
# cluster (n_tumors_DE >= 1, not in FSU/FSD) is visible at the plot's center-bottom.
# Yang 2022 exact per-tumor DE thresholds (KPTracer-release/reproducibility/Figure3_S3/Figure3_S3_diffexp.ipynb)
FDR_CUT, LFC_CUT = 0.10, np.log2(1.25)   # T4 pipeline per-tumor thresholds (matches modules + M3_merged)
FSU_set = {g.strip() for g in (FS/'FitnessPosSignature_All_mouse.txt').read_text().splitlines() if g.strip()}
FSD_set = {g.strip() for g in (FS/'FitnessNegSignature_All_mouse.txt').read_text().splitlines() if g.strip()}

pergene = pd.read_csv(FS/'FitnessSignature_per_gene_per_sample.csv.gz')
pergene = pergene[pergene['padj'] < FDR_CUT]
pergene = pergene[pergene['log2FC'].abs() > LFC_CUT]
up = pergene[pergene['log2FC'] > 0].groupby('gene').size()
dn = pergene[pergene['log2FC'] < 0].groupby('gene').size()
meta = pd.DataFrame({'up': up, 'dn': dn}).fillna(0).astype(int).reset_index()
meta['sign_consistency'] = meta['up'] - meta['dn']
meta['n_tumors_DE']      = meta['up'] + meta['dn']
meta = meta[meta['n_tumors_DE'] >= 1]
# borrow mean_log2FC (for label ranking); recompute per gene for genes not in original meta.csv
lfc_avg = pergene.groupby('gene')['log2FC'].mean().rename('mean_log2FC')
freq_maj = np.maximum(meta['up'].values, meta['dn'].values)
meta['freq_maj']    = freq_maj
meta['mean_log2FC'] = meta['gene'].map(lfc_avg)

# --- Yang 2022 EXACT classification (from KPTracer-release Figure3_S3_diffexp.ipynb) ---
#   SIGN_CONSISTENCY = 2       # |up - dn| >= 2
#   SIGN_TO_FREQ     = 0.5     # |up - dn| / max(1, n_tumors_DE) > 0.5
SIGN_CONSISTENCY = 2
SIGN_TO_FREQ     = 0.5
meta['sign_to_freq'] = meta['sign_consistency'].abs() / meta['n_tumors_DE'].clip(lower=1)
def _cls(row):
    if abs(row['sign_consistency']) < SIGN_CONSISTENCY: return 'other'
    if row['sign_to_freq'] <= SIGN_TO_FREQ: return 'other'
    if row['sign_consistency'] > 0: return 'FSU'
    if row['sign_consistency'] < 0: return 'FSD'
    return 'other'
meta['cls'] = meta.apply(_cls, axis=1)

n_up = int((meta.cls == 'FSU').sum())
n_dn = int((meta.cls == 'FSD').sum())
n_ot = int((meta.cls == 'other').sum())

# ── figure geometry ──────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5.5, 5.5))

# integer positions (no jitter)
sub_other = meta[meta.cls == 'other']
sub_fsu   = meta[meta.cls == 'FSU']
sub_fsd   = meta[meta.cls == 'FSD']

# 2-axis jitter — without it 1,494 genes collapse onto ~35 integer bins and the V-shape is invisible
rng = np.random.default_rng(0)
def _xj(n): return rng.uniform(-0.30, 0.30, size=n)
def _yj(n): return rng.uniform(-0.30, 0.30, size=n)

ax.scatter(sub_other['sign_consistency'] + _xj(len(sub_other)),
           sub_other['n_tumors_DE']      + _yj(len(sub_other)),
           c='#c9ccd1', s=5, alpha=0.5,  edgecolors='none', rasterized=True)
ax.scatter(sub_fsd['sign_consistency']   + _xj(len(sub_fsd)),
           sub_fsd['n_tumors_DE']        + _yj(len(sub_fsd)),
           c='#4a7fbf', s=6, alpha=0.75, edgecolors='none', rasterized=True)
ax.scatter(sub_fsu['sign_consistency']   + _xj(len(sub_fsu)),
           sub_fsu['n_tumors_DE']        + _yj(len(sub_fsu)),
           c='#c44e52', s=6, alpha=0.75, edgecolors='none', rasterized=True)

# ── limits ───────────────────────────────────────────────────────────────────
xmax = int(np.ceil(meta['sign_consistency'].abs().max())) + 1   # +1 tick padding
ymax = int(meta['n_tumors_DE'].max())                            # no padding
ax.set_xlim(-xmax, xmax)
ax.set_ylim(0, ymax)

# ── spines + labels ──────────────────────────────────────────────────────────
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.axvline(0, color='#888', linewidth=0.6, alpha=0.5)
ax.set_xlabel('Sign consistency', fontsize=11)
ax.set_ylabel('Number of tumors as\ndifferentially expressed', fontsize=11)
ax.tick_params(axis='both', labelsize=10)

# ── top labels (ABOVE plot area, via fig.text) ───────────────────────────────
# Compute axes bounding box in figure coords AFTER layout would be ideal,
# but we can position relative to the axes with axes-fraction transforms + offset.
ax.text(-0.5 * xmax, ymax * 1.03,
        f'Consensus low fitness\ngene (N={n_dn:,})',
        ha='center', va='bottom', fontsize=10, color='#000',
        transform=ax.transData, clip_on=False)
ax.text( 0.5 * xmax, ymax * 1.03,
        f'Consensus high fitness\ngene (N={n_up:,})',
        ha='center', va='bottom', fontsize=10, color='#000',
        transform=ax.transData, clip_on=False)

# (labeled boxes + representative gene names removed per user request — plot area kept clean)
top_up = sub_fsu.sort_values(['freq_maj', 'mean_log2FC'],
                             ascending=[False, False]).head(6)['gene'].tolist()
top_dn = sub_fsd.sort_values(['freq_maj', 'mean_log2FC'],
                             ascending=[False, True]).head(6)['gene'].tolist()

plt.tight_layout()

# ── save ─────────────────────────────────────────────────────────────────────
pdf = OUT / 'Yang2022_signconsistency.pdf'
png = OUT / 'Yang2022_signconsistency.png'
plt.savefig(pdf, dpi=300, bbox_inches='tight')
plt.savefig(png, dpi=300, bbox_inches='tight')
plt.close()

print(f'saved: {pdf}')
print(f'saved: {png}')
print(f'FSU={n_up:,} · FSD={n_dn:,} · other={n_ot:,}')
print(f'xmax={xmax}  ymax={ymax}')
print(f'top FSU (right, red):  {top_up}')
print(f'top FSD (left,  blue): {top_dn}')

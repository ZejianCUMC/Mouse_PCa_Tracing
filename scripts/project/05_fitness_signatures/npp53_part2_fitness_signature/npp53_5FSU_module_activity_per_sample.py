#!/usr/bin/env python
"""
Fig 9.6 — Per-sample activity of the 5 FSU modules (25 tumors × 5 modules).

Input : analysis/fitness_signatures_survival/threshold_sweep/T4/merged5/percell_scores_merged5.csv
        (per-cell AUCell M1..M5 already computed by T4_merged5_generate.py)

Output:
  analysis/fitness_signatures_survival/
      npp53_5FSU_module_activity_per_sample.png
      npp53_5FSU_module_activity_per_sample.pdf
      npp53_5FSU_module_activity_per_sample.tsv           (samples x modules, per-sample-mean AUCell)
      npp53_5FSU_module_activity_per_sample_zcol.tsv       (samples x modules, column z-score)

Method
------
- per_sample_activity = mean_over_cells(AUCell) — one number per (sample, module).
  AUCell is a per-cell score, so cell-count-weighted mean is the standard summary.
- Between-module scale differences are large (M2 ~0.4, M1 ~0.07); to make between-module
  variation across samples visually comparable, plot the column z-score (per module,
  centre & scale across the 25 samples). Raw matrix is also saved.
- Rows grouped by treatment arm (Intact / Castration / Cas+Enza) and sorted within-arm
  by M3 (the FSU-signature adverse-prognosis module — highest cross-cohort KM signal,
  see §9.4 Table 9.1) so the ordering reflects the biological hypothesis of interest.
"""
import warnings; warnings.filterwarnings('ignore')
import shutil, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['pdf.fonttype']  = 42
mpl.rcParams['ps.fonttype']   = 42
mpl.rcParams['font.family']   = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial']

ROOT = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
SRC  = ROOT/'analysis/fitness_signatures_survival/threshold_sweep/T4/merged5/percell_scores_merged5.csv'
OUT  = ROOT/'analysis/fitness_signatures_survival'
DROP = Path('/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/fitness_signatures_survival')
DROP.mkdir(parents=True, exist_ok=True)

TP_COLORS = {'Intact':'#3d7dcc', 'Castration':'#e1a83d', 'Castration + Enzalutamide':'#c3403a'}
TP_ORDER  = ['Intact','Castration','Castration + Enzalutamide']

MODS = [f'M{i}_AUCell' for i in range(1,6)]
MOD_LABELS = {
    'M1_AUCell':'M1\nStromal /\nECM-myCAF',
    'M2_AUCell':'M2\nGlycolysis /\nhypoxia',
    'M3_AUCell':'M3  (FSU winner)\nProliferation /\nsplicing',
    'M4_AUCell':'M4\nLum-Squamous /\nkeratin',
    'M5_AUCell':'M5\nMembrane-lipid /\nstress',
}
SORT_MOD = 'M3_AUCell'    # rank within each arm by M3 (FSU-signature adverse-prognosis)

# ── 1. load per-cell scores + aggregate to per-sample mean AUCell ────────────
print(f'[{time.strftime("%H:%M:%S")}] loading per-cell FSU-5 AUCell ...')
df = pd.read_csv(SRC)
print(f'  {len(df):,} cells; sample× treatment counts:')
print(df.groupby(['treatment','jz_sample']).size().unstack(fill_value=0))

per_sample = df.groupby('jz_sample')[MODS].mean()                     # (25, 5)
tp_by_samp = df.groupby('jz_sample')['treatment'].agg(lambda x: x.value_counts().index[0])
per_sample['treatment'] = tp_by_samp.reindex(per_sample.index)
print(f'  per-sample activity matrix: {per_sample.shape}')

per_sample.to_csv(OUT/'npp53_5FSU_module_activity_per_sample.tsv', sep='\t')
print(f'  wrote npp53_5FSU_module_activity_per_sample.tsv')

# ── 2. column z-score across the 25 samples ─────────────────────────────────
mat  = per_sample[MODS].values.astype(float)
mu   = mat.mean(axis=0, keepdims=True)
sd   = mat.std(axis=0, keepdims=True) + 1e-9
matZ = (mat - mu) / sd
zdf  = pd.DataFrame(matZ, index=per_sample.index, columns=MODS)
zdf['treatment'] = per_sample['treatment']
zdf.to_csv(OUT/'npp53_5FSU_module_activity_per_sample_zcol.tsv', sep='\t')

# ── 3. row order: treatment arm, then M3 z-score descending within arm ──────
row_meta = pd.DataFrame({'sample':per_sample.index,
                         'treatment':per_sample['treatment'].values,
                         'M3z':matZ[:, MODS.index(SORT_MOD)]})
row_meta['tp_rank'] = row_meta['treatment'].map({t:i for i,t in enumerate(TP_ORDER)})
row_meta = row_meta.sort_values(['tp_rank','M3z'], ascending=[True, False]).reset_index(drop=True)
row_order = row_meta['sample'].tolist()

matZ_ord = zdf.loc[row_order, MODS].values

# ── 4. figure ───────────────────────────────────────────────────────────────
print(f'[{time.strftime("%H:%M:%S")}] rendering ...')
fig = plt.figure(figsize=(8.2, 9.0))
# GridSpec: left arm-strip | main heatmap | colorbar
gs = fig.add_gridspec(nrows=1, ncols=3,
                      width_ratios=[0.22, 5.0, 0.22],
                      wspace=0.05,
                      left=0.11, right=0.94, top=0.90, bottom=0.19)
ax_left = fig.add_subplot(gs[0,0])
ax_main = fig.add_subplot(gs[0,1])
ax_cbar = fig.add_subplot(gs[0,2])

lim = float(np.nanmax(np.abs(matZ_ord)))
lim = min(lim, 2.5)
im  = ax_main.imshow(matZ_ord, aspect='auto', cmap='RdBu_r',
                     vmin=-lim, vmax=lim, interpolation='nearest')

# annotate each cell with the z-score
for i in range(matZ_ord.shape[0]):
    for j in range(matZ_ord.shape[1]):
        v = matZ_ord[i, j]
        txtc = 'white' if abs(v) > 1.2 else '#222'
        ax_main.text(j, i, f'{v:+.1f}', ha='center', va='center',
                     fontsize=7, color=txtc)

# module labels along x
ax_main.set_xticks(range(len(MODS)))
ax_main.set_xticklabels([MOD_LABELS[m] for m in MODS],
                        fontsize=8.0, rotation=0, ha='center')
ax_main.tick_params(axis='x', length=2, pad=6)

# no y ticks on the main axis (sample labels come from the left arm-strip axis)
ax_main.set_yticks([])
for spine in ax_main.spines.values():
    spine.set_linewidth(0.6); spine.set_edgecolor('#333')

# horizontal separators between arms
row_arm = row_meta['treatment'].tolist()
for i in range(1, len(row_arm)):
    if row_arm[i] != row_arm[i-1]:
        ax_main.axhline(i - 0.5, color='#111', lw=1.2)

# left arm strip + sample labels
arm_rgb = np.array([mpl.colors.to_rgb(TP_COLORS[t]) for t in row_arm])
ax_left.imshow(arm_rgb.reshape(-1, 1, 3), aspect='auto', interpolation='nearest')
ax_left.set_xticks([])
ax_left.set_yticks(range(len(row_order)))
ax_left.set_yticklabels(row_order, fontsize=8)
ax_left.tick_params(axis='y', length=0, pad=2)
for spine in ax_left.spines.values():
    spine.set_linewidth(0.6); spine.set_edgecolor('#333')

# colorbar
cbar = fig.colorbar(im, cax=ax_cbar)
cbar.set_label('per-sample AUCell (column z-score)', fontsize=9)
cbar.ax.tick_params(labelsize=8)

# arm legend below
leg_ax = fig.add_axes([0.05, 0.03, 0.92, 0.04])
leg_ax.axis('off')
counts = pd.Series(row_arm).value_counts()
handles = [mpl.patches.Patch(color=TP_COLORS[t], label=f'{t} (n={counts.get(t,0)})') for t in TP_ORDER]
leg_ax.legend(handles=handles, loc='center', frameon=False, ncol=3, fontsize=9)

fig.suptitle('Per-tumor activity of the 5 FSU modules  ·  25 tumors × 5 Hotspot modules  ·  column z-score  ·  within-arm sort by M3',
             x=0.5, y=0.965, fontsize=9.5, fontweight='bold', color='#0f2747')

png = OUT/'npp53_5FSU_module_activity_per_sample.png'
pdf = OUT/'npp53_5FSU_module_activity_per_sample.pdf'
fig.savefig(png, dpi=200, bbox_inches='tight')
fig.savefig(pdf, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'  wrote {png.name}  ({png.stat().st_size/1024:.0f} KB)')
print(f'  wrote {pdf.name}  ({pdf.stat().st_size/1024:.0f} KB)')

# ── 5. Dropbox sync ────────────────────────────────────────────────────────
for f in (png, pdf,
          OUT/'npp53_5FSU_module_activity_per_sample.tsv',
          OUT/'npp53_5FSU_module_activity_per_sample_zcol.tsv'):
    shutil.copy2(f, DROP/f.name)
print(f'[{time.strftime("%H:%M:%S")}] synced to {DROP}')

# ── 6. quick summary: per-arm top-2 modules (mean z) ───────────────────────
print()
print('=== per-arm top-2 modules (mean column-z) ===')
for tp in TP_ORDER:
    samples_in_arm = row_meta.loc[row_meta['treatment']==tp, 'sample'].tolist()
    arm_mean_z = zdf.loc[samples_in_arm, MODS].mean(axis=0).sort_values(ascending=False)
    print(f'  {tp:26s} (n={len(samples_in_arm)}) top 2:')
    for m, v in arm_mean_z.head(2).items():
        print(f'     {m:12s}  z={v:+.2f}')
print()
print(f'[{time.strftime("%H:%M:%S")}] DONE.')

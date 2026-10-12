#!/usr/bin/env python
"""
Reorganize T4 (9 modules) into 5 merged modules based on Hotspot local-correlation
adjacency (heatmap-derived), then regenerate:
  (0) new gene_module_assignment_merged5.csv + module_{1..5}_{genes,human}.txt
  (a) hotspot local-corr heatmap (5-module coloring)
  (b) 5x3 joint UMAP with AR panels + M1..M5 AUCell (inferno_r)
  (c) per-tumor pseudobulk AUCell heatmap (row-z)
  (d) per-state z + per-treatment z CSVs
  (e) per-cell scores CSV
  (f) module functional annotation table (csv + md)
Merge plan:
  M1_merged  = old M1 + M2 + M3       (heatmap block 1)
  M2_merged  = old M4                   (standalone)
  M3_merged  = old M5 + M6              (heatmap block; user preference)
  M4_merged  = old M7                   (standalone)
  M5_merged  = old M8 + M9              (heatmap block)
Outputs go into T4/merged5/ — nothing in T4/ is overwritten.
"""
import warnings; warnings.filterwarnings('ignore')
import re, time, json
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt, matplotlib as mpl
from scipy.cluster.hierarchy import linkage, leaves_list

mpl.rcParams['pdf.fonttype']=42; mpl.rcParams['ps.fonttype']=42
mpl.rcParams['font.family']='sans-serif'

ROOT = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
T4   = ROOT/'analysis/fitness_signatures_survival/threshold_sweep/T4'
OUT  = T4/'merged5'
OUT.mkdir(parents=True, exist_ok=True)
H5AD = ROOT/'data/npp53_merge/npp53_edited.h5ad'
SCFP = ROOT/'analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
AR_SCORES = ROOT/'analysis/aucell_arhi/aucell_scores.csv'
HOM_RPT = Path('/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/22Science_Tang_CRPC_ATAC/resources/orthology/HOM_MouseHumanSequence.rpt')

# ── MERGE PLAN ────────────────────────────────────────────────────────────────
MERGE_PLAN = {
    1: [1, 2, 3],
    2: [4],
    3: [5, 6],
    4: [7],
    5: [8, 9],
}
NEW_IDS = sorted(MERGE_PLAN.keys())
OLD_TO_NEW = {old_id: new_id for new_id, olds in MERGE_PLAN.items() for old_id in olds}
print(f'MERGE_PLAN: {MERGE_PLAN}')
print(f'OLD_TO_NEW: {OLD_TO_NEW}')

# ── STEP 1: build new gene assignment + per-module gene lists ───────────────
print('\n[1] Building merged5 gene assignment ...')
ga = pd.read_csv(T4/'gene_module_assignment.csv', index_col=0)
print(f'  original ga: {ga.shape}, old modules: {sorted(ga["module"].unique().tolist())}')
ga_merged = ga.copy()
ga_merged['module'] = ga['module'].map(lambda x: OLD_TO_NEW.get(int(x), -1) if x > 0 else -1)
ga_merged.to_csv(OUT/'gene_module_assignment_merged5.csv')
print(f'  saved gene_module_assignment_merged5.csv · new counts: {ga_merged["module"].value_counts().sort_index().to_dict()}')

# Build MGI mouse→human ortholog map
print('\nBuilding MGI mouse→human ortholog map ...')
hom = pd.read_csv(HOM_RPT, sep='\t', dtype=str)
hom.columns = [c.strip() for c in hom.columns]
mouse = hom[hom['Common Organism Name']=='mouse, laboratory'][['DB Class Key','Symbol']].rename(columns={'Symbol':'mouse'})
human = hom[hom['Common Organism Name']=='human'][['DB Class Key','Symbol']].rename(columns={'Symbol':'human'})
map_df = mouse.merge(human, on='DB Class Key')
# preserve first mapping per mouse gene
m2h = {}
for _, r in map_df.iterrows():
    m2h.setdefault(r['mouse'], r['human'])
print(f'  MGI map size: {len(m2h):,} mouse→human pairs')

# Concatenate old per-module mouse lists, preserving order
def read_mouse_mod(old_id):
    return [g.strip() for g in (T4/f'module_{old_id}_genes.txt').read_text().splitlines() if g.strip()]

MODS_M = {}
MODS_H = {}
for new_id in NEW_IDS:
    genes = []
    for old_id in MERGE_PLAN[new_id]:
        genes.extend(read_mouse_mod(old_id))
    # dedupe preserving order
    seen = set(); uniq = []
    for g in genes:
        if g not in seen:
            uniq.append(g); seen.add(g)
    MODS_M[new_id] = uniq
    # human orthologs via MGI
    h_seen = set(); h_uniq = []
    for g in uniq:
        h = m2h.get(g)
        if h and h not in h_seen:
            h_uniq.append(h); h_seen.add(h)
    MODS_H[new_id] = h_uniq
    (OUT/f'module_{new_id}_genes.txt').write_text('\n'.join(uniq) + '\n')
    (OUT/f'module_{new_id}_human.txt').write_text('\n'.join(h_uniq) + '\n')
    print(f'  M{new_id}_merged (old {MERGE_PLAN[new_id]}): n_mouse={len(uniq)}, n_human={len(h_uniq)}')

# ── STEP 3: HOTSPOT LOCAL-CORR HEATMAP (merged5 palette) ──────────────────────
print('\n[3] Building hotspot local-corr heatmap (merged5) ...')
lcz = pd.read_csv(T4/'local_correlation_z.csv', index_col=0)
common = [g for g in lcz.index if g in ga_merged.index]
lcz_c  = lcz.loc[common, common]
mods   = ga_merged.loc[common, 'module']
palette = plt.get_cmap('tab10')(np.arange(10) / 10)  # tab10 first 5 colors used
UNASSIGNED_COLOR = np.array([0.85, 0.85, 0.85, 1.0])

ordered_genes = []; mod_colors_row = []
for mid in NEW_IDS:
    mg = mods[mods == mid].index.tolist()
    if len(mg) > 1:
        try:
            order = leaves_list(linkage(1 - np.corrcoef(lcz_c.loc[mg, mg].values), method='average'))
            mg = [mg[i] for i in order]
        except Exception:
            pass
    ordered_genes.extend(mg)
    mod_colors_row.extend([palette[int(mid)-1]] * len(mg))
unassigned = mods[mods == -1].index.tolist()
if unassigned:
    ordered_genes.extend(unassigned)
    mod_colors_row.extend([UNASSIGNED_COLOR] * len(unassigned))
mat = lcz_c.loc[ordered_genes, ordered_genes].values
n = len(ordered_genes)

fig, (ax_bar, ax_hm) = plt.subplots(2, 1, figsize=(11, 12),
                                     gridspec_kw={'height_ratios':[0.4, 12]},
                                     constrained_layout=True)
ax_bar.imshow(np.array(mod_colors_row).reshape(1, n, -1), aspect='auto')
ax_bar.set_yticks([]); ax_bar.set_xticks([])
cursor = 0
for mid in NEW_IDS:
    block = int((mods == mid).sum())
    if block > 0:
        ax_bar.text(cursor + block/2 - 0.5, 0, f'M{int(mid)}', ha='center', va='center',
                    color='white', fontsize=11, fontweight='bold')
        cursor += block
if unassigned:
    ax_bar.text(cursor + len(unassigned)/2 - 0.5, 0, f'unassigned (n={len(unassigned)})',
                ha='center', va='center', color='black', fontsize=8, fontweight='bold')
im = ax_hm.imshow(mat, aspect='auto', cmap='RdBu_r', vmin=-8, vmax=8, interpolation='nearest')
ax_hm.set_xticks([]); ax_hm.set_yticks([])
ax_hm.set_xlabel(f'{n} genes (ordered by MERGED module 1..5, then {len(unassigned)} unassigned in grey)')
ax_hm.set_ylabel('genes')
plt.colorbar(im, ax=ax_hm, fraction=0.02, label='local-correlation Z')
ax_hm.set_title('Hotspot local-correlation heatmap · 5 MERGED modules (heatmap-adjacency-based)  ·  min_gene_threshold=20')
fp = OUT/'hotspot_local_corr_heatmap_merged5.pdf'
plt.savefig(fp, dpi=180, bbox_inches='tight')
plt.savefig(str(fp).replace('.pdf','.png'), dpi=140, bbox_inches='tight')
plt.close('all')
print(f'  saved {fp.name} + .png')

# ── STEP 2: LOAD H5AD + AUCELL ────────────────────────────────────────────────
print('\n[2/5] Loading h5ad + metadata for AUCell ...')
A = ad.read_h5ad(H5AD)
sc_fp = pd.read_csv(SCFP)
A.obs['cellBC_clean'] = [re.sub(r'^[A-Z0-9]+_', '', x).split('-')[0] for x in A.obs_names]
obs = A.obs.reset_index().rename(columns={'index':'obs_name'})
mrg = sc_fp.merge(obs[['obs_name','sample_id','cellBC_clean']],
                  left_on=['sample','cellBC'], right_on=['sample_id','cellBC_clean'], how='left')
mrg = mrg.drop_duplicates(subset='obs_name').set_index('obs_name')
A.obs['jz_sample'] = mrg.reindex(A.obs_names)['sample']
A.obs['fitness']   = mrg.reindex(A.obs_names)['fitness']
A.obs['plasticity']= mrg.reindex(A.obs_names)['plasticity']
A.obs['treatment'] = mrg.reindex(A.obs_names)['treatment']
A = A[A.obs['jz_sample'].notna(), :].copy()

scvi_cand = [k for k in A.obsm.keys() if 'scvi' in k.lower() or 'X_sc' in k.lower()]
if not scvi_cand:
    scvi_df = pd.read_csv(ROOT/'analysis/evolution_fitness_plasticity/shared_scvi_v10.csv', index_col=0)
    common_c = A.obs_names.intersection(scvi_df.index)
    A = A[common_c, :].copy()
    A.obsm['X_scVI'] = scvi_df.reindex(A.obs_names).values
else:
    A.obsm['X_scVI'] = A.obsm[scvi_cand[0]]
print(f'  {A.shape[0]:,} cells × {A.shape[1]:,} genes')

if 'X_umap' not in A.obsm:
    sc.pp.neighbors(A, use_rep='X_scVI', n_neighbors=30)
    sc.tl.umap(A)

A.obs['cell_type'] = A.obs['putative_state'].astype(str)
def major(x):
    if x.startswith('Luminal-'): return 'Luminal'
    if x.startswith('Mes-'):     return 'Mesenchymal'
    if x.startswith('Basal'):    return 'Basal'
    if x == 'Partial-EMT':       return 'Partial-EMT'
    if x == 'Cycling':           return 'Cycling'
    if x == 'Ionocyte':          return 'Ionocyte'
    return 'Other'
A.obs['major_cell_type'] = A.obs['cell_type'].map(major)

ar_df = pd.read_csv(AR_SCORES, index_col=0)
AR_COLS = ['Dehm_AR', 'HALLMARK_ANDROGEN_RESPONSE', 'Tang_CRPC_AR', 'Hieronymus_AR_signaling']
for col in AR_COLS:
    A.obs[col] = A.obs_names.map(ar_df[col])

print('\nComputing AUCell for M1..M5 merged ...')
if 'counts' in A.layers:
    X_raw = A.layers['counts']
else:
    X_raw = A.X
X = X_raw.toarray() if sp.issparse(X_raw) else X_raw
n_cells, n_genes = X.shape
top_n = int(0.05 * n_genes)
print(f'  ranking {n_cells:,} cells × {n_genes:,} genes; AUCell top-N = {top_n}')

t0 = time.time()
ranks_desc = (-X).argsort(axis=1).argsort(axis=1) + 1
print(f'  ranks done in {time.time()-t0:.1f}s')

var_index = A.var_names
def aucell_score(genes, ranks_desc, top_n):
    gene_idx = np.array([var_index.get_loc(g) for g in genes if g in var_index])
    n_set = len(gene_idx)
    if n_set == 0: return np.zeros(ranks_desc.shape[0])
    set_ranks = ranks_desc[:, gene_idx]
    sorted_ranks = np.sort(set_ranks, axis=1)
    auc = np.zeros(ranks_desc.shape[0])
    for i in range(ranks_desc.shape[0]):
        r = sorted_ranks[i]
        r_clip = np.minimum(r, top_n)
        auc[i] = ((top_n - r_clip).sum()) / (n_set * top_n)
    return auc

MODS_IN = {}
for i in NEW_IDS:
    genes = MODS_M[i]
    genes_in = [g for g in genes if g in var_index]
    MODS_IN[i] = genes_in
    t0 = time.time()
    A.obs[f'M{i}_AUCell'] = aucell_score(genes_in, ranks_desc, top_n)
    print(f'  M{i}_merged: {len(genes)} → {len(genes_in)} present · AUCell in {time.time()-t0:.1f}s')

# ── STEP 2c: per-state / per-treatment z-scores + per-cell CSV ────────────────
print('\nPer-state + per-treatment z-scores ...')
mods_cols = [f'M{i}_AUCell' for i in NEW_IDS]
pss_mean = A.obs.groupby('putative_state', observed=True)[mods_cols].mean().reset_index()
ptr_mean = A.obs.groupby('treatment',      observed=True)[mods_cols].mean().reset_index()
def zscore_cols(df, cols):
    d = df.copy()
    for c in cols:
        v = d[c].astype(float)
        d[c] = (v - v.mean()) / (v.std(ddof=0) + 1e-12)
    return d
pss_z = zscore_cols(pss_mean, mods_cols)
ptr_z = zscore_cols(ptr_mean, mods_cols)
pss_z.to_csv(OUT/'per_state_z_AUCell_merged5.csv', index=False)
ptr_z.to_csv(OUT/'per_treatment_z_AUCell_merged5.csv', index=False)
pss_mean.to_csv(OUT/'per_state_mean_AUCell_merged5.csv', index=False)
ptr_mean.to_csv(OUT/'per_treatment_mean_AUCell_merged5.csv', index=False)
print(f'  saved per_state_z, per_state_mean, per_treatment_z, per_treatment_mean')

sco = A.obs[['jz_sample','treatment','cell_type','major_cell_type',
             'fitness','plasticity',
             'Dehm_AR','HALLMARK_ANDROGEN_RESPONSE','Tang_CRPC_AR','Hieronymus_AR_signaling',
             ] + mods_cols].copy()
sco.to_csv(OUT/'percell_scores_merged5.csv')
print(f'  saved percell_scores_merged5.csv ({A.n_obs:,} rows)')

# ── STEP 5: 5×3 joint UMAP ────────────────────────────────────────────────────
print('\n[5] Building 5×3 joint UMAP ...')
fig, axes = plt.subplots(5, 3, figsize=(15, 21), constrained_layout=True)

# R1
ct = A.obs['cell_type'].astype(str)
ct_cats = sorted(ct.unique())
palette_ct = plt.get_cmap('tab20')(np.linspace(0, 1, len(ct_cats)))
ct_color = dict(zip(ct_cats, palette_ct))
ax = axes[0,0]
c = ct.map(ct_color).values
ax.scatter(A.obsm['X_umap'][:,0], A.obsm['X_umap'][:,1], s=1.2, c=c, rasterized=True)
ax.set_title(f'cell-type (n={A.n_obs:,})', fontsize=11); ax.set_xticks([]); ax.set_yticks([])
for lbl in ct_cats:
    ax.scatter([], [], s=15, c=[ct_color[lbl]], label=lbl)
ax.legend(fontsize=6, loc='upper right', frameon=False, ncol=2, markerscale=1.2, borderpad=0.2)

mct_cats = ['Luminal','Basal','Mesenchymal','Partial-EMT','Cycling','Ionocyte','Other']
mct_colors = {'Luminal':'#1f5fa8','Basal':'#c44e52','Mesenchymal':'#2ca02c',
              'Partial-EMT':'#E69F00','Cycling':'#8b3a8b','Ionocyte':'#a3a3a3','Other':'#e0e0e0'}
ax = axes[0,1]
c = A.obs['major_cell_type'].map(mct_colors).values
ax.scatter(A.obsm['X_umap'][:,0], A.obsm['X_umap'][:,1], s=1.2, c=c, rasterized=True)
ax.set_title('major-cell-type', fontsize=11); ax.set_xticks([]); ax.set_yticks([])
for lbl in mct_cats:
    n_l = (A.obs['major_cell_type']==lbl).sum()
    if n_l>0: ax.scatter([], [], s=15, c=mct_colors[lbl], label=f'{lbl} ({n_l:,})')
ax.legend(fontsize=8, loc='upper right', frameon=False)

tp_colors = {'Intact':'#3d7dcc', 'Castration':'#e1a83d', 'Castration + Enzalutamide':'#c3403a'}
ax = axes[0,2]
c = A.obs['treatment'].map(tp_colors).values
ax.scatter(A.obsm['X_umap'][:,0], A.obsm['X_umap'][:,1], s=1.2, c=c, rasterized=True)
ax.set_title('Treatment', fontsize=11); ax.set_xticks([]); ax.set_yticks([])
for lbl in ['Intact','Castration','Castration + Enzalutamide']:
    n_l=(A.obs['treatment']==lbl).sum()
    ax.scatter([], [], s=15, c=tp_colors[lbl], label=f'{lbl} ({n_l:,})')
ax.legend(fontsize=8, loc='upper right', frameon=False)

def plot_continuous(ax, values, title, cmap='RdBu_r'):
    v = values.astype(float)
    mask = np.isfinite(v)
    vmin, vmax = np.percentile(v[mask], [5, 95])
    sc_p = ax.scatter(A.obsm['X_umap'][mask,0], A.obsm['X_umap'][mask,1], s=1.2,
                      c=v[mask], cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
    if (~mask).any():
        ax.scatter(A.obsm['X_umap'][~mask,0], A.obsm['X_umap'][~mask,1], s=1.2,
                   c='#e5e5e5', rasterized=True, zorder=0)
    ax.set_title(title, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(sc_p, ax=ax, fraction=0.04)

# R2: jungle fitness | plasticity | HALLMARK_AR
plot_continuous(axes[1,0], A.obs['fitness'].values, 'jungle fitness (LBI)', 'RdBu_r')
plot_continuous(axes[1,1], A.obs['plasticity'].values, 'plasticity', 'viridis_r')
plot_continuous(axes[1,2], A.obs['HALLMARK_ANDROGEN_RESPONSE'].values,
                'HALLMARK_ANDROGEN_RESPONSE (AUCell)', 'RdBu_r')

# R3: Dehm_AR | Tang_CRPC_AR | Hieronymus_AR
plot_continuous(axes[2,0], A.obs['Dehm_AR'].values, 'Dehm_AR (AUCell)', 'RdBu_r')
plot_continuous(axes[2,1], A.obs['Tang_CRPC_AR'].values, 'Tang_CRPC_AR (AUCell)', 'RdBu_r')
plot_continuous(axes[2,2], A.obs['Hieronymus_AR_signaling'].values,
                'Hieronymus_AR_signaling (AUCell)', 'RdBu_r')

# R4: M1|M2|M3 ; R5: M4|M5|blank
def plot_module(ax, mid):
    v = A.obs[f'M{mid}_AUCell'].values
    vmin, vmax = np.percentile(v, [5, 95])
    sc_p = ax.scatter(A.obsm['X_umap'][:,0], A.obsm['X_umap'][:,1], s=1.2, c=v,
                      cmap='inferno_r', vmin=vmin, vmax=vmax, rasterized=True)
    ax.set_title(f'AUCell · M{mid}_merged (n={len(MODS_IN[mid])})', fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(sc_p, ax=ax, fraction=0.04)

plot_module(axes[3,0], 1); plot_module(axes[3,1], 2); plot_module(axes[3,2], 3)
plot_module(axes[4,0], 4); plot_module(axes[4,1], 5)
axes[4,2].axis('off')

fig.suptitle('Npp53 T4 (heatmap-merged 5 modules) · joint UMAP with AR panels', fontsize=14, y=1.01)
fp = OUT/'joint_UMAP_merged5.pdf'
fig.savefig(fp, dpi=180, bbox_inches='tight')
fig.savefig(str(fp).replace('.pdf','.png'), dpi=140, bbox_inches='tight')
plt.close(fig)
print(f'  saved {fp.name} + .png')

# ── STEP 4: pseudobulk AUCell heatmap ─────────────────────────────────────────
print('\n[4] Per-tumor pseudobulk AUCell heatmap ...')
pb = pd.DataFrame({'sample': A.obs['jz_sample'].values,
                   'treatment': A.obs['treatment'].values})
for i in NEW_IDS:
    pb[f'M{i}'] = A.obs[f'M{i}_AUCell'].values
pb = pb.groupby(['treatment','sample']).mean().reset_index()
pb_sorted = pb.sort_values(['treatment','M1'], ascending=[True, False]).reset_index(drop=True)
mods_cols_lbl = [f'M{i}' for i in NEW_IDS]
mat_pb = pb_sorted[mods_cols_lbl].values.T
mat_z = (mat_pb - mat_pb.mean(axis=1, keepdims=True)) / (mat_pb.std(axis=1, keepdims=True) + 1e-9)

fig, ax = plt.subplots(figsize=(max(6, len(pb_sorted)*0.35), 3.2), constrained_layout=True)
lim = min(2.5, np.max(np.abs(mat_z)))
im = ax.imshow(mat_z, aspect='auto', cmap='RdBu_r', vmin=-lim, vmax=lim, interpolation='nearest')
ax.set_yticks(range(len(NEW_IDS)))
ax.set_yticklabels([f'M{i}_merged (n={len(MODS_IN[i])})' for i in NEW_IDS])
ax.set_xticks(range(len(pb_sorted)))
ax.set_xticklabels([f'{r.treatment[:6]} · {r["sample"]}' for _, r in pb_sorted.iterrows()],
                   rotation=90, fontsize=7)
for i, r in pb_sorted.iterrows():
    ax.add_patch(mpl.patches.Rectangle((i-0.5, len(NEW_IDS)+0.4), 1, 0.4,
                 color=tp_colors[r['treatment']], clip_on=False))
plt.colorbar(im, ax=ax, fraction=0.03, label='AUCell z-score (row-scaled)')
ax.set_title('Per-tumor pseudobulk AUCell (row-z) · 5 MERGED modules (sorted within-treatment by M1)')
fp = OUT/'pseudobulk_AUCell_heatmap_merged5.pdf'
fig.savefig(fp, dpi=180, bbox_inches='tight')
fig.savefig(str(fp).replace('.pdf','.png'), dpi=140, bbox_inches='tight')
plt.close(fig)
print(f'  saved {fp.name} + .png')

# ── STEP 8: functional annotation table (deferred KM enrichment done later)  ─
# We save an intermediate table WITHOUT KM columns; second script fills KM.
print('\n[8-pre] Functional annotation (KM columns will be filled after KM run) ...')
autoc = pd.read_csv(T4/'autocorrelations.csv').set_index('Gene')
pss_z_i = pss_z.set_index('putative_state')
ptr_z_i = ptr_z.set_index('treatment')

# Hand-curated one-sentence interpretations (based on old-module composition)
INTERP = {
    1: 'Stromal / ECM-myCAF program (old M1+M2+M3) — collagens (Col3a1/5a1/6a1/6a2/8a1/16a1), Dcn, Ptn, Timp1/2, Myl9, Loxl2, Nes, Adam19, Msn, S100a6, Cd47, Nrp2, Adam12, Twist1.',
    2: 'Glycolysis + hypoxia + translation module (old M4) — Gapdh, Pgk1, Pgam1, Pfkl, Tpi1, Eno1, Bnip3, Egln1, Vegfa, ribosomal proteins.',
    3: 'Proliferation-housekeeping combined (old M5+M6) — G2M mitosis + splicing + proteostasis (Aspm, Bub1, Ccna2, Ezh2, Hells, Diaph3, Pclaf, Ska1, snRNPs, Psm subunits, Ran/Ranbp1) plus ubiquitin/transcription regulators (Ankrd12, Kansl1, Ubr2, Setd7, Ranbp9, Bcl2l11).',
    4: 'Luminal-Squamous / keratinization / stress-response (old M7) — Krt23, Krt80, Sprr2d, Cst6, Elf5, Rarb, Lgals3, Nupr1, S100a11, Il33/Il18r1, Akr1b3.',
    5: 'Membrane-lipid + chaperone + stress-cytokine module (old M8+M9) — Plin2, Cyb5a, Hsp90b1, Runx3, Lamtor2, Myl6/12b, Ost4, Cox14, plus Gdf15, Angptl4, Ephx1, Podxl, Lefty1, Ifngr1.',
}

rows_pre = []
for i in NEW_IDS:
    genes_m = MODS_M[i]; genes_h = MODS_H[i]
    z_series = autoc.loc[[g for g in genes_m if g in autoc.index], 'Z'].sort_values(ascending=False)
    top20 = z_series.head(20).index.tolist()
    col = f'M{i}_AUCell'
    if col in pss_z_i.columns:
        state_z = pss_z_i[col].sort_values(ascending=False)
        top_states = [f'{s} (z={state_z.loc[s]:+.2f})' for s in state_z.head(3).index]
    else:
        top_states = []
    if col in ptr_z_i.columns:
        tz = ptr_z_i[col].sort_values(ascending=False)
        top_tr_name = tz.index[0]
        top_tr = f'{top_tr_name} (z={tz.iloc[0]:+.2f})'
    else:
        top_tr = ''
    rows_pre.append({
        'module': f'M{i}_merged',
        'old_modules': ','.join(f'M{o}' for o in MERGE_PLAN[i]),
        'n_mouse': len(genes_m),
        'n_human': len(genes_h),
        'top20_genes_by_Z': ', '.join(top20),
        'top3_states': '; '.join(top_states),
        'top_treatment': top_tr,
        'interpretation': INTERP.get(i,''),
    })

pd.DataFrame(rows_pre).to_csv(OUT/'module_functional_table_merged5_pre_KM.csv', index=False)
print(f'  saved module_functional_table_merged5_pre_KM.csv (KM columns to be filled)')

# save autocorr Z with new module for downstream xlsx
long_rows = []
for i in NEW_IDS:
    mg = MODS_M[i]; hg = MODS_H[i]
    z_full = pd.Series({g: autoc.loc[g,'Z'] if g in autoc.index else np.nan for g in mg})
    fdr_full = pd.Series({g: autoc.loc[g,'FDR'] if g in autoc.index else np.nan for g in mg})
    order = z_full.sort_values(ascending=False).index.tolist()
    rank_map = {g: r+1 for r,g in enumerate(order)}
    # map mouse→human via MGI m2h
    for pos, g in enumerate(mg):
        h = m2h.get(g, '')
        long_rows.append({
            'module': f'M{i}_merged',
            'gene_mouse': g,
            'gene_human': h,
            'autocorr_Z': float(z_full[g]) if pd.notna(z_full[g]) else np.nan,
            'autocorr_FDR': float(fdr_full[g]) if pd.notna(fdr_full[g]) else np.nan,
            'rank_in_module': rank_map[g],
            'is_representative': rank_map[g] <= 30,
        })
pd.DataFrame(long_rows).to_csv(OUT/'all_genes_long_merged5.csv', index=False)
print(f'  saved all_genes_long_merged5.csv ({len(long_rows)} rows)')

print('\nAll Python assets (Step 1-5, 8-pre) generated.')

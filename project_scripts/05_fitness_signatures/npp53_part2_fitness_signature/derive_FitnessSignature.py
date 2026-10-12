#!/usr/bin/env python
"""Derive FitnessSignature following Yang et al. 2022 (Cell) methodology.

Pipeline (per Yang 2022 STAR methods, DOI 10.1016/j.cell.2022.04.015):
  1. Within each tumor sample separately, run per-gene linear regression:
         log(1 + expr_ij) ~ fitness_i + size_factor_i
     where size_factor_i = number of genes detected in cell i.
  2. LRT vs null model (size_factor only). FDR-BH correct per sample.
  3. log2FC = mean expression in top 10% fitness cells vs bottom 10% fitness cells.
  4. Retain genes expressed in >10 cells per sample.
  5. Meta-vote across samples: gene qualifies for FitnessSignature if
       (a) it is DE (padj<0.05 AND |log2FC|>log2(1.5)) in >=2 tumors
       (b) with consistent direction — ratio (frequency / consistency) > 0.5
     where consistency = #same-direction / total-DE-calls
  6. Split into FSU (up-with-fitness) and FSD (down-with-fitness).
  7. Score every cell with sc.tl.score_genes(FSU) and (FSD).
  8. Compute per-sample mean signature scores.

Vectorized OLS per sample (10-30s / sample; ~10 min total for 25 samples × ~20K genes).
"""
import warnings; warnings.filterwarnings('ignore')
import re, time
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy import stats
import anndata as ad
from statsmodels.stats.multitest import multipletests

M = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
H5AD  = M/'1-mouse_based/Npp53_tracing/data/npp53_merge/npp53_edited.h5ad'
SCFP  = M/'1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT   = M/'1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/FitnessSignature'
OUT.mkdir(parents=True, exist_ok=True)

MIN_CELLS_GENE = 10       # Yang: "genes appearing in more than 10 cells"
Q_TAIL         = 0.10     # Yang: "top vs bottom 10th percentile"
FDR_CUT        = 0.05
LFC_CUT        = np.log2(1.5)
META_MIN_FREQ  = 2        # Yang: "appeared to show up at least 2 times"
META_RATIO_MIN = 0.5      # Yang: "ratio between frequency and consistency was greater than 0.5"

# ─── 1. Load ────────────────────────────────────────────────────────────────
print(f'[{time.strftime("%H:%M:%S")}] Loading h5ad ...')
A = ad.read_h5ad(H5AD)
X_counts = A.layers['counts']   # sparse csr counts

sc_fp = pd.read_csv(SCFP)
# strip run prefix + -1 suffix on obs_names to align with sc_fp['cellBC']
A.obs['cellBC_clean'] = [re.sub(r'^[A-Z0-9]+_', '', x).split('-')[0] for x in A.obs_names]
obs = A.obs.reset_index().rename(columns={'index':'obs_name'})
merge = sc_fp.merge(obs[['obs_name','sample_id','cellBC_clean']],
                    left_on=['sample','cellBC'], right_on=['sample_id','cellBC_clean'], how='left')
assert merge['obs_name'].notna().all(), 'cell mapping failed'
merge = merge.drop_duplicates(subset='obs_name').set_index('obs_name')
# align to A.obs_names order
A.obs['fitness'] = merge.reindex(A.obs_names)['fitness']
A.obs['jz_sample'] = merge.reindex(A.obs_names)['sample']
n_ok = A.obs['fitness'].notna().sum()
print(f'  h5ad {A.shape}; cells with jungle fitness: {n_ok:,} / {A.shape[0]:,}')

# subset to jungle-scored cells
mask = A.obs['fitness'].notna().values
A = A[mask, :].copy()
X_counts = A.layers['counts'].tocsr()

# size_factor per Yang = number of genes detected per cell
size_factor = np.asarray((X_counts > 0).sum(axis=1)).ravel().astype(float)
A.obs['size_factor'] = size_factor
samples = sorted(A.obs['jz_sample'].unique())
print(f'  {len(samples)} tumor samples ready. Genes: {A.shape[1]:,}')

# ─── 2. Per-sample per-gene regression (vectorized) ─────────────────────────
def per_sample_regress(s):
    """Return DataFrame with gene, log2FC, p_lrt, padj, n_cells_expr for sample s."""
    sub_idx = np.where(A.obs['jz_sample'].values == s)[0]
    n = len(sub_idx)
    if n < 20: return None
    y_counts = X_counts[sub_idx, :]                             # cells × genes sparse
    f = A.obs['fitness'].values[sub_idx].astype(float)
    sf = A.obs['size_factor'].values[sub_idx].astype(float)
    # gene-detection filter — expressed in > MIN_CELLS_GENE
    n_cells_expr = np.asarray((y_counts > 0).sum(axis=0)).ravel()
    gene_mask = n_cells_expr > MIN_CELLS_GENE
    if gene_mask.sum() == 0: return None
    y = np.log1p(y_counts[:, gene_mask].toarray())              # (n cells, g genes)
    n_genes = y.shape[1]
    # design matrices
    X0 = np.column_stack([np.ones(n), sf])                      # null: intercept + size_factor
    X1 = np.column_stack([np.ones(n), sf, f])                   # full: + fitness
    # closed-form OLS: β = (XᵀX)⁻¹ Xᵀy  → residuals per gene
    def rss(X, Y):
        beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
        resid = Y - X @ beta
        return (resid**2).sum(axis=0), beta
    rss0, _   = rss(X0, y)
    rss1, b1  = rss(X1, y)
    dof0, dof1 = n - X0.shape[1], n - X1.shape[1]
    # F test comparing nested models (df numerator=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        F = ((rss0 - rss1) / 1) / (rss1 / dof1)
    p = 1.0 - stats.f.cdf(F, 1, dof1)
    p = np.where(np.isfinite(p), p, 1.0)
    # log2FC via top-10% vs bottom-10% fitness (in expression space; Yang used log-normalized expr)
    q_lo, q_hi = np.quantile(f, Q_TAIL), np.quantile(f, 1 - Q_TAIL)
    mask_lo, mask_hi = f <= q_lo, f >= q_hi
    if mask_lo.sum() < 2 or mask_hi.sum() < 2:
        lfc = np.full(n_genes, np.nan)
    else:
        # normalize to median UMI count (like Yang) then log2 mean+1
        counts_lo = y_counts[mask_lo, :][:, gene_mask].toarray()
        counts_hi = y_counts[mask_hi, :][:, gene_mask].toarray()
        med_umi = np.median(sf) if np.median(sf) > 0 else 1.0
        norm_lo = counts_lo / (sf[mask_lo, None] / med_umi + 1e-9)
        norm_hi = counts_hi / (sf[mask_hi, None] / med_umi + 1e-9)
        lfc = np.log2(norm_hi.mean(axis=0) + 1e-9) - np.log2(norm_lo.mean(axis=0) + 1e-9)
    # FDR-BH per sample
    _, padj, _, _ = multipletests(p, method='fdr_bh')
    return pd.DataFrame({
        'gene': A.var_names[gene_mask],
        'log2FC': lfc,
        'p': p, 'padj': padj,
        'n_cells_expr': n_cells_expr[gene_mask],
        'sample': s,
    })

print(f'[{time.strftime("%H:%M:%S")}] Running per-sample regression on {len(samples)} tumors ...')
per_sample_dfs = []
for i, s in enumerate(samples):
    t0 = time.time()
    df = per_sample_regress(s)
    if df is not None:
        per_sample_dfs.append(df)
        n_de = int(((df['padj']<FDR_CUT)&(df['log2FC'].abs()>LFC_CUT)).sum())
        print(f'  [{i+1:2d}/{len(samples)}] {s}: n_cells={int((A.obs["jz_sample"]==s).sum())}, '
              f'genes_tested={len(df):5d}, DE(FDR<{FDR_CUT}, |lfc|>{LFC_CUT:.2f})={n_de:4d} '
              f'({time.time()-t0:.1f}s)')
per_gene_all = pd.concat(per_sample_dfs, ignore_index=True)
per_gene_all.to_csv(OUT/'FitnessSignature_per_gene_per_sample.csv.gz', index=False, compression='gzip')

# ─── 3. Meta-vote (Yang 2022 STAR "Meta-analysis and derivation of the FitnessSignature") ─
print(f'\n[{time.strftime("%H:%M:%S")}] Meta-vote across tumors ...')
de = per_gene_all[(per_gene_all['padj']<FDR_CUT) & (per_gene_all['log2FC'].abs()>LFC_CUT)].copy()
de['direction'] = np.where(de['log2FC'] > 0, 'up', 'dn')
grp = de.groupby(['gene','direction']).size().unstack(fill_value=0)
grp['freq_total'] = grp.sum(axis=1)
grp['majority_dir'] = np.where(grp['up'] >= grp['dn'], 'up', 'dn')
grp['consistency'] = np.maximum(grp['up'], grp['dn']) / grp['freq_total']
# Yang criterion: freq >= 2 in the majority direction; ratio(freq/consistency) > 0.5
grp['freq_maj'] = np.where(grp['majority_dir']=='up', grp['up'], grp['dn'])
grp['ratio']    = grp['freq_maj'] / grp['consistency']
signature_hit = (grp['freq_maj'] >= META_MIN_FREQ) & (grp['ratio'] > META_RATIO_MIN)
sig = grp[signature_hit].copy()
sig['mean_log2FC'] = de.groupby('gene')['log2FC'].mean().reindex(sig.index).values
sig.to_csv(OUT/'FitnessSignature_meta.csv')

FSU = sig[sig['majority_dir']=='up'].index.tolist()
FSD = sig[sig['majority_dir']=='dn'].index.tolist()
Path(OUT/'FitnessPosSignature_mouse.txt').write_text('\n'.join(FSU))
Path(OUT/'FitnessNegSignature_mouse.txt').write_text('\n'.join(FSD))
print(f'  FitnessSignature: FSU={len(FSU)} genes, FSD={len(FSD)} genes; total {len(sig)}')

# ─── 4. Score cells + per-sample mean ───────────────────────────────────────
print(f'\n[{time.strftime("%H:%M:%S")}] Scoring cells with sc.tl.score_genes ...')
# reload as non-backed AnnData with log-normalized in .X for scoring
A2 = A.copy()
A2.X = X_counts.astype(np.float32)
sc.pp.normalize_total(A2, target_sum=1e4)
sc.pp.log1p(A2)
sc.tl.score_genes(A2, gene_list=[g for g in FSU if g in A2.var_names], score_name='FSU_score', use_raw=False)
sc.tl.score_genes(A2, gene_list=[g for g in FSD if g in A2.var_names], score_name='FSD_score', use_raw=False)
A2.obs['FS_diff'] = A2.obs['FSU_score'] - A2.obs['FSD_score']

per_cell = A2.obs[['jz_sample','fitness','FSU_score','FSD_score','FS_diff']].copy()
per_cell.index.name = 'obs_name'
per_cell.to_csv(OUT/'FitnessSignature_per_cell_scores.csv.gz', compression='gzip')

per_sample = (per_cell.groupby('jz_sample')
              .agg(ncells=('fitness','count'),
                   mean_FSU=('FSU_score','mean'),
                   mean_FSD=('FSD_score','mean'),
                   mean_FS_diff=('FS_diff','mean'),
                   mean_jungle=('fitness','mean'),
                   median_jungle=('fitness','median'))
              .reset_index().rename(columns={'jz_sample':'sample'}))
# n_cells confound check for the new scalar
def rho(a, b):
    ok = np.isfinite(a)&np.isfinite(b)
    return stats.spearmanr(a[ok], b[ok])
print('\n── FitnessSignature n_cells confound check ──')
for col in ['mean_FSU','mean_FSD','mean_FS_diff']:
    r, p = rho(np.log10(per_sample['ncells'].astype(float)), per_sample[col].astype(float))
    flag = ' ✓' if abs(r)<0.30 else (' ⚠' if abs(r)<0.55 else ' ✗')
    print(f'  Spearman ρ(log10 n_cells, {col:14s}) = {r:+.2f}  p={p:.3g}{flag}')

per_sample.to_csv(OUT/'FitnessSignature_per_sample_summary.csv', index=False)
print(f'\n[{time.strftime("%H:%M:%S")}] DONE. Outputs in {OUT}')
print(per_sample.to_string(index=False))

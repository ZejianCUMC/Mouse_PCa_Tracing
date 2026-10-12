#!/usr/bin/env python
"""4-way threshold sweep for M8-like module discovery.

Steps per threshold combo (T1..T4):
  1. Meta-vote FSU/FSD from the pre-computed per-gene per-sample DE table
     (reuse FitnessSignature_per_gene_per_sample.csv.gz — 25 tumors × 33,676 genes).
     Apply noNC regex filter + MGI ortholog conversion.
  2. Run Hotspot on FSU (h5ad + scVI latent, danb, k=30, top-500 sig-by-Z,
     create_modules(min_gene_threshold=20, core_only=True, fdr_threshold=0.05)).
  3. Compute per-cell AUCell for each module → per-state, per-treatment z-scores,
     and Cas/CasEnz-enriched flags (z > 0.5).

Threshold table:
  T1 : FDR<0.05, |LFC|>log2(1.25), |sign_consistency|>=2, no consistency filter
  T2 : same + sign_to_freq > 0.5
  T3 : FDR<0.10, |LFC|>log2(1.25), |sign_consistency|>=2, no consistency filter
  T4 : same + sign_to_freq > 0.5

  sign_consistency = up - dn (net directional agreement across per-tumor DE calls)
  sign_to_freq     = |sign_consistency| / max(1, n_tumors_DE)
"""
import warnings; warnings.filterwarnings('ignore')
import re, time, json, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, hotspot

M   = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
FS  = M/'analysis/fitness_signatures_survival'
SWEEP = FS/'threshold_sweep'

PER_GENE_FP = FS/'FitnessSignature/FitnessSignature_per_gene_per_sample.csv.gz'
HOM_FP      = '/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/22Science_Tang_CRPC_ATAC/resources/orthology/HOM_MouseHumanSequence.rpt'

# -------------------- threshold combos --------------------
THRESHOLDS = {
    'T1': dict(fdr=0.05, lfc=np.log2(1.25), min_abs_sign=2, use_sign_to_freq=False),
    'T2': dict(fdr=0.05, lfc=np.log2(1.25), min_abs_sign=2, use_sign_to_freq=True),
    'T3': dict(fdr=0.10, lfc=np.log2(1.25), min_abs_sign=2, use_sign_to_freq=False),
    'T4': dict(fdr=0.10, lfc=np.log2(1.25), min_abs_sign=2, use_sign_to_freq=True),
}

# noNC regex (per task): Rik$, ^Gm\d+, ^mt-, ^[A-Z]{2}\d{6}$
NC_PATTERNS = [re.compile(r'Rik$'),
               re.compile(r'^Gm\d+'),
               re.compile(r'^mt-'),
               re.compile(r'^[A-Z]{2}\d{6}$')]
def is_nc(g):
    return any(p.search(g) for p in NC_PATTERNS)

def stamp(): return time.strftime('%H:%M:%S')
def log(msg): print(f'[{stamp()}] {msg}', flush=True)

# =========================================================
# LOAD SHARED DATA (once)
# =========================================================
log('=== Loading per-gene per-sample DE table ===')
PER = pd.read_csv(PER_GENE_FP)
log(f'  rows={len(PER):,}  tumors={PER["sample"].nunique()}  genes={PER["gene"].nunique()}')

log('=== Loading MGI HOM mouse→human table ===')
hom = pd.read_csv(HOM_FP, sep='\t')
mouse_h = hom[hom['Common Organism Name']=='mouse, laboratory'][['DB Class Key','Symbol']].rename(columns={'Symbol':'mouse_symbol'})
human_h = hom[hom['Common Organism Name']=='human']            [['DB Class Key','Symbol']].rename(columns={'Symbol':'human_symbol'})
mh = mouse_h.merge(human_h, on='DB Class Key').drop_duplicates('mouse_symbol')
def to_human(mouse_genes):
    return mh[mh.mouse_symbol.isin(mouse_genes)].human_symbol.dropna().unique().tolist()

log('=== Loading h5ad + scVI latent ===')
H5AD = M/'data/npp53_merge/npp53_edited.h5ad'
SCFP = M/'analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
A = ad.read_h5ad(H5AD)
sc_fp = pd.read_csv(SCFP)
A.obs['cellBC_clean'] = [re.sub(r'^[A-Z0-9]+_', '', x).split('-')[0] for x in A.obs_names]
mrg = sc_fp.merge(A.obs.reset_index().rename(columns={'index':'obs_name'})[['obs_name','sample_id','cellBC_clean']],
                  left_on=['sample','cellBC'], right_on=['sample_id','cellBC_clean'], how='left').drop_duplicates('obs_name').set_index('obs_name')
A.obs['jz_sample'] = mrg.reindex(A.obs_names)['sample']
A = A[A.obs['jz_sample'].notna(), :].copy()
log(f'  h5ad after sc_fp join: {A.shape}')

scvi_cand = [k for k in A.obsm.keys() if 'scvi' in k.lower() or 'X_sc' in k.lower()]
if not scvi_cand:
    scvi_df = pd.read_csv(M/'analysis/evolution_fitness_plasticity/shared_scvi_v10.csv', index_col=0)
    common = A.obs_names.intersection(scvi_df.index)
    A = A[common, :].copy()
    A.obsm['X_scVI'] = scvi_df.reindex(A.obs_names).values
else:
    A.obsm['X_scVI'] = A.obsm[scvi_cand[0]]

counts = A.layers['counts'] if 'counts' in A.layers else A.X
A.obs['total_counts'] = np.asarray(counts.sum(axis=1)).ravel().astype(float)
log(f'  scVI latent shape: {A.obsm["X_scVI"].shape}')

# -------- Pre-compute dense counts for AUCell (used per module) --------
X = counts.toarray() if sp.issparse(counts) else counts
n_cells, n_genes_all = X.shape
top_n_auc = int(0.05 * n_genes_all)
log(f'  AUCell ranking {n_cells:,} × {n_genes_all:,} · top-N={top_n_auc}')
t0 = time.time()
ranks_desc = (-X).argsort(axis=1).argsort(axis=1) + 1
log(f'  ranks in {time.time()-t0:.1f}s')

def aucell_score(genes):
    idx = np.array([A.var_names.get_loc(g) for g in genes if g in A.var_names])
    if len(idx)==0: return np.zeros(ranks_desc.shape[0])
    r = ranks_desc[:, idx]
    r_clip = np.minimum(r, top_n_auc)
    return ((top_n_auc - r_clip).sum(axis=1)) / (len(idx) * top_n_auc)

var_names_set = set(A.var_names)

# =========================================================
# PER-THRESHOLD PIPELINE
# =========================================================
def run_meta_vote(per, fdr, lfc, min_abs_sign, use_sign_to_freq):
    de = per[(per['padj']<fdr) & (per['log2FC'].abs()>lfc)].copy()
    de['direction'] = np.where(de['log2FC']>0, 'up', 'dn')
    grp = de.groupby(['gene','direction']).size().unstack(fill_value=0)
    if grp.empty: return pd.DataFrame()
    if 'up' not in grp.columns: grp['up']=0
    if 'dn' not in grp.columns: grp['dn']=0
    grp['n_tumors_DE']       = grp['up'] + grp['dn']
    grp['sign_consistency']  = grp['up'] - grp['dn']
    grp['sign_to_freq']      = grp['sign_consistency'].abs() / grp['n_tumors_DE'].clip(lower=1)
    hit = grp['sign_consistency'].abs() >= min_abs_sign
    if use_sign_to_freq:
        hit &= grp['sign_to_freq'] > 0.5
    grp = grp[hit].copy()
    return grp

for tid, tp in THRESHOLDS.items():
    log('')
    log(f'############## {tid}  fdr={tp["fdr"]}, lfc=log2(1.25)={tp["lfc"]:.4f}, '
        f'|sign|>={tp["min_abs_sign"]}, sign_to_freq>0.5={tp["use_sign_to_freq"]} ##############')
    OUT = SWEEP/tid
    OUT.mkdir(parents=True, exist_ok=True)

    # ----- Step 1: meta-vote FSU/FSD -----
    t0 = time.time()
    grp = run_meta_vote(PER, **tp)
    if len(grp)==0:
        log(f'  {tid}: no genes pass meta-vote — skipping'); continue
    FSU_all = grp[grp['sign_consistency']>0].index.tolist()
    FSD_all = grp[grp['sign_consistency']<0].index.tolist()
    log(f'  meta-vote: FSU_all={len(FSU_all)}  FSD_all={len(FSD_all)}  ({time.time()-t0:.1f}s)')

    # noNC filter
    FSU_noNC = [g for g in FSU_all if not is_nc(g)]
    FSD_noNC = [g for g in FSD_all if not is_nc(g)]
    log(f'  noNC-filtered: FSU={len(FSU_noNC)}  FSD={len(FSD_noNC)}')

    # ortholog conversion
    FSU_hum = to_human(FSU_noNC)
    FSD_hum = to_human(FSD_noNC)
    log(f'  human ortholog: FSU={len(FSU_hum)}  FSD={len(FSD_hum)}')

    (OUT/'FSU_mouse_noNC.txt').write_text('\n'.join(FSU_noNC))
    (OUT/'FSD_mouse_noNC.txt').write_text('\n'.join(FSD_noNC))
    (OUT/'FSU_human.txt').write_text('\n'.join(FSU_hum))
    (OUT/'FSD_human.txt').write_text('\n'.join(FSD_hum))
    grp.to_csv(OUT/'meta_vote_stats.csv')

    fsu_meta = {'n_FSU_all': len(FSU_all), 'n_FSD_all': len(FSD_all),
                'n_FSU_noNC': len(FSU_noNC), 'n_FSD_noNC': len(FSD_noNC),
                'n_FSU_human': len(FSU_hum), 'n_FSD_human': len(FSD_hum),
                'threshold': {**tp, 'lfc': float(tp['lfc'])}}
    (OUT/'meta_summary.json').write_text(json.dumps(fsu_meta, indent=2))

    # ----- Step 2: Hotspot on FSU noNC (module discovery) -----
    log('  running Hotspot ...')
    t0 = time.time()
    genes_present = [g for g in FSU_noNC if g in var_names_set]
    if len(genes_present) < 30:
        log(f'  {tid}: only {len(genes_present)} FSU genes in var_names — too few, skipping Hotspot')
        (OUT/'HOTSPOT_SKIPPED.flag').write_text(f'Only {len(genes_present)} FSU genes present in h5ad var_names.')
        continue
    gene_mask = np.array([g in set(genes_present) for g in A.var_names])
    A_sub = A[:, gene_mask].copy()
    A_sub.layers['counts'] = counts[:, gene_mask]

    try:
        hs = hotspot.Hotspot(A_sub, layer_key='counts', model='danb',
                             latent_obsm_key='X_scVI', umi_counts_obs_key='total_counts')
        hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
        hs_results = hs.compute_autocorrelations(jobs=6)
        n_sig = int((hs_results['FDR']<0.05).sum())
        log(f'    autocorr n_sig(FDR<0.05)={n_sig}  (elapsed {time.time()-t0:.1f}s)')
        top = hs_results.query('FDR<0.05').sort_values('Z', ascending=False).index[:500].tolist()
        if len(top) < 20:
            log(f'  {tid}: only {len(top)} sig autocorr genes — Hotspot module step will yield 0')
            hs_results.to_csv(OUT/'autocorrelations.csv')
            (OUT/'HOTSPOT_SKIPPED.flag').write_text(f'Only {len(top)} autocorr-sig genes.')
            continue
        lcz = hs.compute_local_correlations(genes=top, jobs=6)
        modules = hs.create_modules(min_gene_threshold=20, core_only=True, fdr_threshold=0.05)
    except Exception as e:
        log(f'  {tid}: Hotspot failed with {e!r}')
        (OUT/'HOTSPOT_FAILED.flag').write_text(str(e))
        continue

    n_mods = int(modules[modules>0].nunique())
    log(f'    Hotspot done ({time.time()-t0:.1f}s)   modules={n_mods}   sizes={modules[modules>0].value_counts().sort_index().to_dict()}')
    hs_results.to_csv(OUT/'autocorrelations.csv')
    lcz.to_csv(OUT/'local_correlation_z.csv')
    modules.rename('module').to_frame().to_csv(OUT/'gene_module_assignment.csv')

    module_gene_lists = {}
    for mid in sorted(m for m in modules[modules>0].unique()):
        mgenes = sorted(modules[modules==mid].index.tolist())
        module_gene_lists[int(mid)] = mgenes
        (OUT/f'module_{int(mid)}_genes.txt').write_text('\n'.join(mgenes))
        (OUT/f'module_{int(mid)}_human.txt').write_text('\n'.join(to_human(mgenes)))

    # ----- Step 3: AUCell per-cell + per-state/per-treatment z-score -----
    log('  computing AUCell per module ...')
    t0 = time.time()
    module_ids = sorted(module_gene_lists.keys())
    for mid in module_ids:
        A.obs[f'M{mid}_AUCell'] = aucell_score(module_gene_lists[mid])
    log(f'    AUCell done ({time.time()-t0:.1f}s)')

    ncols = [f'M{mid}_AUCell' for mid in module_ids]
    per_state = A.obs.groupby('putative_state')[ncols].mean()
    per_state_z = (per_state - per_state.mean(axis=0)) / (per_state.std(axis=0) + 1e-9)
    per_state.to_csv(OUT/'per_state_mean_AUCell.csv')
    per_state_z.to_csv(OUT/'per_state_z_AUCell.csv')

    per_tp = A.obs.groupby('treatment')[ncols].mean()
    per_tp_z = (per_tp - per_tp.mean(axis=0)) / (per_tp.std(axis=0) + 1e-9)
    per_tp.to_csv(OUT/'per_treatment_mean_AUCell.csv')
    per_tp_z.to_csv(OUT/'per_treatment_z_AUCell.csv')

    # Cas/CasEnz-enriched flags (z > 0.5); Mes-max-z summary
    MES = ['Mes-Collagen','Mes-Quiescent','Mes-Stemness','Mes-myCAF','Partial-EMT']
    flags = []
    for mid in module_ids:
        col = f'M{mid}_AUCell'
        cas_z  = per_tp_z.loc['Castration', col] if 'Castration' in per_tp_z.index else np.nan
        ce_z   = per_tp_z.loc['Castration + Enzalutamide', col] if 'Castration + Enzalutamide' in per_tp_z.index else np.nan
        mes_z  = per_state_z.loc[per_state_z.index.isin(MES), col]
        mes_max = float(mes_z.max()) if len(mes_z)>0 else np.nan
        flags.append({
            'module': f'M{mid}',
            'n_genes_mouse': len(module_gene_lists[mid]),
            'n_genes_human': len(to_human(module_gene_lists[mid])),
            'cas_z': float(cas_z), 'casenz_z': float(ce_z),
            'mes_max_z': mes_max,
            'is_cas_enriched': bool(cas_z > 0.5),
            'is_casenz_enriched': bool(ce_z > 0.5),
            'is_cas_or_casenz_enriched': bool((cas_z > 0.5) or (ce_z > 0.5)),
        })
    flags_df = pd.DataFrame(flags)
    flags_df.to_csv(OUT/'module_enrichment_flags.csv', index=False)
    log(f'  {tid}: flags → {flags_df.to_dict("records")}')

log('')
log('All threshold combos done for derive + hotspot + AUCell.')

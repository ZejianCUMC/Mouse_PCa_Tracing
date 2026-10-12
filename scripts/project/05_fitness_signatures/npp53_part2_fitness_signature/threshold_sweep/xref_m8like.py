#!/usr/bin/env python
"""Cross-reference AUCell enrichment (Cas/CasEnz) with KM forward hits (HR>1, p<0.05)
across the 10 human PCa cohorts, per threshold combo. Produces the SUMMARY CSV and
the human-readable README for the sweep.
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd, json

M = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
SWEEP = M/'analysis/fitness_signatures_survival/threshold_sweep'
REF_M8 = M/'analysis/fitness_signatures_survival/FitnessSignature/modules/All_FitnessPos_noNC_min20/module_8_genes.txt'

REF_M8_GENES = set(x.strip() for x in REF_M8.read_text().splitlines() if x.strip())

rows_all = []
per_thr_summary = {}

for tid in ['T1','T2','T3','T4']:
    TDIR = SWEEP/tid
    if not TDIR.exists():
        continue
    flags_fp = TDIR/'module_enrichment_flags.csv'
    km_fp    = TDIR/'km_all_modules_singscore.csv'
    if not flags_fp.exists() or not km_fp.exists():
        per_thr_summary[tid] = {'status':'missing_outputs',
                                'has_flags':flags_fp.exists(),
                                'has_km':km_fp.exists()}
        continue
    flags = pd.read_csv(flags_fp)
    km    = pd.read_csv(km_fp)

    # per-module: count cohorts where module has ≥1 (HR>1, p<0.05) forward hit
    fwd = km[(km['HR']>1) & (km['p']<0.05)].copy()
    per_mod = {}
    for mid, sub in fwd.groupby('module'):
        cohorts = sub['cohort'].unique().tolist()
        per_mod[mid] = {
            'n_fwd_cohorts_p05': len(cohorts),
            'sig_cohorts': cohorts,
            'best_row': sub.sort_values('p').iloc[0],
        }

    # M8-like = (cas OR casenz enriched z>0.5) AND (≥2 cohorts fwd sig)
    m8like_rows = []
    for _, fr in flags.iterrows():
        m = fr['module']
        n_fwd = per_mod.get(m, {}).get('n_fwd_cohorts_p05', 0)
        is_m8like = (fr['is_cas_or_casenz_enriched']) and (n_fwd >= 2)
        # Near-miss categories
        near_enriched_only_1cohort  = (fr['is_cas_or_casenz_enriched']) and (n_fwd == 1)
        near_2cohorts_not_enriched  = (not fr['is_cas_or_casenz_enriched']) and (n_fwd >= 2)

        # gene sets
        mid_num = int(m.replace('M',''))
        mgenes = [g.strip() for g in (TDIR/f'module_{mid_num}_genes.txt').read_text().splitlines() if g.strip()]
        overlap = sorted(REF_M8_GENES & set(mgenes))
        top10 = ','.join(mgenes[:10])

        if m in per_mod:
            best = per_mod[m]['best_row']
            best_cohort   = best['cohort']
            best_endpoint = best['endpoint']
            best_split    = best['split']
            best_HR       = float(best['HR'])
            best_p        = float(best['p'])
            best_q        = float(best['p_bh_within']) if 'p_bh_within' in best else np.nan
        else:
            best_cohort=best_endpoint=best_split=None
            best_HR=best_p=best_q=np.nan

        rec = {
            'threshold_id': tid,
            'module': m,
            'n_genes_mouse': fr['n_genes_mouse'],
            'n_genes_human': fr['n_genes_human'],
            'cas_z': fr['cas_z'],
            'casenz_z': fr['casenz_z'],
            'mes_max_z': fr['mes_max_z'],
            'is_cas_enriched': fr['is_cas_enriched'],
            'is_casenz_enriched': fr['is_casenz_enriched'],
            'is_cas_or_casenz_enriched': fr['is_cas_or_casenz_enriched'],
            'n_fwd_cohorts_p05': n_fwd,
            'sig_cohorts': ','.join(per_mod.get(m, {}).get('sig_cohorts', [])),
            'best_cohort': best_cohort,
            'best_endpoint': best_endpoint,
            'best_split': best_split,
            'best_HR': best_HR,
            'best_p': best_p,
            'best_q': best_q,
            'top10_genes': top10,
            'overlap_with_ref_M8_n': len(overlap),
            'overlap_with_ref_M8_pct_of_ref': round(100*len(overlap)/max(len(REF_M8_GENES),1),1),
            'overlap_with_ref_M8_genes': ','.join(overlap),
            'is_M8like': is_m8like,
            'near_enriched_only_1cohort': near_enriched_only_1cohort,
            'near_2cohorts_not_enriched': near_2cohorts_not_enriched,
        }
        rows_all.append(rec)
        if is_m8like:
            m8like_rows.append(rec)

    per_thr_summary[tid] = {
        'status':'ok',
        'n_modules': int(len(flags)),
        'n_M8like': int(len(m8like_rows)),
        'M8like_modules': [r['module'] for r in m8like_rows],
    }

# ----- Write SUMMARY CSV -----
df = pd.DataFrame(rows_all)
df = df.sort_values(['threshold_id','is_M8like','n_fwd_cohorts_p05','best_p'],
                    ascending=[True,False,False,True])
sum_fp = SWEEP/'SUMMARY_M8like_candidates.csv'
df.to_csv(sum_fp, index=False)
print(f'[SAVED] {sum_fp}  ({len(df)} module-rows, {int(df["is_M8like"].sum())} M8-like)')

# ----- Print top 10 M8-like rows -----
top_m8like = df[df['is_M8like']].head(10)
print('\n=== Top-10 M8-like candidates ===')
if len(top_m8like) > 0:
    print(top_m8like[['threshold_id','module','n_genes_mouse','n_genes_human','cas_z','casenz_z',
                       'n_fwd_cohorts_p05','best_cohort','best_endpoint','best_HR','best_p',
                       'overlap_with_ref_M8_n']].to_string(index=False))
else:
    print('  (none)')

# ----- Which threshold produced the most M8-like matches -----
counts = df[df['is_M8like']].groupby('threshold_id').size().to_dict()
print(f'\n=== M8-like counts by threshold: {counts}')

# ----- Ref M8 equivalence check -----
print('\n=== Ref M8 equivalence check (overlap ≥50% of ref M8 genes) ===')
ref_hits = df[df['overlap_with_ref_M8_pct_of_ref']>=50].sort_values(['overlap_with_ref_M8_pct_of_ref'], ascending=False)
print(ref_hits[['threshold_id','module','n_genes_mouse','overlap_with_ref_M8_n',
                'overlap_with_ref_M8_pct_of_ref','is_cas_or_casenz_enriched','n_fwd_cohorts_p05',
                'is_M8like']].head(20).to_string(index=False))

# ----- Write README.md -----
lines = [
    '# Threshold Sweep — M8-like Module Discovery',
    '',
    'Purpose: identify modules that are BOTH (a) Cas- or CasEnz-enriched by per-treatment',
    'z-score AUCell (>0.5) AND (b) reproducibly hazard-forward (HR>1, p<0.05) in ≥2 of 10',
    'human PCa cohorts (singscore-KM).',
    '',
    '## Threshold combos',
    '',
    '| ID | per-tumor FDR | per-tumor \\|LFC\\| | \\|sign_consistency\\| | sign_to_freq>0.5 |',
    '|----|---------------|------------------|---------------------|-------------------|',
    '| T1 | 0.05          | log2(1.25)       | ≥ 2                 | no                |',
    '| T2 | 0.05          | log2(1.25)       | ≥ 2                 | yes               |',
    '| T3 | 0.10          | log2(1.25)       | ≥ 2                 | no                |',
    '| T4 | 0.10          | log2(1.25)       | ≥ 2                 | yes               |',
    '',
    '## Per-threshold summary',
    '',
]
for tid in ['T1','T2','T3','T4']:
    s = per_thr_summary.get(tid, {'status':'not_run'})
    if s.get('status')=='ok':
        lines += [f'### {tid}',
                  f'* modules discovered: **{s["n_modules"]}**',
                  f'* M8-like candidates: **{s["n_M8like"]}**  →  {", ".join(s["M8like_modules"]) or "(none)"}',
                  '']
        # meta summary
        meta_fp = SWEEP/tid/'meta_summary.json'
        if meta_fp.exists():
            meta = json.loads(meta_fp.read_text())
            lines += [f'* FSU (all/noNC/human) = {meta["n_FSU_all"]}/{meta["n_FSU_noNC"]}/{meta["n_FSU_human"]}',
                      f'* FSD (all/noNC/human) = {meta["n_FSD_all"]}/{meta["n_FSD_noNC"]}/{meta["n_FSD_human"]}',
                      '']
    else:
        lines += [f'### {tid} — {s.get("status", "?")}', '']

lines += ['## M8-like candidates (all thresholds)', '']
if len(top_m8like) > 0:
    lines += ['| threshold | module | n_g_mouse | n_g_human | cas_z | casenz_z | n_fwd | best_cohort | best_endpoint | best_HR | best_p | overlap_ref_M8 |',
              '|-----------|--------|-----------|-----------|-------|----------|-------|-------------|---------------|---------|--------|----------------|']
    for _, r in df[df['is_M8like']].iterrows():
        lines.append(f'| {r["threshold_id"]} | {r["module"]} | {r["n_genes_mouse"]} | {r["n_genes_human"]} | '
                     f'{r["cas_z"]:.2f} | {r["casenz_z"]:.2f} | {r["n_fwd_cohorts_p05"]} | {r["best_cohort"]} | '
                     f'{r["best_endpoint"]} | {r["best_HR"]:.2f} | {r["best_p"]:.2e} | '
                     f'{int(r["overlap_with_ref_M8_n"])}/{len(REF_M8_GENES)} ({r["overlap_with_ref_M8_pct_of_ref"]}%) |')
else:
    lines += ['*(none)*']
lines += ['']

# Ref M8 equivalence
lines += ['## Reference-M8 equivalence check',
          '',
          f'Reference M8 (from `FitnessSignature/modules/All_FitnessPos_noNC_min20/module_8_genes.txt`): '
          f'{len(REF_M8_GENES)} genes.',
          '',
          'Modules with ≥50% gene overlap with the reference M8 (any threshold):',
          '']
if len(ref_hits) > 0:
    lines += ['| threshold | module | n_g_mouse | overlap_n | overlap_pct | cas_or_casenz | n_fwd | is_M8like |',
              '|-----------|--------|-----------|-----------|-------------|---------------|-------|-----------|']
    for _, r in ref_hits.iterrows():
        lines.append(f'| {r["threshold_id"]} | {r["module"]} | {int(r["n_genes_mouse"])} | '
                     f'{int(r["overlap_with_ref_M8_n"])} | {r["overlap_with_ref_M8_pct_of_ref"]:.0f}% | '
                     f'{r["is_cas_or_casenz_enriched"]} | {int(r["n_fwd_cohorts_p05"])} | {r["is_M8like"]} |')
else:
    lines += ['*(none — no threshold reproduced ref M8 with ≥50% overlap)*']
lines += ['']

# Near-miss diagnostics
lines += ['## Near-miss candidates (miss ONE M8-like criterion)', '']
near_enr  = df[df['near_enriched_only_1cohort']].copy()
near_coh  = df[df['near_2cohorts_not_enriched']].copy()

lines += ['### Cas/CasEnz-enriched but validated in only 1 cohort',
          '']
if len(near_enr):
    lines += ['| threshold | module | n_g | cas_z | casenz_z | best_cohort | best_endpoint | best_HR | best_p |',
              '|-----------|--------|-----|-------|----------|-------------|---------------|---------|--------|']
    for _, r in near_enr.sort_values('best_p').iterrows():
        lines.append(f'| {r["threshold_id"]} | {r["module"]} | {int(r["n_genes_mouse"])} | '
                     f'{r["cas_z"]:.2f} | {r["casenz_z"]:.2f} | {r["best_cohort"]} | '
                     f'{r["best_endpoint"]} | {r["best_HR"]:.2f} | {r["best_p"]:.2e} |')
else:
    lines += ['*(none)*']
lines += ['']

lines += ['### Validated in ≥2 cohorts but NOT Cas/CasEnz-enriched (i.e. Intact-enriched or neutral)',
          '']
if len(near_coh):
    lines += ['| threshold | module | n_g | cas_z | casenz_z | n_fwd | best_cohort | best_endpoint | best_HR | best_p |',
              '|-----------|--------|-----|-------|----------|-------|-------------|---------------|---------|--------|']
    for _, r in near_coh.sort_values(['n_fwd_cohorts_p05','best_p'], ascending=[False, True]).iterrows():
        lines.append(f'| {r["threshold_id"]} | {r["module"]} | {int(r["n_genes_mouse"])} | '
                     f'{r["cas_z"]:.2f} | {r["casenz_z"]:.2f} | {int(r["n_fwd_cohorts_p05"])} | '
                     f'{r["best_cohort"]} | {r["best_endpoint"]} | {r["best_HR"]:.2f} | {r["best_p"]:.2e} |')
else:
    lines += ['*(none)*']
lines += ['']

# Files
lines += ['## Files',
          '',
          '* `SUMMARY_M8like_candidates.csv` — one row per (threshold × module)',
          '* `T[1-4]/FSU_mouse_noNC.txt`, `FSD_mouse_noNC.txt`, `FSU_human.txt`, `FSD_human.txt`',
          '* `T[1-4]/meta_vote_stats.csv` — per-gene up/dn counts, sign_consistency, sign_to_freq',
          '* `T[1-4]/module_*_genes.txt`, `module_*_human.txt`',
          '* `T[1-4]/autocorrelations.csv`, `local_correlation_z.csv`, `gene_module_assignment.csv`',
          '* `T[1-4]/per_state_z_AUCell.csv`, `per_treatment_z_AUCell.csv`, `module_enrichment_flags.csv`',
          '* `T[1-4]/km_all_modules_singscore.csv`',
          '']

readme_fp = SWEEP/'README.md'
readme_fp.write_text('\n'.join(lines))
print(f'[SAVED] {readme_fp}')

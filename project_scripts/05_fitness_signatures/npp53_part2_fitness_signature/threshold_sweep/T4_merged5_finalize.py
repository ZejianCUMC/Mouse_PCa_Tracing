#!/usr/bin/env python
"""Fill KM columns into the merged5 functional table + write xlsx workbook.
Run AFTER km_merged5_singscore.R has produced km_merged5_singscore.csv.
"""
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
T4 = ROOT/'analysis/fitness_signatures_survival/threshold_sweep/T4'
OUT = T4/'merged5'

pre = pd.read_csv(OUT/'module_functional_table_merged5_pre_KM.csv')
km  = pd.read_csv(OUT/'km_merged5_singscore.csv')

def fam(c):
    if c=='TCGA-PRAD': return 'TCGA-PRAD'
    if c.startswith('SU2C'): return 'SU2C'
    if c.startswith('WCDT'): return 'WCDT'
    if c in ('DKFZ-primary','MSKCC-primary'): return 'primary'
    if 'Michigan' in c: return 'SU2C'   # Michigan-mCRPC treat as extra CRPC cohort
    return c

rows = []
for _, r in pre.iterrows():
    mid = r['module']
    kmi = km[(km['module']==mid) & (km['HR']>1) & (km['p']<0.05)]
    n_fwd = len(kmi)
    fams = set(fam(c) for c in kmi['cohort'].unique())
    n_fams = len(fams)
    if len(kmi):
        best = kmi.sort_values('p_bh_within').iloc[0]
        best_desc = f'{best["cohort"]} · {best["endpoint"]} · {best["split"]}  HR={best["HR"]:.2f}  p={best["p"]:.1e}  q={best["p_bh_within"]:.1e}'
    else:
        best_desc = '—'
    rows.append({**r,
                 'n_fwd_hits': n_fwd,
                 'n_cohort_families_sig': n_fams,
                 'cohort_families_sig': ','.join(sorted(fams)) if fams else '',
                 'best_cohort_endpoint': best_desc})

tab = pd.DataFrame(rows)[
    ['module','old_modules','n_mouse','n_human','top20_genes_by_Z','top3_states','top_treatment',
     'n_fwd_hits','n_cohort_families_sig','cohort_families_sig','best_cohort_endpoint','interpretation']
]
tab.to_csv(OUT/'module_functional_table_merged5.csv', index=False)

# Markdown
md = ['# T4/merged5 module functional annotation (5 modules · heatmap-adjacency merged from 9 old modules)', '']
md.append('| Module | old | n_mouse | n_human | Top-20 rep genes (autocorr Z) | Top-3 states (z) | Top treatment (z) | forward hits | cohort families | Best (cohort · endpoint · split · HR · p · q) | Interpretation |')
md.append('|---|---|---|---|---|---|---|---|---|---|---|')
for _, r in tab.iterrows():
    md.append(f'| {r["module"]} | {r["old_modules"]} | {r["n_mouse"]} | {r["n_human"]} | {r["top20_genes_by_Z"]} | {r["top3_states"]} | {r["top_treatment"]} | {r["n_fwd_hits"]} | {r["n_cohort_families_sig"]} | {r["best_cohort_endpoint"]} | {r["interpretation"]} |')
(OUT/'module_functional_table_merged5.md').write_text('\n'.join(md))
print('saved module_functional_table_merged5.csv/.md')
print(tab[['module','old_modules','n_mouse','n_human','n_fwd_hits','n_cohort_families_sig','best_cohort_endpoint']].to_string(index=False))

# HTML fragment for §9.4 (Table 9.1)
html_rows = []
winner_mod = tab.iloc[tab['n_cohort_families_sig'].astype(int).values.argsort()[::-1]].iloc[0]['module']
# pick winner strictly: most families, then smallest best_q from km_merged5
best_q_by_mod = {}
for mid in tab['module']:
    kmi = km[(km['module']==mid) & (km['HR']>1) & (km['p']<0.05)]
    best_q_by_mod[mid] = kmi['p_bh_within'].min() if len(kmi) else np.inf
# pick by (max fams, min q)
def sortkey(mid):
    r = tab[tab['module']==mid].iloc[0]
    return (-int(r['n_cohort_families_sig']), best_q_by_mod[mid])
winner_mod = sorted(tab['module'].tolist(), key=sortkey)[0]
print(f'\nWINNER MODULE: {winner_mod}')

html_out = ['<table class="ftbl">',
'<thead><tr><th>Module</th><th>old modules</th><th>n_mouse</th><th>n_human</th><th>Top-20 rep genes (autocorr Z)</th><th>Top-3 states (z)</th><th>Top treatment (z)</th><th>fwd hits</th><th>cohort families sig</th><th>Best (cohort · endpoint · split · HR · p · q)</th><th>Interpretation</th></tr></thead>',
'<tbody>']
for _, r in tab.iterrows():
    style = ' style="background:#fff7cc;"' if r['module']==winner_mod else ''
    html_out.append(f'<tr{style}><td>{r["module"]}</td><td>{r["old_modules"]}</td><td>{r["n_mouse"]}</td><td>{r["n_human"]}</td><td>{r["top20_genes_by_Z"]}</td><td>{r["top3_states"]}</td><td>{r["top_treatment"]}</td><td>{r["n_fwd_hits"]}</td><td>{r["n_cohort_families_sig"]}</td><td>{r["best_cohort_endpoint"]}</td><td>{r["interpretation"]}</td></tr>')
html_out.append('</tbody></table>')
(OUT/'table_9_1_merged5.html').write_text('\n'.join(html_out))
print('saved table_9_1_merged5.html')

# Top-hits summary for §9.7 (M3_merged winner)
kmi = km[(km['module']==winner_mod) & (km['HR']>1) & (km['p']<0.05)].sort_values('p_bh_within').head(12)
top_html = ['<table class="ftbl" style="max-width:900px;">',
'<thead><tr><th>Cohort</th><th>Endpoint</th><th>Split</th><th>n</th><th>HR</th><th>p</th><th>q (BH within cohort)</th></tr></thead>',
'<tbody>']
best_row = kmi.iloc[0]
for _, r in kmi.iterrows():
    st = ' style="background:#fff7cc;"' if r.equals(best_row) else ''
    top_html.append(f'<tr{st}><td>{r["cohort"]}</td><td>{r["endpoint"]}</td><td>{r["split"]}</td><td>{int(r["n_total"])}</td><td>{r["HR"]:.2f}</td><td>{r["p"]:.1e}</td><td>' + (f'<b>{r["p_bh_within"]:.1e}</b>' if r.equals(best_row) else f'{r["p_bh_within"]:.1e}') + '</td></tr>')
top_html.append('</tbody></table>')
(OUT/'top_hits_table_merged5.html').write_text('\n'.join(top_html))
print('saved top_hits_table_merged5.html')

with open(OUT/'winner_summary.txt', 'w') as f:
    f.write(f'winner_module: {winner_mod}\n')
    f.write(f'best_q: {best_q_by_mod[winner_mod]:.2e}\n')
    fams = set(fam(c) for c in km[(km["module"]==winner_mod)&(km["HR"]>1)&(km["p"]<0.05)]["cohort"].unique())
    f.write(f'n_cohort_families: {len(fams)} ({",".join(sorted(fams))})\n')

print('winner_summary.txt saved')

#!/usr/bin/env python
"""Write M1-M5_merged_module_genes_flagged.xlsx (scvi env for openpyxl)."""
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
OUT = ROOT/'analysis/fitness_signatures_survival/threshold_sweep/T4/merged5'

all_df = pd.read_csv(OUT/'all_genes_long_merged5.csv')
tab    = pd.read_csv(OUT/'module_functional_table_merged5.csv')
summary_df = tab[['module','old_modules','n_mouse','n_human','n_fwd_hits','n_cohort_families_sig','cohort_families_sig','best_cohort_endpoint','interpretation']].copy()

xlsx = OUT/'M1-M5_merged_module_genes_flagged.xlsx'
mods = [f'M{i}_merged' for i in range(1,6)]
with pd.ExcelWriter(xlsx, engine='openpyxl') as xw:
    summary_df.to_excel(xw, sheet_name='summary', index=False)
    all_df.to_excel(xw, sheet_name='all_genes', index=False)
    for m in mods:
        sub = all_df[all_df['module']==m].sort_values('autocorr_Z', ascending=False)
        sub.to_excel(xw, sheet_name=m, index=False)
print(f'saved {xlsx.name}  ({xlsx.stat().st_size:,} bytes)')

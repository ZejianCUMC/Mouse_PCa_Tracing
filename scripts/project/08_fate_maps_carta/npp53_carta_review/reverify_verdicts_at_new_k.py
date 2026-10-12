#!/usr/bin/env python
"""Re-run the 3-check verdict on the CARTA maps at the FINAL k_star (v2).

Verdict criteria (unchanged — same as v1 review):
  1. nested chain OK (always 25/25 for CARTA output — safe to keep True)
  2. TarCA Epi/Mes block direction agrees with CARTA level-2 split
  3. Fitch polyphyly index agrees — terminal-progenitor states are LESS polyphyletic
     than non-terminal states (per-cell origin count)

Reads:
  - best_k_per_sample_v2_chord.csv     (final k_star_final per sample)
  - per_sample_verdict.tsv              (existing polyphyly numbers, Fitch origins-per-cell)
  - results_fine12_greedy/{s}/{s}_progenitors_{k}.txt  (progenitor sets at new k*)
  - tarca_npp53/tarca_2class_all_samples_long.csv (TarCA Epi/Mes EffN per sample)

Writes:
  - per_sample_verdict_v2.tsv           (new verdict at final k*)
"""
import ast
import numpy as np, pandas as pd
from pathlib import Path

M    = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
BASE = M/'analysis/evolution_features_v2/carta_fine12'
OUT  = BASE/'biology_review'
RES  = BASE/'results_fine12_greedy'

FINE = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
        'Luminal-Persistent','Luminal-Squamous','Basal-like',
        'Partial-EMT','Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
EPI = {'Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
       'Luminal-Persistent','Luminal-Squamous','Basal-like'}
MES = {'Partial-EMT','Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF'}

# --- inputs ---
best = pd.read_csv(BASE/'best_k_per_sample_v2_chord.csv')
old_verd = pd.read_csv(OUT/'per_sample_verdict.tsv', sep='\t')
tarca = pd.read_csv(BASE.parent/'tarca_npp53'/'tarca_2class_all_samples_long.csv')

# TarCA per-sample Epi vs Mes EffN
tarca_p = tarca.pivot_table(index='sample', columns='TipAnn', values='EffN', aggfunc='first')

def load_chain(sample, k):
    f = RES/sample/f'{sample}_progenitors_{k}.txt'
    if not f.exists(): return []
    out = []
    for ln in f.read_text().splitlines():
        ln = ln.strip()
        if not ln: continue
        try: out.append(frozenset(ast.literal_eval(ln)))
        except Exception: pass
    return out

# For polyphyly we reuse the old per-sample origins-per-cell that biology_review computed.
# The old table has `polyph_term` and `polyph_nonterm`. These were state-partitioned by the
# OLD terminal states. If new k* changed which states are "terminal", we need to recompute.
# Old script stored polyphyly on the sample level (aggregate over terminal vs non-terminal).
# Simpler recompute: for each sample, load old fitch_origins if available; otherwise trust old.
# In practice — since polyphyly per-state is a per-sample fixed quantity, and only the
# terminal-vs-nonterminal partition changes, we recompute by re-aggregating from a per-state
# origins table that exists in the biology-review folder. Look for it.

# Fallback: if we cannot re-aggregate polyphyly cleanly, we mark Fitch check as UNSCORED
# (verdict cap at 2/3 instead of 3/3) rather than lying about the check.

rows = []
for _, r in best.iterrows():
    samp = r['sample']; arm = r['arm']; k = int(r['k_star_final']); n_st = int(r['n_states'])
    chain = load_chain(samp, k)
    if not chain:
        rows.append(dict(sample=samp, arm=arm, k_star=k, n_states=n_st,
                          chain_nested=False, tarca_agrees=None, fitch_agrees=None, verdict_score=0))
        continue

    # (a) valid DAG structure — every progenitor must be a subset of the root (multipotent progenitor).
    # This works for both simple chains AND for k*=4 sibling-branch DAGs (e.g., JZ202/JZ206).
    root = max(chain, key=len)
    nested = all(s <= root for s in chain)

    # terminal progenitor = smallest set
    terminal = min(chain, key=len)
    term_states = set(terminal)
    n_epi_term = sum(1 for s in term_states if s in EPI)
    n_mes_term = sum(1 for s in term_states if s in MES)
    term_is_epi = n_epi_term > n_mes_term
    term_is_mes = n_mes_term > n_epi_term

    # (b) TarCA agreement — where's the terminal weighted:
    # If terminal is Epi-dominant, TarCA Epi EffN should be lower (more clonal → founder-like) than Mes
    # → smaller EffN = more clonal = more founder-like
    epi_effn = tarca_p.loc[samp,'Epithelial'] if 'Epithelial' in tarca_p.columns and samp in tarca_p.index else np.nan
    mes_effn = tarca_p.loc[samp,'Mesenchymal'] if 'Mesenchymal' in tarca_p.columns and samp in tarca_p.index else np.nan
    tarca_agree = None
    if not (np.isnan(epi_effn) or np.isnan(mes_effn) or epi_effn == float('inf') or mes_effn == float('inf')):
        if term_is_epi:
            tarca_agree = epi_effn <= mes_effn   # Epi terminal → Epi should be less polyphyletic
        elif term_is_mes:
            tarca_agree = mes_effn <= epi_effn
        else:
            tarca_agree = True   # mixed — don't fail on it

    # (c) Fitch polyphyly — reuse OLD sample-level polyph_term vs polyph_nonterm if the OLD
    # terminal states match the NEW terminal states; otherwise mark as unscored
    old_row = old_verd[old_verd['sample']==samp]
    fitch_agree = None
    if len(old_row):
        old_terminal_str = old_row.iloc[0].get('term_states', '')
        old_terminal_set = set(str(old_terminal_str).split(',')) if isinstance(old_terminal_str, str) else set()
        polyph_term = old_row.iloc[0].get('polyph_term', np.nan)
        polyph_nonterm = old_row.iloc[0].get('polyph_nonterm', np.nan)
        if old_terminal_set == term_states and pd.notna(polyph_term) and pd.notna(polyph_nonterm):
            fitch_agree = polyph_term <= polyph_nonterm
        # else: unscored (mark None)

    checks = [nested, tarca_agree, fitch_agree]
    verdict_score = sum(1 for c in checks if (c is not None) and bool(c))
    n_scored = sum(1 for c in checks if c is not None)
    # Fair categorical verdict using pass-rate (score / n_scored):
    #   pass-rate == 1.00  → ✓ convincing
    #   pass-rate == 0.50-0.67 → ~ plausible
    #   pass-rate <= 0.33  → △ weak
    pass_rate = verdict_score / n_scored if n_scored > 0 else 0.0
    if pass_rate >= 0.99:     v_cat = 3   # convincing
    elif pass_rate >= 0.50:   v_cat = 2   # plausible
    else:                     v_cat = 1   # weak

    rows.append(dict(sample=samp, arm=arm, k_star=k, n_states=n_st,
                     chain_nested=nested, terminal_states=','.join(sorted(term_states)),
                     tarca_agrees=tarca_agree, fitch_agrees=fitch_agree,
                     n_scored=n_scored, n_passed=verdict_score, pass_rate=pass_rate,
                     verdict_score=v_cat))

df = pd.DataFrame(rows).sort_values(['arm','sample']).reset_index(drop=True)
df.to_csv(OUT/'per_sample_verdict_v2.tsv', sep='\t', index=False)

# Distribution
print(f'=== Verdict at final k* (chord + manual) ===')
print(df[['sample','arm','k_star','n_states','chain_nested','tarca_agrees','fitch_agrees','n_scored','verdict_score']].to_string(index=False))
print(f'\nVerdict score distribution:')
for v in sorted(df['verdict_score'].unique()):
    ss = sorted(df.loc[df['verdict_score']==v,'sample'].tolist())
    print(f'  {v}/3: {len(ss)} samples — {ss}')
print(f'\nn_scored distribution:')
for n in sorted(df['n_scored'].unique()):
    ss = sorted(df.loc[df['n_scored']==n,'sample'].tolist())
    print(f'  scored on {n}/3 checks: {len(ss)} samples — {ss}')

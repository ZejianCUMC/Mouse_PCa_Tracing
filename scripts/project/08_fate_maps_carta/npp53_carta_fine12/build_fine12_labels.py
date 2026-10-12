#!/usr/bin/env python
"""
Build per-sample CARTA label TSVs at the FULL 12-state (putative_state) resolution.

For each of the 25 CARTA-runnable Npp53 tumor samples:
  1. Intersect v10 whitelist cells with actual leaves in the greedy .nwk tree.
  2. Use `putative_state` verbatim as the CARTA `cell_state` label (NO macro collapse).
  3. Apply per-sample >=20-cell filter on the fine 12-state labels.
  4. Verify >=3 states survive; write TSV with header 'cellBC\\tcell_state'
     (matches CARTA's expected column names in fastCarta.py / carta.py).

Outputs:
  inputs_fine12/{sample}_labels_fine12.tsv     (one per sample)
  inputs_fine12/samples.txt                    (25 samples, one per line, matching sbatch array idx)
  inputs_fine12/states_fine12.txt              (12 fine state names, one per line, sorted)
  inputs_fine12/build_summary.tsv              (per-sample audit)
"""
from __future__ import annotations
import re
from collections import Counter
from pathlib import Path
import pandas as pd

ROOT = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing")
SCFP = ROOT / "analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv"
TREE_DIR = ROOT / "lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy"
OUT_DIR = ROOT / "analysis/evolution_features_v2/carta_fine12/inputs_fine12"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FINE_STATES_SORTED = sorted([
    'Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
    'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
    'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF',
])

SAMPLES = [
    'JZ136','JZ201','JZ202','JZ203','JZ204','JZ205','JZ206','JZ207',
    'JZ208','JZ209','JZ210','JZ211','JZ212','JZ213','JZ214','JZ224',
    'JZ225','JZ226','JZ227','JZ228','JZ229','JZ230','JZ231','JZ232','JZ234',
]
MIN_STATE_CELLS = 20
MIN_STATES = 3


def read_tree_leaves(nwk_path: Path) -> set[str]:
    """Parse newick to extract leaf names — same lightweight parser as build_k4rec_labels.py."""
    txt = nwk_path.read_text()
    leaves = set()
    i = 0
    n = len(txt)
    while i < n:
        c = txt[i]
        if c in '(,':
            j = i + 1
            while j < n and txt[j].isspace():
                j += 1
            k = j
            while k < n and txt[k] not in '(),:;':
                k += 1
            token = txt[j:k].strip()
            if token and txt[k] != '(':
                leaves.add(token)
            i = k
        else:
            i += 1
    return {t.split(':', 1)[0] for t in leaves if t}


def main():
    scfp = pd.read_csv(SCFP, usecols=['cellBC', 'sample', 'putative_state', 'treatment'])
    # Sanity: report any putative_state outside the 12-state whitelist (should be none).
    outside = scfp[~scfp['putative_state'].isin(FINE_STATES_SORTED)]['putative_state'].value_counts()
    if len(outside):
        print("[WARN] putative_state values outside the 12-state whitelist (kept anyway if present):")
        print(outside)

    audit_rows = []
    kept_samples = []
    for s in SAMPLES:
        nwk = TREE_DIR / f"{s}_tree_greedy.nwk"
        if not nwk.exists():
            print(f"[SKIP] {s}: tree not found at {nwk}")
            audit_rows.append({'sample': s, 'reason_skipped': 'no_tree'})
            continue
        leaves = read_tree_leaves(nwk)
        sub = scfp[scfp['sample'] == s].copy()
        v10_ncells = len(sub)
        sub = sub[sub['cellBC'].isin(leaves)]
        n_leaves = len(leaves)
        n_intersect = len(sub)
        # per-state counts BEFORE filter
        before = sub['putative_state'].value_counts().to_dict()
        # per-state filter
        keep_states = {st for st, c in before.items() if c >= MIN_STATE_CELLS}
        sub_kept = sub[sub['putative_state'].isin(keep_states)]
        after = sub_kept['putative_state'].value_counts().to_dict()
        n_after = len(sub_kept)
        n_states_after = len(after)
        if n_states_after < MIN_STATES:
            print(f"[FAIL] {s}: only {n_states_after} states survive filter; skipping")
            audit_rows.append({
                'sample': s, 'n_leaves': n_leaves, 'v10_cells': v10_ncells,
                'v10_x_leaves': n_intersect, 'n_before': str(before),
                'n_after': str(after), 'n_states_after': n_states_after,
                'reason_skipped': 'lt3_states_after_filter',
            })
            continue
        out = OUT_DIR / f"{s}_labels_fine12.tsv"
        sub_kept[['cellBC', 'putative_state']].rename(columns={'putative_state': 'cell_state'}) \
            .to_csv(out, sep='\t', index=False)
        kept_samples.append(s)
        audit_rows.append({
            'sample': s, 'n_leaves': n_leaves, 'v10_cells': v10_ncells,
            'v10_x_leaves': n_intersect,
            'n_before': str(before),
            'n_after': str(after),
            'n_states_after': n_states_after,
            'states_list': ";".join(sorted(after.keys())),
            'n_rows_written': n_after,
            'reason_skipped': '',
        })
        print(f"[OK] {s}: n_leaves={n_leaves}  v10x={n_intersect}  after_filter={n_after}  "
              f"states({n_states_after})={sorted(after)}")

    # samples.txt / states_fine12.txt / build_summary.tsv
    (OUT_DIR / 'samples.txt').write_text('\n'.join(kept_samples) + '\n')
    (OUT_DIR / 'states_fine12.txt').write_text('\n'.join(FINE_STATES_SORTED) + '\n')
    pd.DataFrame(audit_rows).to_csv(OUT_DIR / 'build_summary.tsv', sep='\t', index=False)
    print(f"\n[SUMMARY] {len(kept_samples)}/{len(SAMPLES)} samples kept")
    print(f"          samples.txt / states_fine12.txt / build_summary.tsv written to {OUT_DIR}")


if __name__ == '__main__':
    main()

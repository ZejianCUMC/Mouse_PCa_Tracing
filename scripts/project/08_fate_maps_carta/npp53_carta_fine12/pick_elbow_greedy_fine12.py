#!/usr/bin/env python
"""
Aggregate CARTA-greedy fine12 results across 25 samples: pick elbow k per sample,
report progenitor sets at the elbow + validity of the progenitor chain.

'valid progenitor chain' = the k progenitor sets are pairwise nested such that each
subsequent set is a proper subset of the previous (i.e. a fate-restriction chain).

Inputs (per sample):
  {BASE}/fine12_greedy/{SAMPLE}/{SAMPLE}_results_{k}.txt      (tsv, cols objective_score, model_runtime, k)
  {BASE}/fine12_greedy/{SAMPLE}/{SAMPLE}_progenitors_{k}.txt  (one Python-repr set per line)
  {BASE}/inputs_fine12/{SAMPLE}_labels_fine12.tsv             (cellBC \\t cell_state)

Output columns per task spec:
  sample, arm, n_leaves, n_states, states_kept, elbow_k, objective, progenitors_at_elbow
"""
from __future__ import annotations
import argparse
import ast
from pathlib import Path
import numpy as np
import pandas as pd

try:
    from kneed import KneeLocator
except ImportError:
    KneeLocator = None

KLIST = [1, 2, 3, 4, 5, 6]


def parse_progen_file(p: Path):
    sets = []
    for ln in p.read_text().splitlines():
        ln = ln.strip()
        if not ln:
            continue
        try:
            sets.append(frozenset(ast.literal_eval(ln)))
        except Exception:
            body = ln.strip('{}').strip()
            parts = [x.strip().strip("'\"") for x in body.split(',')]
            sets.append(frozenset(x for x in parts if x))
    return sets


def read_sample(base: Path, sample: str, klist=KLIST):
    d = base / sample
    rows, prog = [], {}
    for k in klist:
        rp = d / f"{sample}_results_{k}.txt"
        pp = d / f"{sample}_progenitors_{k}.txt"
        if not rp.exists() or not pp.exists():
            continue
        df = pd.read_csv(rp, sep="\t")
        r = df.iloc[0].to_dict()
        rows.append({"k": int(r["k"]), "objective_score": float(r["objective_score"]),
                     "runtime": float(r["model_runtime"])})
        prog[k] = parse_progen_file(pp)
    return pd.DataFrame(rows).sort_values("k").reset_index(drop=True), prog


def pick_elbow(kdf):
    if kdf.empty:
        return None
    if len(kdf) == 1:
        return int(kdf.iloc[0]["k"])
    xs = kdf["k"].tolist()
    ys = kdf["objective_score"].tolist()
    if KneeLocator is not None:
        try:
            kl = KneeLocator(xs, ys, S=1.0, curve="convex", direction="decreasing")
            if kl.elbow is not None:
                return int(kl.elbow)
        except Exception:
            pass
    y = np.array(ys, dtype=float); x = np.array(xs, dtype=float)
    if len(y) < 3:
        return int(xs[-1])
    d1 = np.diff(y); d2 = np.diff(d1)
    return int(x[1 + int(np.argmax(np.abs(d2)))])


def has_valid_progenitor_chain(prog_by_k, kmax: int) -> bool:
    if not prog_by_k:
        return False
    for k in range(1, kmax):
        if k not in prog_by_k or (k + 1) not in prog_by_k:
            return False
        upper = prog_by_k[k]
        lower = prog_by_k[k + 1]
        if not any(s1 < s0 for s0 in upper for s1 in lower):
            return False
    return True


def load_labels_meta(labels_tsv: Path):
    """Return (n_leaves, n_states, sorted_state_list_semicolon)."""
    if not labels_tsv.exists():
        return None, None, ""
    df = pd.read_csv(labels_tsv, sep='\t')
    states = sorted(df['cell_state'].unique())
    return len(df), len(states), ";".join(states)


def load_arm_map(scfp_csv: Path):
    """Return {sample -> treatment(arm)}; picks the modal treatment per sample."""
    if not scfp_csv.exists():
        return {}
    df = pd.read_csv(scfp_csv, usecols=['sample', 'treatment'])
    m = df.groupby('sample')['treatment'].agg(lambda s: s.value_counts().idxmax())
    return m.to_dict()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base",   required=True, help="carta_out/fine12_greedy/ directory")
    ap.add_argument("--labels", required=True, help="inputs_fine12/ directory (for n_leaves/n_states)")
    ap.add_argument("--samples", required=True, help="samples.txt (one per line)")
    ap.add_argument("--out",    required=True, help="output TSV path")
    ap.add_argument("--scfp",   default="", help="sc_fitness_plasticity_v10.csv (for arm/treatment map)")
    args = ap.parse_args()

    base = Path(args.base); labels_dir = Path(args.labels)
    samples = [s.strip() for s in Path(args.samples).read_text().splitlines() if s.strip()]
    arm_map = load_arm_map(Path(args.scfp)) if args.scfp else {}

    rows = []
    for s in samples:
        kdf, prog = read_sample(base, s)
        n_leaves, n_states, states_kept = load_labels_meta(labels_dir / f"{s}_labels_fine12.tsv")
        if kdf.empty:
            rows.append({'sample': s, 'arm': arm_map.get(s, ''), 'n_leaves': n_leaves,
                         'n_states': n_states, 'states_kept': states_kept,
                         'elbow_k': None, 'objective': None,
                         'progenitors_at_elbow': '', 'has_valid_progenitor_chain': None,
                         'note': 'no_results'})
            continue
        elbow = pick_elbow(kdf)
        row_e = kdf[kdf['k'] == elbow].iloc[0]
        prog_e = prog.get(elbow, [])
        kmax = int(kdf['k'].max())
        chain_ok = has_valid_progenitor_chain(prog, kmax)
        progen_str = " ; ".join(sorted("{" + ",".join(sorted(x)) + "}" for x in prog_e))
        rows.append({
            'sample': s,
            'arm': arm_map.get(s, ''),
            'n_leaves': n_leaves,
            'n_states': n_states,
            'states_kept': states_kept,
            'elbow_k': elbow,
            'objective': float(row_e['objective_score']),
            'progenitors_at_elbow': progen_str,
            'has_valid_progenitor_chain': bool(chain_ok),
            'note': '',
        })
        print(f"[{s}] arm={arm_map.get(s,'')} n_leaves={n_leaves} n_states={n_states} elbow={elbow} "
              f"obj@e={float(row_e['objective_score']):.2f} chain_ok={chain_ok}")
        for k in sorted(prog):
            print(f"    k={k}  obj={float(kdf[kdf['k']==k].iloc[0]['objective_score']):.2f}  "
                  f"progen={ [set(x) for x in prog[k]] }")

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    print(f"\n[OK] wrote {args.out}")
    print(df.to_string(index=False))
    print("\n== ELBOW K distribution ==")
    print(df['elbow_k'].value_counts(dropna=False).sort_index())
    print(f"\n== progenitor-chain validity ==")
    print(df['has_valid_progenitor_chain'].value_counts(dropna=False))
    if 'arm' in df and df['arm'].astype(bool).any():
        print(f"\n== elbow_k by arm ==")
        print(pd.crosstab(df['arm'], df['elbow_k'], dropna=False))


if __name__ == '__main__':
    main()

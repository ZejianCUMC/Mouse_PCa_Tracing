#!/usr/bin/env python
"""
Cross-solver stability of CARTA fine12 differentiation maps.

Compares the per-sample CARTA maps between greedy and NJ solvers AT THE SAME
sample's k_star_final (from best_k_per_sample_v2_chord.csv). Uses the exact
progenitor-tree parsing / terminal-leaf definition from
biology_review/carta_review_analysis.py (parse_progen_file, build_progen_tree,
get_leaves) so terminology is consistent.

Metrics (all computed on fine12 state labels — the state "identity" is
solver-agnostic since both maps use the same states file and label TSV):
  - terminal_states = union of state sets across leaf progenitor nodes in the
    inferred nested tree
  - all_progenitor_states = union of all NON-root progenitor sets (states that
    ever appear as a progenitor, i.e., NOT split off directly from root)
  - jaccard_terminal = |A ∩ B| / |A ∪ B|
  - jaccard_all_progenitors = same, but on the all-progenitor union
Also emits sizes and an early-split indicator per sample.

Writes: biology_review/cross_solver_stability.tsv
"""
from __future__ import annotations
import ast
from pathlib import Path
import pandas as pd

ROOT = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing")
CARTA_DIR = ROOT / "analysis/evolution_features_v2/carta_fine12"
GREEDY_DIR = CARTA_DIR / "results_fine12_greedy"
NJ_DIR = CARTA_DIR / "results_fine12_nj"
BEST_K = CARTA_DIR / "best_k_per_sample_v2_chord.csv"
OUT = CARTA_DIR / "biology_review" / "cross_solver_stability.tsv"


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


def build_progen_tree(sets_ordered):
    n = len(sets_ordered)
    parents = [-1] * n
    for i in range(1, n):
        cand_idx, cand_size = 0, len(sets_ordered[0])
        for j in range(i):
            if sets_ordered[j] > sets_ordered[i]:
                if len(sets_ordered[j]) < cand_size or (len(sets_ordered[j]) == cand_size and j < cand_idx):
                    cand_idx, cand_size = j, len(sets_ordered[j])
        parents[i] = cand_idx
    return parents


def get_leaves(parents, n):
    is_parent = set(parents[1:])
    return [i for i in range(n) if i not in is_parent and i > 0]


def analyze_map(res_dir: Path, sample: str, k: int):
    """Return (terminal_union, all_progen_union) as frozensets of state labels."""
    p = res_dir / sample / f"{sample}_progenitors_{k}.txt"
    if not p.exists():
        return None, None
    sets_ord = parse_progen_file(p)
    n = len(sets_ord)
    if n == 0:
        return frozenset(), frozenset()
    parents = build_progen_tree(sets_ord)
    leaves = get_leaves(parents, n)
    if not leaves and n > 1:
        leaves = [n - 1]
    term_union = frozenset().union(*(sets_ord[i] for i in leaves)) if leaves else frozenset()
    all_non_root = frozenset().union(*(sets_ord[i] for i in range(1, n))) if n > 1 else frozenset()
    return term_union, all_non_root


def jaccard(a: frozenset, b: frozenset) -> float:
    if not a and not b:
        return 1.0  # both empty = trivially identical
    u = a | b
    return len(a & b) / len(u) if u else 0.0


def main():
    best = pd.read_csv(BEST_K)
    rows = []
    for _, r in best.iterrows():
        s = r["sample"]
        arm = r["arm"]
        k_star = int(r["k_star_final"])
        gterm, gprog = analyze_map(GREEDY_DIR, s, k_star)
        nterm, nprog = analyze_map(NJ_DIR, s, k_star)
        if gterm is None or nterm is None:
            rows.append({
                "sample": s, "arm": arm, "k_star_final": k_star,
                "terminal_greedy_states": "MISSING" if gterm is None else "",
                "terminal_nj_states": "MISSING" if nterm is None else "",
                "jaccard_terminal": None, "jaccard_all_progenitors": None,
            })
            continue
        jt = jaccard(gterm, nterm)
        jp = jaccard(gprog, nprog)
        rows.append({
            "sample": s,
            "arm": arm,
            "k_star_final": k_star,
            "terminal_greedy_states": ";".join(sorted(gterm)),
            "terminal_nj_states": ";".join(sorted(nterm)),
            "n_terminal_greedy": len(gterm),
            "n_terminal_nj": len(nterm),
            "n_terminal_intersect": len(gterm & nterm),
            "jaccard_terminal": round(jt, 4),
            "all_progenitor_states_greedy": ";".join(sorted(gprog)),
            "all_progenitor_states_nj": ";".join(sorted(nprog)),
            "n_all_prog_greedy": len(gprog),
            "n_all_prog_nj": len(nprog),
            "jaccard_all_progenitors": round(jp, 4),
        })
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"[OK] {len(out)} samples -> {OUT}")
    # summary
    jt_vals = out["jaccard_terminal"].dropna().astype(float).values
    jp_vals = out["jaccard_all_progenitors"].dropna().astype(float).values
    if len(jt_vals):
        import numpy as np
        print(f"\nJaccard(terminal) : median={np.median(jt_vals):.3f}  Q1={np.percentile(jt_vals,25):.3f}  Q3={np.percentile(jt_vals,75):.3f}  min={jt_vals.min():.3f}  max={jt_vals.max():.3f}")
        print(f"Jaccard(all_prog) : median={np.median(jp_vals):.3f}  Q1={np.percentile(jp_vals,25):.3f}  Q3={np.percentile(jp_vals,75):.3f}  min={jp_vals.min():.3f}  max={jp_vals.max():.3f}")
        stable   = out[out["jaccard_terminal"] > 0.7]["sample"].tolist()
        unstable = out[out["jaccard_terminal"] < 0.5]["sample"].tolist()
        print(f"\nSTABLE (jaccard_terminal > 0.7, n={len(stable)}): {stable}")
        print(f"UNSTABLE (jaccard_terminal < 0.5, n={len(unstable)}): {unstable}")
        # {c1, c2, c7} = {Luminal-Distal, Luminal-Proximal, Basal-like} — the shared
        # luminal-epithelial terminal-progenitor claimed on greedy.
        C127 = frozenset({"Luminal-Distal", "Luminal-Proximal", "Basal-like"})
        def contains_c127(cell: str) -> bool:
            if not isinstance(cell, str) or not cell:
                return False
            return C127.issubset(set(cell.split(";")))
        n_g = out["terminal_greedy_states"].map(contains_c127).sum()
        n_n = out["terminal_nj_states"].map(contains_c127).sum()
        print(f"\nSamples with {{c1,c2,c7}} (Luminal-Distal, Luminal-Proximal, Basal-like) SUBSET IN terminal:")
        print(f"  GREEDY: {int(n_g)}/{len(out)} samples")
        print(f"  NJ    : {int(n_n)}/{len(out)} samples")
        modal_g = out["terminal_greedy_states"].value_counts().head(5)
        modal_n = out["terminal_nj_states"].value_counts().head(5)
        print(f"\nTop-5 terminal combos GREEDY:\n{modal_g.to_string()}")
        print(f"\nTop-5 terminal combos NJ:\n{modal_n.to_string()}")


if __name__ == "__main__":
    main()

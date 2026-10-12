#!/usr/bin/env python3
# =============================================================================
# Annotate greedy Cassiopeia trees with edit-distance branch lengths
# =============================================================================
# Description : Load each per-sample greedy Cassiopeia tree (currently all edges
#               = 1) plus its allele table -> character matrix, run Camin-Sokal
#               parsimony ancestral character reconstruction, then set each edge
#               length = Hamming edit distance between parent and child character
#               states (ignoring missing '-1' entries).  Write the re-annotated
#               tree as a newick file with branch lengths.
# Input       : lineage_tracing/allele_tables_per_sample/{sample}_allele_table.csv
#               lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy/{sample}_tree_greedy.nwk
# Output      : analysis/evolution_features_v2/PATH_out/trees_editdist/{sample}_tree_greedy_editdist.nwk
#               analysis/evolution_features_v2/PATH_out/trees_editdist/editdist_summary.csv
# Conda env   : Cassiopeia_v2 (/home/jrchristin/miniconda3/envs/Cassiopeia_v2/bin/python)
# Key deps    : cassiopeia 2.0.0, ete3, pandas, numpy
# Validated   : 2026-09-08
# Parameters  : See DEFAULTS below
# =============================================================================

import argparse, os, sys, time
import numpy as np
import pandas as pd
from ete3 import Tree
import cassiopeia as cas

# ── DEFAULTS ────────────────────────────────────────────────────────────────
DEFAULTS = dict(
    proj      = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing",
    tree_dir  = "lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy",
    tree_suf  = "_tree_greedy.nwk",
    allele_dir= "lineage_tracing/allele_tables_per_sample",
    allele_suf= "_allele_table.csv",
    outdir    = "analysis/evolution_features_v2/PATH_out/trees_editdist",
    # Pseudocount for zero-Hamming edges (parent == child on all reconstructed
    # sites): must be > 0 for PATH's rate-matrix eigendecomposition to be stable.
    # 1e-3 (matching the existing multi2di zero-edge fix) preserves the edit-
    # distance geometry (zero-count edges are ~1000x shorter than a single-
    # mutation edge).
    min_edge  = 1e-3,
    samples   = ["JZ136","JZ201","JZ202","JZ203","JZ204","JZ205","JZ206","JZ207",
                 "JZ208","JZ209","JZ210","JZ211","JZ212","JZ213","JZ214","JZ224",
                 "JZ225","JZ226","JZ227","JZ228","JZ229","JZ230","JZ231","JZ232","JZ234"],
)

# ── HELPERS ─────────────────────────────────────────────────────────────────
def hamming(a, b):
    """Hamming distance between two character-state vectors; missing (-1) ignored."""
    a = np.asarray(a); b = np.asarray(b)
    m = (a != -1) & (b != -1)
    return int((a[m] != b[m]).sum())

def run_sample(sample, p):
    t0 = time.time()
    allele_csv = os.path.join(p["proj"], p["allele_dir"], sample + p["allele_suf"])
    tree_nwk   = os.path.join(p["proj"], p["tree_dir"],   sample + p["tree_suf"])
    if not os.path.exists(allele_csv):
        return dict(sample=sample, status="no_allele", elapsed=None, n_leaves=None, n_edges=None)
    if not os.path.exists(tree_nwk):
        return dict(sample=sample, status="no_tree", elapsed=None, n_leaves=None, n_edges=None)

    # 1) Build character matrix
    at = pd.read_csv(allele_csv)
    cm, _, _ = cas.pp.convert_alleletable_to_character_matrix(at, missing_data_state=-1)
    # 2) Load tree, filter char matrix to tree leaves
    newick = open(tree_nwk).read().strip()
    et = Tree(newick, format=1)
    leaves = [n.name for n in et.get_leaves()]
    have = [x for x in leaves if x in cm.index]
    if len(have) < 30 or len(have) < 0.5 * len(leaves):
        return dict(sample=sample, status="insufficient_char_matrix_overlap",
                    elapsed=round(time.time()-t0, 2), n_leaves=len(leaves), n_edges=None,
                    n_chars=None, n_zero_edges=None, edit_dist_mean=None, edit_dist_median=None,
                    edit_dist_max=None, edit_dist_sum=None)
    cm_sub = cm.loc[have]

    # 3) Load into CassiopeiaTree + reconstruct ancestral characters (Camin-Sokal)
    tree = cas.data.CassiopeiaTree(tree=newick, character_matrix=cm_sub)
    tree.reconstruct_ancestral_characters()

    # 4) Compute per-edge Hamming, set branch lengths (pseudocount for zeros)
    edges = list(tree.edges)
    dists = np.array([hamming(tree.get_character_states(par), tree.get_character_states(chi))
                      for par, chi in edges], dtype=float)
    n_zero = int((dists == 0).sum())
    for i, (par, chi) in enumerate(edges):
        tree.set_branch_length(par, chi, max(float(dists[i]), p["min_edge"]))

    # 5) Write out
    out_nwk = os.path.join(p["outdir"], sample + "_tree_greedy_editdist.nwk")
    os.makedirs(p["outdir"], exist_ok=True)
    with open(out_nwk, "w") as f:
        f.write(tree.get_newick(record_branch_lengths=True) + "\n")

    return dict(sample=sample, status="ok",
                elapsed=round(time.time()-t0, 2),
                n_leaves=tree.n_cell,
                n_chars=tree.n_character,
                n_edges=len(edges),
                n_zero_edges=n_zero,
                edit_dist_mean=float(dists.mean()),
                edit_dist_median=float(np.median(dists)),
                edit_dist_max=float(dists.max()),
                edit_dist_sum=float(dists.sum()))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default=None, help="Run only one sample; default = all 25")
    args = ap.parse_args()
    p = dict(DEFAULTS)
    p["outdir"] = os.path.join(p["proj"], p["outdir"])
    os.makedirs(p["outdir"], exist_ok=True)
    samples = [args.sample] if args.sample else p["samples"]

    rows = []
    for s in samples:
        r = run_sample(s, p)
        rows.append(r)
        print(f"[editdist] {s}: {r['status']}  "
              f"elapsed={r['elapsed']}  "
              f"leaves={r.get('n_leaves')}  edges={r.get('n_edges')}  "
              f"zero_edges={r.get('n_zero_edges')}  "
              f"dist mean={r.get('edit_dist_mean')}  median={r.get('edit_dist_median')}  "
              f"max={r.get('edit_dist_max')}")

    df = pd.DataFrame(rows)
    out_csv = os.path.join(p["outdir"], "editdist_summary.csv")
    df.to_csv(out_csv, index=False)
    n_ok = (df["status"] == "ok").sum()
    print(f"\nsaved: {out_csv}  ({n_ok}/{len(df)} samples OK)")

if __name__ == "__main__":
    main()

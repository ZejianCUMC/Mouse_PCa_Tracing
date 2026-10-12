#!/usr/bin/env python3
# =============================================================================
# Annotate greedy Cassiopeia trees with SCALED PHYLOVELO PSEUDOTIME branch lengths
# =============================================================================
# Description : Third branch-length variant (unit-length, edit-distance already
#               computed). Load each per-sample greedy Cassiopeia tree, attach
#               per-leaf pseudotime from sc_fitness_plasticity_v10.csv
#               (pseudotime_umap column — canonical v10 pseudotime), then
#               reconstruct internal-node pseudotime by continuous parsimony:
#                 - Post-order: internal node pt = MEAN of its (already inferred)
#                   children's pt values.
#                 - Pre-order refinement: enforce monotonicity by clipping any
#                   internal node with pt > min(children pt) down to
#                   min(children pt); leaves are fixed observations and NEVER
#                   moved.
#               Per-edge branch length = Δpt = pt(child) − pt(parent).  Any
#               non-positive Δ is pseudocounted to `min_edge` (default 1e-3),
#               matching the convention used by the unit-length and
#               edit-distance PATH runners.
#               All resulting edge lengths are then divided by the per-sample
#               median-nonzero Δ so that the median edge = 1 pseudotime unit
#               (cross-sample comparability).  The per-sample scale factor is
#               saved to allow reversibility.
# Input       : lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy/{sample}_tree_greedy.nwk
#               + analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv
# Output      : analysis/evolution_features_v2/PATH_out/trees_pseudotime/{sample}_tree_greedy_pseudotime.nwk
#               + analysis/evolution_features_v2/PATH_out/trees_pseudotime/pseudotime_annotate_summary.csv
# Conda env   : hotspot (/home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python)
# Key deps    : ete3, pandas, numpy
# Validated   : 2026-09-08
# =============================================================================

import os, sys, time
import numpy as np
import pandas as pd
from ete3 import Tree

# ── DEFAULTS ────────────────────────────────────────────────────────────────
DEFAULTS = dict(
    proj      = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing",
    tree_dir  = "lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy",
    tree_suf  = "_tree_greedy.nwk",
    scfp_csv  = "analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv",
    pt_col    = "pseudotime_umap",  # v10 canonical pseudotime (per §9 method)
    outdir    = "analysis/evolution_features_v2/PATH_out/trees_pseudotime",
    out_suf   = "_tree_greedy_pseudotime.nwk",
    min_edge  = 1e-3,  # pseudocount for zero / negative Δ pseudotime edges
    samples   = ["JZ136","JZ201","JZ202","JZ203","JZ204","JZ205","JZ206","JZ207",
                 "JZ208","JZ209","JZ210","JZ211","JZ212","JZ213","JZ214","JZ224",
                 "JZ225","JZ226","JZ227","JZ228","JZ229","JZ230","JZ231","JZ232","JZ234"],
)


def annotate_tree(sample: str, pt_map_all: dict, p: dict) -> dict:
    """
    Load greedy topology and annotate branch lengths with scaled pseudotime.
    Returns a dict of QC/summary numbers.
    """
    t0 = time.time()
    tree_nwk = os.path.join(p["proj"], p["tree_dir"], sample + p["tree_suf"])
    if not os.path.exists(tree_nwk):
        return dict(sample=sample, status="no_tree_file")

    # Load newick (format=1 permits internal-node names)
    t = Tree(open(tree_nwk).read().strip(), format=1)

    leaves = t.get_leaves()
    n_leaves_raw = len(leaves)

    # Per-sample pseudotime map (cellBC → pt)
    pt_map = pt_map_all.get(sample, {})
    if len(pt_map) == 0:
        return dict(sample=sample, status="no_pt_map_for_sample", n_leaves_raw=n_leaves_raw)

    n_leaves_with_pt = 0
    n_leaves_without_pt = 0
    to_prune = []
    for leaf in leaves:
        if leaf.name in pt_map:
            leaf.add_feature("pt", float(pt_map[leaf.name]))
            n_leaves_with_pt += 1
        else:
            to_prune.append(leaf.name)
            n_leaves_without_pt += 1

    # Prune leaves without pseudotime (mirrors PATH runners which drop
    # unlabeled tips downstream, so we do the same here to keep tree consistent)
    if to_prune:
        keep = [l for l in leaves if l.name not in set(to_prune)]
        if len(keep) < 2:
            return dict(sample=sample, status="insufficient_leaves_with_pt",
                        n_leaves_raw=n_leaves_raw, n_leaves_with_pt=n_leaves_with_pt,
                        n_leaves_without_pt=n_leaves_without_pt)
        t.prune(keep, preserve_branch_length=False)

    # Post-order MEAN — internal node pt = mean of children pt
    for node in t.traverse("postorder"):
        if node.is_leaf():
            continue
        child_pts = [c.pt for c in node.children]
        node.add_feature("pt", float(np.mean(child_pts)))

    # Pre-order refinement — enforce monotonicity
    # For each internal node (root→leaves), clip node.pt down to min(children.pt)
    # so parent.pt <= all its children.pt.  Leaves are NEVER modified.
    # Iterate to ~convergence (usually 1 pass is enough, but 3 caps pathological cases).
    for _ in range(3):
        changed = 0
        for node in t.traverse("preorder"):
            if node.is_leaf() or not node.children:
                continue
            min_child = min(c.pt for c in node.children)
            if node.pt > min_child:
                node.pt = min_child
                changed += 1
        if changed == 0:
            break

    # Compute per-edge Δpt.  Root has no parent edge.
    edges_raw = []
    for node in t.traverse():
        if node.up is None:  # root
            continue
        dpt = node.pt - node.up.pt
        # We should have dpt >= 0 after monotonicity refinement, but sanity-check.
        # Any residual non-positive value → pseudocount.
        edges_raw.append(max(dpt, p["min_edge"]))

    edges_raw = np.array(edges_raw, dtype=float)
    n_pseudocount = int((edges_raw <= p["min_edge"]).sum())
    n_edges = len(edges_raw)

    # Median-scale: divide by median of non-pseudocount edges so median-real-edge = 1
    non_ps_edges = edges_raw[edges_raw > p["min_edge"]]
    if non_ps_edges.size < 2:
        # Degenerate — everything at floor. Fall back to overall median (= min_edge).
        scale = float(np.median(edges_raw))
    else:
        scale = float(np.median(non_ps_edges))
    if scale <= 0:
        scale = 1.0

    edges_scaled = edges_raw / scale

    # Re-assign scaled branch lengths back to the tree.  Traverse in same order.
    i = 0
    for node in t.traverse():
        if node.up is None:
            continue
        node.dist = float(edges_scaled[i])
        i += 1

    # Write
    os.makedirs(p["outdir_abs"], exist_ok=True)
    out_nwk = os.path.join(p["outdir_abs"], sample + p["out_suf"])
    t.write(format=1, outfile=out_nwk)

    return dict(
        sample=sample,
        status="ok",
        elapsed=round(time.time() - t0, 2),
        n_leaves_raw=n_leaves_raw,
        n_leaves_with_pt=n_leaves_with_pt,
        n_leaves_without_pt=n_leaves_without_pt,
        n_edges=n_edges,
        n_pseudocount_edges=n_pseudocount,
        median_dpt_raw=float(np.median(non_ps_edges)) if non_ps_edges.size else 0.0,
        mean_dpt_raw=float(edges_raw.mean()),
        max_dpt_raw=float(edges_raw.max()),
        scale_factor=scale,
        out_nwk=out_nwk,
    )


def main():
    p = dict(DEFAULTS)
    p["outdir_abs"] = os.path.join(p["proj"], p["outdir"])
    os.makedirs(p["outdir_abs"], exist_ok=True)

    scfp = pd.read_csv(os.path.join(p["proj"], p["scfp_csv"]))
    print(f"[annotate_pt] sc_fp = {len(scfp):,} rows / {scfp['sample'].nunique()} samples "
          f"/ pt col={p['pt_col']} range=[{scfp[p['pt_col']].min():.3f},"
          f"{scfp[p['pt_col']].max():.3f}] NaN={scfp[p['pt_col']].isna().sum()}")

    # Build per-sample cellBC → pt map for O(1) lookup
    pt_map_all = {}
    for sample, sub in scfp.groupby("sample"):
        pt_map_all[sample] = dict(zip(sub["cellBC"].astype(str), sub[p["pt_col"]].astype(float)))

    rows = []
    for s in p["samples"]:
        r = annotate_tree(s, pt_map_all, p)
        rows.append(r)
        # Pretty log line
        keys = ["sample","status","elapsed","n_leaves_raw","n_leaves_with_pt",
                "n_leaves_without_pt","n_edges","n_pseudocount_edges",
                "median_dpt_raw","scale_factor"]
        msg = " ".join(f"{k}={r.get(k)}" for k in keys)
        print(f"[annotate_pt] {msg}")

    df = pd.DataFrame(rows)
    csv_out = os.path.join(p["outdir_abs"], "pseudotime_annotate_summary.csv")
    df.to_csv(csv_out, index=False)
    n_ok = int((df["status"] == "ok").sum())
    print(f"\nsaved: {csv_out} ({n_ok}/{len(df)} samples OK)")


if __name__ == "__main__":
    main()

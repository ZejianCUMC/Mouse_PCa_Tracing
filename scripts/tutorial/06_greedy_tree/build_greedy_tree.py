#!/usr/bin/env python
# =============================================================================
# 06 · Per-mouse lineage tree: allele table + cell states -> VanillaGreedy tree
# =============================================================================
# Description : The TSL library has no hashtag, so the run-level allele table
#               mixes all mice. Mouse identity + tumour state come from the
#               scRNA-seq side (step 04) by joining on the 16-bp cell barcode.
#               For each mouse:
#                 1. keep allele rows of that mouse's tumour cells,
#                 2. drop cells flagged exclude_from_tree (Cycling / Ionocyte),
#                 3. character matrix: cells x (intBC x r1/r2/r3); 0 = uncut,
#                    -1 = missing, 1..n = indel states
#                    (convert_alleletable_to_character_matrix, allele_rep_thresh 0.99),
#                 4. Cassiopeia VanillaGreedySolver (the lab's primary solver;
#                    NJ / MaxCut / Spectral / Percolation only as robustness checks),
#                 5. write Newick + character matrix + leaf states + QC + figure.
#               Greedy trees have unit branch lengths and many polytomies — expected.
# Input       : ${OUT_CAS}/${TSL_LIB}_allele_table.csv, ${OUT_ST}/cell_states.csv
# Output      : ${OUT_TREE}/<SAMPLE>/  <SAMPLE>_tree_greedy.nwk, character_matrix.csv,
#               leaf_states.csv, tree_qc.json, <SAMPLE>_tree_greedy.png
# Conda env   : HPC ${ENVS}/cassiopeia_env  (Cassiopeia 2.1)
# Validated   : 2026-09-27
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import os, json, argparse, time
import numpy as np, pandas as pd
import cassiopeia as cas
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
DEFAULTS = {
    "allele_rep_thresh": 0.99,    # drop intBCs where one allele covers >99 % of cells (uninformative)
    "min_cells": 30,              # skip a mouse with fewer lineage-covered cells
    "solver": "greedy",           # greedy (primary) | nj | maxcut | spectral | percolation
    "drop_excluded_states": True, # remove Cycling / Ionocyte before solving
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def get_solver(name):
    return {"greedy": lambda: cas.solver.VanillaGreedySolver(),
            "nj": lambda: cas.solver.NeighborJoiningSolver(add_root=True),
            "maxcut": lambda: cas.solver.MaxCutGreedySolver(),
            "spectral": lambda: cas.solver.SpectralSolver(),
            "percolation": lambda: cas.solver.PercolationSolver(joining_solver=cas.solver.VanillaGreedySolver()),
            }[name]()

def matrix_qc(cm, at):
    """Missing / uncut / edited fractions of the character matrix + per-cut-site edit rates
    (per-site rates come from the allele table: '[None]' = uncut)."""
    v = cm.values
    q = {"cells": int(cm.shape[0]), "characters_kept": int(cm.shape[1]),
         "intBCs_in_allele_table": int(at["intBC"].nunique()),
         "frac_missing": float((v == -1).mean()), "frac_uncut": float((v == 0).mean()),
         "frac_edited": float((v > 0).mean())}
    for site in ["r1", "r2", "r3"]:
        x = at[site].dropna().astype(str)
        q[f"{site}_edited_of_observed"] = float((~x.str.contains(r"\[None\]")).mean()) if len(x) else None
    return q

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str, output_dir: str, **kwargs) -> dict:
    """input_path = allele table; kwargs must include states=<cell_states.csv>, sample=<JZ id>."""
    p = {**DEFAULTS, **kwargs}
    sample = p["sample"]
    out = os.path.join(output_dir, sample); os.makedirs(out, exist_ok=True)

    at = pd.read_csv(input_path)
    at["cellBC"] = at["cellBC"].str.split("-").str[0]              # 16-bp join key
    st = pd.read_csv(p["states"])
    st = st[st["sample_id"] == sample]
    n_rna = len(st)
    if p["drop_excluded_states"]:
        st = st[~st["exclude_from_tree"].astype(bool)]
    at = at.merge(st[["cellBC", "sample_id", "treatment", "state"]], on="cellBC", how="inner")
    n_cells = at["cellBC"].nunique()
    print(f"{sample}: tumour cells (scRNA) {n_rna} | after state filter {len(st)} | with lineage {n_cells}")
    if n_cells < p["min_cells"]:
        raise SystemExit(f"{sample}: only {n_cells} cells with lineage data — not solving")

    cm, priors, state2indel = cas.pp.convert_alleletable_to_character_matrix(
        at, allele_rep_thresh=p["allele_rep_thresh"])
    cm.to_csv(os.path.join(out, "character_matrix.csv"))
    qc = matrix_qc(cm, at)
    qc.update({"sample": sample, "tumour_cells_scRNA": n_rna, "cells_with_lineage": n_cells,
               "lineage_recovery": n_cells / max(len(st), 1)})

    tree = cas.data.CassiopeiaTree(character_matrix=cm, priors=priors)
    t0 = time.time()
    get_solver(p["solver"]).solve(tree, logfile=os.path.join(out, "solver.log"))
    qc["solve_seconds"] = round(time.time() - t0, 1)
    nwk_fp = os.path.join(out, f"{sample}_tree_{p['solver']}.nwk")
    with open(nwk_fp, "w") as f:
        f.write(tree.get_newick(record_branch_lengths=True))
    leaves = at.drop_duplicates("cellBC").set_index("cellBC").loc[tree.leaves, ["sample_id", "treatment", "state"]]
    leaves.to_csv(os.path.join(out, "leaf_states.csv"))
    qc.update({"leaves": len(tree.leaves), "internal_nodes": len(tree.internal_nodes)})
    depths = [len(tree.get_all_ancestors(l)) for l in tree.leaves]
    qc["max_depth"], qc["mean_depth"] = int(max(depths)), float(np.mean(depths))
    json.dump(qc, open(os.path.join(out, "tree_qc.json"), "w"), indent=2)
    print(json.dumps(qc, indent=2))

    try:                                                            # quick-look figure
        tree.cell_meta = leaves[["state"]].astype(str)
        fig, ax = cas.pl.plot_matplotlib(tree, meta_data=["state"], figsize=(10, 10))
        ax.set_title(f"{sample} {p['solver']} tree ({len(tree.leaves)} cells)")
        plt.savefig(os.path.join(out, f"{sample}_tree_{p['solver']}.png"), dpi=150, bbox_inches="tight")
        plt.close("all")
    except Exception as e:
        print("plotting skipped:", e)
    return {"output_path": nwk_fp, "params_used": p, "summary": qc}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--allele_table", default=os.path.join(os.environ.get("OUT_CAS", "."),
                    f"{os.environ.get('TSL_LIB', '')}_allele_table.csv"))
    ap.add_argument("--states", default=os.path.join(os.environ.get("OUT_ST", "."), "cell_states.csv"))
    ap.add_argument("--sample", required=True)
    ap.add_argument("--solver", default=DEFAULTS["solver"])
    ap.add_argument("--outdir", default=os.environ.get("OUT_TREE"))
    a = ap.parse_args()
    run_pipeline(a.allele_table, a.outdir, states=a.states, sample=a.sample, solver=a.solver)

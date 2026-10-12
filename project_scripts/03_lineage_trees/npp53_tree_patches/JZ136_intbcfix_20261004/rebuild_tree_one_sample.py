#!/usr/bin/env python
"""[intbcfix 2026-10-04] copy of analysis/hpc_tree_rebuild/rebuild_tree_one_sample.py (v8); ONLY input/output paths changed: input = error-intBC-free JZ136 table, outputs -> lineage_tracing/JZ136_intbcfix_20261004/. Methods/params identical.
Rebuild a Cassiopeia tree for ONE sample with ONE solver.
   Ionocyte + Cycling cellBCs are filtered from the allele table before solving.
   Usage: rebuild_tree_one_sample.py <SAMPLE> <SOLVER> [<OUTDIR>]
   SOLVER ∈ {greedy, maxcut, spectral, percolation, nj}"""
import os, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')

# Auto-detect env (local vs HPC)
LOCAL_ROOT = '/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing'
HPC_ROOT   = '/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing'
ROOT = LOCAL_ROOT if os.path.isdir(LOCAL_ROOT) else HPC_ROOT
EXCLUDE_STATES = {'Ionocyte','Cycling'}

def build_character_matrix(alt_df, min_intbc_thresh=0.05):
    """Convert allele table to character matrix using Cassiopeia's data preprocess helper."""
    import cassiopeia as cas
    cm, priors, state_to_indel = cas.pp.convert_alleletable_to_character_matrix(
        alt_df, allele_rep_thresh=0.99)
    return cm, priors, state_to_indel

def solve(cm, priors, solver_name):
    import cassiopeia as cas
    tree = cas.data.CassiopeiaTree(character_matrix=cm, priors=priors)
    if   solver_name == 'greedy':      s = cas.solver.VanillaGreedySolver()
    elif solver_name == 'maxcut':      s = cas.solver.MaxCutGreedySolver()
    elif solver_name == 'spectral':    s = cas.solver.SpectralSolver()
    elif solver_name == 'percolation': s = cas.solver.PercolationSolver(joining_solver=cas.solver.VanillaGreedySolver())
    elif solver_name == 'nj':          s = cas.solver.NeighborJoiningSolver(add_root=True)
    else: raise ValueError(f'unknown solver {solver_name}')
    s.solve(tree, logfile='/dev/null')
    return tree

def main():
    if len(sys.argv) < 3:
        print('Usage: rebuild_tree_one_sample.py SAMPLE SOLVER [OUTDIR]', file=sys.stderr); sys.exit(1)
    sample, solver = sys.argv[1], sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else f'{ROOT}/cassiopeia_trees_v8_rebuilt/{solver}'
    Path(outdir).mkdir(parents=True, exist_ok=True)

    at_fp = f'/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/lineage_tracing/JZ136_intbcfix_20261004/input/{sample}_allele_table.csv'
    print(f'[env] ROOT={ROOT}', flush=True)
    print(f'[env] SAMPLE={sample} SOLVER={solver} OUTDIR={outdir}', flush=True)
    print(f'[load] allele table: {at_fp}', flush=True)
    at = pd.read_csv(at_fp)
    n_cells_before = at['cellBC'].nunique()
    at = at[~at['putative_state'].isin(EXCLUDE_STATES)]
    n_cells = at['cellBC'].nunique()
    print(f'[filter] cellBC: {n_cells_before} → {n_cells} (Ionocyte+Cycling removed)', flush=True)
    if n_cells < 30:
        print(f'[abort] too few cells after filter ({n_cells})', file=sys.stderr); sys.exit(2)

    print(f'[build] character matrix', flush=True)
    t0 = time.time()
    cm, priors, _ = build_character_matrix(at)
    print(f'  cm: {cm.shape}  ({time.time()-t0:.1f}s)', flush=True)

    print(f'[solve] {solver}', flush=True)
    t0 = time.time()
    tr = solve(cm, priors, solver)
    print(f'  solve: {time.time()-t0:.1f}s', flush=True)

    print(f'[save] newick', flush=True)
    nwk = tr.get_newick(record_branch_lengths=True)
    out_fp = f'{outdir}/{sample}_tree_{solver}.nwk'
    with open(out_fp, 'w') as f: f.write(nwk)
    print(f'  saved: {out_fp}  ({len(tr.leaves)} leaves)', flush=True)

if __name__ == '__main__':
    main()

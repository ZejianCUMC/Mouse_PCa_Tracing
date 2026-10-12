#!/usr/bin/env python
"""Run jungle LBI fitness on ONE sample's rebuilt Greedy tree.
   Uses cassiopeia.tools.fitness_estimator.LBIJungle (vendored Py3-compatible jungle).
   Outputs: {sample}_lbi_fitness.csv  +  {sample}_mean_fitness.nwk (with mean_fitness attr).
   Usage: run_jungle_one_sample.py <SAMPLE>"""
import os, sys, warnings
warnings.filterwarnings('ignore')

LOCAL_ROOT = '/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing'
HPC_ROOT   = '/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing'
ROOT = LOCAL_ROOT if os.path.isdir(LOCAL_ROOT) else HPC_ROOT

def main():
    sample = sys.argv[1]
    tree_fp = f'{ROOT}/cassiopeia_trees_v8_rebuilt/greedy/{sample}_tree_greedy.nwk'
    outdir  = f'{ROOT}/cassiopeia_trees_v8_rebuilt/jungle_fitness'
    os.makedirs(outdir, exist_ok=True)
    print(f'[env] ROOT={ROOT}  SAMPLE={sample}', flush=True)
    print(f'[load] tree: {tree_fp}', flush=True)

    import cassiopeia as cas
    from cassiopeia.tools.fitness_estimator import LBIJungle
    import pandas as pd

    tree = cas.data.CassiopeiaTree(tree=tree_fp)
    n_leaves = len(tree.leaves)
    print(f'  leaves={n_leaves}, nodes={len(tree.nodes)}', flush=True)

    est = LBIJungle(random_seed=42)
    print(f'[infer] LBIJungle.estimate_fitness ...', flush=True)
    est.estimate_fitness(tree)

    # Extract fitness per node (skip root, per LBIJungle docs)
    rows = []
    for node in tree.nodes:
        try:
            fitness = tree.get_attribute(node, 'fitness')
        except Exception:
            continue
        rows.append(dict(node=str(node), is_leaf=tree.is_leaf(node), fitness=float(fitness)))
    df = pd.DataFrame(rows)
    out_csv = f'{outdir}/{sample}_lbi_fitness.csv'
    df.to_csv(out_csv, index=False)
    print(f'  saved: {out_csv}  ({len(df)} rows)', flush=True)
    leaf_df = df[df['is_leaf']]
    print(f'  leaf fitness stats: min={leaf_df["fitness"].min():.4f}, max={leaf_df["fitness"].max():.4f}, mean={leaf_df["fitness"].mean():.4f}', flush=True)

    # Also save NHX-style newick with fitness attribute (compatible with tree viewer)
    nwk_out = f'{outdir}/{sample}_mean_fitness.nwk'
    with open(nwk_out, 'w') as f:
        f.write(tree.get_newick(record_branch_lengths=True))
    print(f'  saved: {nwk_out}', flush=True)

if __name__ == '__main__':
    main()

#!/usr/bin/env python
# =============================================================================
# Prune Ionocyte tips from ALL Cassiopeia trees per Zejian's directive (2026-09-01)
# =============================================================================
# For every existing tree in:
#   cassiopeia_trees_alt/{sample}_tree_{maxcut,spectral,percolation}.nwk
#   cassiopeia_trees_output/{sample}_tree.nwk        (Greedy solver output)
#   cassiopeia_trees_ilp_v3/{sample}_tree_ilp.nwk
#   cassiopeia_phylovelo_out/{sample}_pseudotime_scaled.newick
# Drop leaves whose cellBC is labeled 'Ionocyte' in sc_fitness_plasticity.csv.
# Original trees preserved with '.ionOK.bak' suffix; new trees replace originals.
# =============================================================================
import os, shutil, sys
from pathlib import Path
import pandas as pd
from ete3 import Tree

M    = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
SCFP = M/'1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity.csv'

TREE_DIRS_PATTERNS = [
    ('cassiopeia_trees_alt', '{sample}_tree_maxcut.nwk'),
    ('cassiopeia_trees_alt', '{sample}_tree_spectral.nwk'),
    ('cassiopeia_trees_alt', '{sample}_tree_percolation.nwk'),
    ('cassiopeia_trees_output', '{sample}_tree.nwk'),
    ('cassiopeia_trees_ilp_v3', '{sample}_tree_ilp.nwk'),
    ('cassiopeia_phylovelo_out', '{sample}_pseudotime_scaled.newick'),
]

def main():
    sc_fp = pd.read_csv(SCFP)
    ion_by_sample = (sc_fp[sc_fp['putative_state'] == 'Ionocyte']
                       .groupby('sample')['cellBC'].apply(set).to_dict())
    samples = sorted(sc_fp['sample'].unique())
    print(f'[main] {len(samples)} samples, Ionocyte counts: {sum(len(v) for v in ion_by_sample.values())} total')
    rows = []
    for sample in samples:
        ion_bcs = ion_by_sample.get(sample, set())
        for subdir, pat in TREE_DIRS_PATTERNS:
            fp = M/subdir/pat.format(sample=sample)
            if not fp.exists(): continue
            try:
                tr = Tree(str(fp), format=1)
            except Exception as e:
                print(f'  {fp.name}: PARSE ERROR {e}'); continue
            n_before = len(tr.get_leaves())
            # find leaves matching Ionocyte cellBCs
            leaves_to_prune = [l.name for l in tr.get_leaves() if l.name in ion_bcs]
            n_prune = len(leaves_to_prune)
            if n_prune == 0:
                rows.append(dict(sample=sample, tree_file=fp.name, n_before=n_before, n_pruned=0, n_after=n_before))
                continue
            # backup original if not already backed up
            bak = fp.with_suffix(fp.suffix + '.ionOK.bak')
            if not bak.exists():
                shutil.copy(fp, bak)
            # ete3 prune: keep = all leaves NOT in ion_bcs
            keep = [l.name for l in tr.get_leaves() if l.name not in ion_bcs]
            if len(keep) < 5:
                print(f'  {fp.name}: SKIP (only {len(keep)} non-Ion leaves)'); continue
            tr.prune(keep, preserve_branch_length=True)
            n_after = len(tr.get_leaves())
            # save
            tr.write(outfile=str(fp), format=1)
            rows.append(dict(sample=sample, tree_file=fp.name, n_before=n_before,
                              n_pruned=n_prune, n_after=n_after))
            print(f'  {sample} · {fp.name}: {n_before} → {n_after} (-{n_prune} Ion)')
    df = pd.DataFrame(rows)
    df.to_csv(M/'1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/tree_prune_ionocyte_summary.csv', index=False)
    print(f'\n[done] pruned {len(df[df.n_pruned>0])} trees; summary saved.')
    print(f'\nTotal Ionocyte tips removed: {df["n_pruned"].sum()}')
    print(f'Trees processed: {len(df)}')

if __name__ == '__main__':
    main()

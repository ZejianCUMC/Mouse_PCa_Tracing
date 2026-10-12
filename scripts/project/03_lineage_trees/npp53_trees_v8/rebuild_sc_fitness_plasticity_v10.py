#!/usr/bin/env python
# Rebuild sc_fitness_plasticity_v10.csv from jungle LBI (v8-rebuilt trees, unchanged)
# + v10 PhyloVelo pseudotime (RECOMPUTED shared UMAP/scVI/diffmap on 18,241 no-cyc-no-ion cells).
# Fitch parsimony plasticity per subclade also recomputed on the rebuilt trees.
# Only pseudotime differs from v9 — fitness + plasticity are identical.
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import pandas as pd, numpy as np
from ete3 import Tree

M = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
OUT_DIR = M/'1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity'
OUT_DIR.mkdir(exist_ok=True)
JUNGLE_DIR   = M/'cassiopeia_trees_v8_rebuilt/jungle_fitness'
PHYLO_DIR    = M/'cassiopeia_trees_v8_rebuilt/phylovelo_out'                # now v10
GREEDY_DIR   = M/'cassiopeia_trees_v8_rebuilt/greedy'
OLD_SCFP     = M/'1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity.csv'
NEW_SCFP     = OUT_DIR/'sc_fitness_plasticity_v10.csv'

old = pd.read_csv(OLD_SCFP)
print(f'old sc_fp: {len(old)} rows, cols={list(old.columns)}')

SAMPLES = sorted([p.stem.replace('_lbi_fitness', '') for p in JUNGLE_DIR.glob('*_lbi_fitness.csv')])
print(f'samples with jungle fitness: {len(SAMPLES)}')

def compute_fitch_plasticity_per_cell(tree, state_map, min_size=20):
    node_set = {}
    for n in tree.traverse('postorder'):
        if n.is_leaf():
            s = state_map.get(n.name)
            node_set[n] = {s} if s else set()
        else:
            child_sets = [node_set[c] for c in n.children if node_set[c]]
            if not child_sets: node_set[n] = set(); continue
            inter = set.intersection(*child_sets)
            node_set[n] = inter if inter else set.union(*child_sets)
    node_state = {}
    for n in tree.traverse('preorder'):
        if n.is_root():
            node_state[n] = next(iter(node_set[n])) if node_set[n] else None
        else:
            ps = node_state.get(n.up)
            if ps and ps in node_set[n]: node_state[n] = ps
            else: node_state[n] = next(iter(node_set[n])) if node_set[n] else ps
    cell_plast_sum, cell_plast_n = {}, {}
    for node in tree.traverse('postorder'):
        if node.is_leaf(): continue
        leaves = node.get_leaves()
        if len(leaves) < min_size: continue
        trans, edges = 0, 0
        for sub in node.traverse():
            if sub is node or sub.up is None: continue
            if node_state[sub] is not None and node_state[sub.up] is not None:
                edges += 1
                if node_state[sub] != node_state[sub.up]: trans += 1
        if edges == 0: continue
        plast = trans / edges
        for l in leaves:
            cell_plast_sum[l.name] = cell_plast_sum.get(l.name, 0) + plast
            cell_plast_n[l.name]   = cell_plast_n.get(l.name, 0) + 1
    return {c: cell_plast_sum[c] / cell_plast_n[c] for c in cell_plast_sum}

rows = []
for sample in SAMPLES:
    tree_fp = GREEDY_DIR/f'{sample}_tree_greedy.nwk'
    if not tree_fp.exists():
        print(f'  {sample}: SKIP — no greedy tree'); continue
    tr = Tree(str(tree_fp), format=1)
    tree_tips = {l.name for l in tr.get_leaves()}
    jf = pd.read_csv(JUNGLE_DIR/f'{sample}_lbi_fitness.csv')
    jf_leaf = jf[jf['is_leaf']].set_index('node')['fitness'].to_dict()
    # v10 phylovelo pseudotime (UMAP-default + diffmap)
    pt_umap_fp = PHYLO_DIR/f'{sample}_phylovelo_pseudotime.csv'
    pt_dm_fp   = PHYLO_DIR/f'{sample}_phylovelo_pseudotime_diffmap.csv'
    pt_umap = pd.read_csv(pt_umap_fp, index_col=0).iloc[:, 0].to_dict() if pt_umap_fp.exists() else {}
    pt_dm   = pd.read_csv(pt_dm_fp,   index_col=0).iloc[:, 0].to_dict() if pt_dm_fp.exists() else {}
    sub = old[old['sample']==sample].drop_duplicates('cellBC').set_index('cellBC')
    state_map = {l.name: sub['putative_state'].get(l.name) for l in tr.get_leaves() if l.name in sub.index}
    cell_plast = compute_fitch_plasticity_per_cell(tr, state_map)
    n_before = len(rows)
    for bc in tree_tips:
        if bc not in sub.index: continue
        rows.append(dict(
            cellBC=bc, sample=sample,
            treatment=sub.loc[bc, 'treatment'],
            putative_state=sub.loc[bc, 'putative_state'],
            fitness=jf_leaf.get(bc, np.nan),
            plasticity=cell_plast.get(bc, np.nan),
            pseudotime_umap=pt_umap.get(bc, np.nan),
            pseudotime_diffmap=pt_dm.get(bc, np.nan)
        ))
    print(f'  {sample}: {len(rows)-n_before} cells matched   (pt_umap: {len(pt_umap)}, pt_dm: {len(pt_dm)})')

new = pd.DataFrame(rows)
new.to_csv(NEW_SCFP, index=False)
print(f'\n[saved] {NEW_SCFP} ({len(new)} rows, {list(new.columns)})')
print(f'\nsummary by treatment:')
print(new.groupby('treatment').size().to_string())
print(f'\nNaN counts: fitness={new.fitness.isna().sum()}, plasticity={new.plasticity.isna().sum()}, '
      f'pt_umap={new.pseudotime_umap.isna().sum()}, pt_dm={new.pseudotime_diffmap.isna().sum()}')

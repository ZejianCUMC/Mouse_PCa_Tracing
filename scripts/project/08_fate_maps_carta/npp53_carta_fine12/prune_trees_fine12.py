#!/usr/bin/env python
"""
For each fine12 label TSV, prune the greedy .nwk tree to only leaves whose cell survived
the per-sample >=20-cell fine-12-state filter. CARTA errors out with KeyError on any tree
leaf not present in the label TSV, so an aligned tree is required.

Uses ete3 for pruning (branch lengths preserved).

Writes: inputs_fine12/{sample}_tree_greedy_fine12.nwk
Also writes: inputs_fine12/prune_summary.tsv (n_original, n_kept, n_pruned)
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from ete3 import Tree

ROOT = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing")
TREE_DIR = ROOT / "lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy"
IN_DIR = ROOT / "analysis/evolution_features_v2/carta_fine12/inputs_fine12"
SAMPLES = (IN_DIR / "samples.txt").read_text().strip().splitlines()

rows = []
for s in SAMPLES:
    src = TREE_DIR / f"{s}_tree_greedy.nwk"
    labels = pd.read_csv(IN_DIR / f"{s}_labels_fine12.tsv", sep='\t')
    keep = set(labels['cellBC'])
    t = Tree(src.read_text(), format=1)
    leaves = t.get_leaves()
    leaf_names = [l.name for l in leaves]
    n_orig = len(leaf_names)
    keep_in_tree = [n for n in leaf_names if n in keep]
    n_kept = len(keep_in_tree)
    n_pruned = n_orig - n_kept
    t.prune(keep_in_tree, preserve_branch_length=True)
    out = IN_DIR / f"{s}_tree_greedy_fine12.nwk"
    out.write_text(t.write(format=1))
    check = Tree(out.read_text(), format=1)
    n_after = len(check.get_leaves())
    ok = (n_after == n_kept)
    print(f"[{s}] orig={n_orig}  kept={n_kept}  pruned={n_pruned}  written_leaves={n_after}  ok={ok}")
    rows.append({'sample': s, 'n_original_leaves': n_orig, 'n_kept_leaves': n_kept,
                 'n_pruned': n_pruned, 'n_leaves_in_pruned_tree': n_after, 'ok': ok})

pd.DataFrame(rows).to_csv(IN_DIR / 'prune_summary.tsv', sep='\t', index=False)
print(f"\n[OK] pruned trees + prune_summary.tsv written to {IN_DIR}")

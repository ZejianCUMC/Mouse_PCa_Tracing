# Compare JZ272 original (558 cells) vs state>=20 test (528 cells): pseudotime spread, clade recovery, LBI.
import json, numpy as np, pandas as pd
from pathlib import Path
from ete3 import Tree
from scipy.stats import spearmanr
R = Path(__file__).resolve().parent; RL = R.parent; S = "JZ272"
OC, OV = RL / "cassiopeia_mjz026_028_20260925_seed42", RL / "tree_viewer_mjz026_028_20260925_seed42"
told = Tree(open(OC / "results" / S / f"{S}_tree_greedy.nwk").read(), format=1); tnew = Tree(open(R / "cassiopeia/results" / S / f"{S}_tree_greedy.nwk").read(), format=1)
shared = set(tnew.get_leaf_names()); told.prune(list(shared), preserve_branch_length=True)
def clades(t, mn=10):
    return [frozenset(n.get_leaf_names()) for n in t.traverse() if not n.is_leaf() and mn <= len(n) < len(t)]
co, cn = clades(told), clades(tnew)
best = [max(len(a & b) / len(a | b) for b in cn) for a in co]
big = [b for a, b in zip(co, best) if len(a) >= 50]
pt_o = pd.read_csv(OV / "results" / S / "PhyloVelo/pseudotime.csv", index_col=0).pseudotime
pt_n = pd.read_csv(R / "viewer/results" / S / "PhyloVelo/pseudotime.csv", index_col=0).pseudotime
def desc(p): q = p.quantile([.01, .25, .5, .75, .99]); s = np.sort(p.values); return dict(n=len(p), unique=p.round(4).nunique(), q01=q[.01], q25=q[.25], q50=q[.5], q75=q[.75], q99=q[.99], IQR=q[.75]-q[.25], gap_bottom=s[1]-s[0])
c = pt_o.index.intersection(pt_n.index)
lo = pd.read_csv(OV / "results" / S / "LBI/lbi_fitness.csv"); ln = pd.read_csv(R / "viewer/results" / S / "LBI/lbi_fitness.csv")
lo, ln = lo[lo.is_leaf].set_index("node").fitness, ln[ln.is_leaf].set_index("node").fitness; cl = lo.index.intersection(ln.index)
cj = json.loads((R / "viewer/results" / S / "PhyloVelo/COMPLETE.json").read_text())
res = dict(cells_old=558, cells_new=len(shared),
           clades_ge10_old=len(co), clades_ge10_new=len(cn), clade_best_jaccard_median=float(np.median(best)), clades_ge10_recovered_J_ge_0_8=float(np.mean(np.array(best) >= .8)),
           clades_ge50_n=len(big), clades_ge50_recovered_J_ge_0_8=float(np.mean(np.array(big) >= .8)) if big else None,
           pt_old=desc(pt_o), pt_new=desc(pt_n), pt_spearman_shared=float(spearmanr(pt_o[c], pt_n[c])[0]),
           lbi_spearman_shared=float(spearmanr(lo[cl], ln[cl])[0]), phylovelo_new_k=cj["n_neighbors"], phylovelo_new_graph=cj["graph_checks"], phylovelo_new_valid=cj["pseudotime_valid"])
(R / "COMPARISON.json").write_text(json.dumps(res, indent=2, default=float)); print(json.dumps(res, indent=1, default=lambda v: round(float(v), 3)))

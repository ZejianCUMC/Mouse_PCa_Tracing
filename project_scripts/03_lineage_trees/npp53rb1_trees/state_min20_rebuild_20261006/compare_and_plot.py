# Current (original MJZ026-028 run) vs state>=20 rebuild for one sample: clade recovery, pseudotime spread, LBI; side-by-side heatmap.
import sys, json, numpy as np, pandas as pd
from pathlib import Path
from ete3 import Tree
from scipy.stats import spearmanr
from PIL import Image, ImageDraw, ImageFont
S = sys.argv[1]; G = Path(__file__).resolve().parent; R = G / S; RL = G.parent
OC, OV = RL / "cassiopeia_mjz026_028_20260925_seed42", RL / "tree_viewer_mjz026_028_20260925_seed42"
told = Tree(open(OC / "results" / S / f"{S}_tree_greedy.nwk").read(), format=1); tnew = Tree(open(R / "cassiopeia/results" / S / f"{S}_tree_greedy.nwk").read(), format=1)
n_old = len(told); shared = set(tnew.get_leaf_names()); told.prune(list(shared), preserve_branch_length=True)
clades = lambda t, mn: [frozenset(n.get_leaf_names()) for n in t.traverse() if not n.is_leaf() and mn <= len(n) < len(t)]
co, cn = clades(told, 10), clades(tnew, 10); best = np.array([max(len(a & b) / len(a | b) for b in cn) for a in co]); big = np.array([b for a, b in zip(co, best) if len(a) >= 50])
pt_o = pd.read_csv(OV / "results" / S / "PhyloVelo/pseudotime.csv", index_col=0).pseudotime; pt_n = pd.read_csv(R / "viewer/results" / S / "PhyloVelo/pseudotime.csv", index_col=0).pseudotime
def desc(p): q = p.quantile([.01, .25, .5, .75, .99]); s = np.sort(p.values); return dict(unique=int(p.round(4).nunique()), q01=q[.01], q25=q[.25], q75=q[.75], q99=q[.99], IQR=q[.75] - q[.25], gap_bottom=s[1] - s[0])
c = pt_o.index.intersection(pt_n.index)
lo = pd.read_csv(OV / "results" / S / "LBI/lbi_fitness.csv"); ln = pd.read_csv(R / "viewer/results" / S / "LBI/lbi_fitness.csv")
lo, ln = lo[lo.is_leaf].set_index("node").fitness, ln[ln.is_leaf].set_index("node").fitness; cl = lo.index.intersection(ln.index)
cj = json.loads((R / "viewer/results" / S / "PhyloVelo/COMPLETE.json").read_text()); oj = json.loads((OV / "results" / S / "PhyloVelo/COMPLETE.json").read_text())
res = dict(sample=S, cells_old=n_old, cells_new=len(shared), clades_ge10_old=len(co), clades_ge10_new=len(cn), clade_best_jaccard_median=float(np.median(best)),
           clades_ge10_recovered_J_ge_0_8=float((best >= .8).mean()), clades_ge50_n=len(big), clades_ge50_recovered_J_ge_0_8=float((big >= .8).mean()) if len(big) else None,
           pt_old=desc(pt_o), pt_new=desc(pt_n), pt_spearman_shared=float(spearmanr(pt_o[c], pt_n[c])[0]), lbi_spearman_shared=float(spearmanr(lo[cl], ln[cl])[0]),
           phylovelo_k_old=oj["n_neighbors"], phylovelo_k_new=cj["n_neighbors"], phylovelo_new_graph=cj["graph_checks"], phylovelo_new_valid=cj["pseudotime_valid"])
(R / "COMPARISON.json").write_text(json.dumps(res, indent=2, default=float))
old_png = next((RL.parent / "figures/barcode_heatmaps_v4_20261005").glob(f"Np53Rb1_*_{S}_barcode_heatmap.png")); new_png = next((R / "heatmap").glob(f"Np53Rb1_*_{S}_barcode_heatmap.png"))
A, B = Image.open(old_png).convert("RGB"), Image.open(new_png).convert("RGB"); h = 1400
A = A.resize((round(A.width * h / A.height), h)); B = B.resize((round(B.width * h / B.height), h))
cv = Image.new("RGB", (A.width + B.width + 40, h + 80), "white"); d = ImageDraw.Draw(cv); f = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 40)
cv.paste(A, (0, 80)); cv.paste(B, (A.width + 40, 80)); d.text((10, 15), f"a  Current ({n_old} cells)", fill="black", font=f); d.text((A.width + 50, 15), f"b  States >= 20 cells ({len(shared)} cells)", fill="black", font=f)
cv.save(R / f"{S}_current_vs_min20.png")
r = lambda v: round(v, 3)
print(S, f"cells {n_old}->{len(shared)} | clades>=50 recovered {res['clades_ge50_recovered_J_ge_0_8']:.2f} (n={len(big)}) | clade J median {r(res['clade_best_jaccard_median'])} | "
      f"pt IQR {r(res['pt_old']['IQR'])}->{r(res['pt_new']['IQR'])} q01 {r(res['pt_old']['q01'])}->{r(res['pt_new']['q01'])} gap {r(res['pt_old']['gap_bottom'])}->{r(res['pt_new']['gap_bottom'])} "
      f"uniq {res['pt_old']['unique']}->{res['pt_new']['unique']} | pt rho {r(res['pt_spearman_shared'])} | LBI rho {r(res['lbi_spearman_shared'])} | k {oj['n_neighbors']}->{cj['n_neighbors']}")

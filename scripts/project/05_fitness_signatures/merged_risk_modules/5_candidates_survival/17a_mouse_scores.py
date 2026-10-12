"""opt5 step 17a — per-tumour mouse module scores for the user-rule passers (usage: 17a_mouse_scores.py <arm>), for plotting (17b).
Same scoring as step 16: scanpy score_genes (log1p CP10k, random_state 42) on the arm's cells, per-tumour mean, z-scored over tumours.
Output <arm>/user_rule_tumour_scores.csv (cand, sample_id, model, treatment, n_cells, z). Conda env: moscot."""
import sys, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
O = Path(__file__).resolve().parent; A = O / sys.argv[1]; MERGED = O.parents[5]
R = pd.read_csv(A / "user_rule_passers.csv"); cells = pd.read_csv(A / "cells_n3000_PRIVATE.csv", index_col=0)
sh = ad.read_h5ad(MERGED / "shared_counts.h5ad"); sh = sh[cells.index.values].copy(); sh.var_names = sh.var.gene_symbol.astype(str).values; sh.var_names_make_unique()
sc.pp.normalize_total(sh, target_sum=1e4); sc.pp.log1p(sh); out = []
for r in R.itertuples():
    sc.tl.score_genes(sh, [g for g in r.genes.split() if g in sh.var_names], score_name="s", random_state=42)
    d = cells.assign(s=sh.obs["s"].values).groupby("sample_id").agg(s=("s", "mean"), n_cells=("s", "size"), treatment=("treatment", "first"), model=("model", "first")).reset_index()
    d["z"] = (d.s - d.s.mean()) / d.s.std(); d["cand"] = r.cand; out.append(d[["cand", "sample_id", "model", "treatment", "n_cells", "z"]])
pd.concat(out).to_csv(A / "user_rule_tumour_scores.csv", index=False); print(sys.argv[1], "tumour scores for", len(R), "modules")

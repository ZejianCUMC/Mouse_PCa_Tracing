"""opt5 step 04 — Hotspot inputs per arm and cell set (usage: 04_hotspot_inputs.py <arm> <cellset>; cellset = Npp53 | Npp53Rb1 |
Npp53Rb1-n2000 | Npp53Rb1-n1000). Cells = that genotype's LBI cells of the arm (cycOUT non-cycling / cycKEPT all), JZ135 LBI from the
named tree (n3000 default). Genes = union of the arm's kept lists of that cell set (+ <arm>/hotspot_inputs/<cellset>_extra_genes.txt).
Plain npz/csv/npy (the hotspot env's anndata 0.8 cannot read the merged h5ad): <cellset>_bundle_PRIVATE.npz (raw counts + gene names,
written atomically), <cellset>_obs_PRIVATE.csv (cell, total_counts = all-gene UMI), <cellset>_latent_PRIVATE.npy (merged X_scVI).
Conda env: moscot."""
import os, sys, json, numpy as np, pandas as pd, anndata as ad
from pathlib import Path
O = Path(__file__).resolve().parent; RV10 = O.parents[1]; MERGED = O.parents[5]; ARM, cs = sys.argv[1], sys.argv[2]; A = O / ARM
D = A / "hotspot_inputs"; D.mkdir(exist_ok=True)
L = json.load(open(A / "fsu_lists.json")); genes = {g for d in L.values() if d["cells"] == cs for g in d["genes"]}
if (D / f"{cs}_extra_genes.txt").exists(): genes |= set((D / f"{cs}_extra_genes.txt").read_text().split())
genes = sorted(genes); geno, v = (cs.split("-") + ["n3000"])[:2]
c = pd.read_csv(A / f"cells_{v}_PRIVATE.csv", index_col=0)
if geno != "all": c = c[c.model == geno]                                      # "all" = both genotypes (pooled meta-analysis lists)
assert c.sample_id.nunique() == {"Npp53": 25, "Npp53Rb1": 22, "all": 47}[geno]
sh = ad.read_h5ad(MERGED / "shared_counts.h5ad"); sh = sh[c.index.values].copy(); sym = sh.var.gene_symbol.astype(str).values
gi = [i for i, x in enumerate(sym) if x in set(genes)]; assert len({sym[i] for i in gi}) == len(gi) == len(genes), "gene mismatch"
Xg = sh.X[:, gi].tocsr(); tmp = D / f"{cs}_bundle_PRIVATE.tmp.npz"
np.savez(tmp, data=Xg.data, indices=Xg.indices, indptr=Xg.indptr, shape=np.array(Xg.shape), genes=np.array(sym[gi], dtype=str)); os.replace(tmp, D / f"{cs}_bundle_PRIVATE.npz")
if not (D / f"{cs}_obs_PRIVATE.csv").exists():
    pd.DataFrame({"cell": c.index, "total_counts": np.asarray(sh.X.sum(1)).ravel()}).to_csv(D / f"{cs}_obs_PRIVATE.csv", index=False)
    g = ad.read_h5ad(RV10 / "named_graph_PRIVATE.h5ad"); np.save(D / f"{cs}_latent_PRIVATE.npy", pd.DataFrame(g.obsm["X_scVI"], index=g.obs_names).loc[c.index].values.astype(np.float64))
print(ARM, cs, "cells", len(c), "tumours", c.sample_id.nunique(), "genes", len(gi))

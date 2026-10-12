"""opt5 step 21 — meta-analysis workflow, step 2 of 3: per-tumour fitness effect sizes with subclone-robust standard errors
(usage: 21_meta_effects.py <arm> [lbi_transform: blom | z]). All cells (cycling kept), Greedy-tree LBI (JZ135 = n3000 tree).
Per tumour and gene (detected in > 10 cells; tracing transgenes excluded):
  log1p(counts) ~ 1 + genes-detected + f,  f = LBI rank-inverse-normal within the tumour (Blom; default) or z-scored LBI
  beta   = expression change per 1 SD of within-tumour fitness (comparable across tumours)
  se_cr  = cluster-robust (CR1) standard error with clusters = Greedy-tree subclones (step 20): cells of one subclone share lineage and
           LBI and are not independent; Var(beta) = K/(K-1) * (n-1)/(n-p) * sum_k (sum_{i in k} a_i e_i)^2, a = row of (X'X)^-1 X'
  se_ols = classical OLS standard error (for comparison only)
Expression (env EXPR): raw (default; log1p of raw UMI, as the original FitnessSignature) | cp10k (Zejian 2026-10-10: "u need to use
normalized counts": log1p(counts / total UMI per cell x 1e4), scanpy normalize_total + log1p; genes-detected covariate kept).
Output <arm>/effects.csv.gz (unit, gene, beta, se_cr, se_ols, n_cells, n_clades, n_expr). Deterministic. Conda env: moscot."""
import os, sys, time, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from pathlib import Path
from multiprocessing import Pool
from scipy import stats
O = Path(__file__).resolve().parent; MERGED = O.parents[5]; A = O / sys.argv[1]; TRANS = sys.argv[2] if len(sys.argv) > 2 else "blom"
TG = {"mCherry", "Cas9", "GFP", "Cre"}
cells = pd.read_csv(A / "cells_n3000_PRIVATE.csv", index_col=0); cl = pd.read_csv(A / "clades_PRIVATE.csv").set_index("cell").clade
sh = ad.read_h5ad(MERGED / "shared_counts.h5ad"); genes = sh.var.gene_symbol.astype(str).values; isTG = np.isin(genes, list(TG))
row = pd.Series(np.arange(sh.n_obs), index=sh.obs_names); X = sh.X.tocsr(); SF = np.asarray((X > 0).sum(1)).ravel().astype(float)
import os; EXPR = os.environ.get("EXPR", "raw"); TOT = np.asarray(X.sum(1)).ravel().astype(float)
units = [(s, cells.index[cells.sample_id == s], cells.loc[cells.sample_id == s, "LBI"].values) for s in sorted(cells.sample_id.unique())]
def tf(f):
    if TRANS == "z": return (f - f.mean()) / f.std()
    r = stats.rankdata(f); return stats.norm.ppf((r - .375) / (len(f) + .25))
def one(u):
    name, idx, f = u; ii = row[idx].values; yc = X[ii]; n = len(ii)
    nexp = np.asarray((yc > 0).sum(0)).ravel(); gm = (nexp > 10) & ~isTG; yd = yc[:, gm].toarray().astype(np.float64)
    if EXPR == "cp10k": yd = yd / TOT[ii][:, None] * 1e4      # library-size normalised (total UMI over all genes)
    y = np.log1p(yd)
    D = np.column_stack([np.ones(n), SF[ii], tf(f)]); p = D.shape[1]
    XtXi = np.linalg.inv(D.T @ D); Bm = XtXi @ D.T @ y; E = y - D @ Bm; a = (XtXi @ D.T)[2]
    codes, uniq = pd.factorize(cl.loc[idx].values); K = len(uniq)
    C = sp.csr_matrix((np.ones(n), (codes, np.arange(n))), shape=(K, n))
    S = C @ (a[:, None] * E); var_cr = (S ** 2).sum(0) * (K / (K - 1)) * ((n - 1) / (n - p))
    var_ols = (E ** 2).sum(0) / (n - p) * XtXi[2, 2]
    return pd.DataFrame(dict(unit=name, gene=genes[gm], beta=Bm[2], se_cr=np.sqrt(var_cr), se_ols=np.sqrt(var_ols), n_cells=n, n_clades=K, n_expr=nexp[gm]))
if __name__ == "__main__":
    t0 = time.time()
    with Pool(int(os.environ.get("NPROC", 10))) as pool: R = pool.map(one, units, chunksize=1)
    D = pd.concat(R, ignore_index=True); D.to_csv(A / "effects.csv.gz", index=False, compression="gzip")
    s = D.groupby("unit").agg(genes=("gene", "size"), n_cells=("n_cells", "first"), n_clades=("n_clades", "first"), se_ratio_median=("se_cr", lambda x: np.nan))
    D["ratio"] = D.se_cr / D.se_ols; s["se_cr_over_ols_median"] = D.groupby("unit").ratio.median().round(2)
    print(s.drop(columns=["se_ratio_median"]).to_string()); print(f"{sys.argv[1]} ({TRANS}): units {D.unit.nunique()}, rows {len(D)}, {time.time()-t0:.0f}s")

#!/usr/bin/env python
# =============================================================================
# Job 0 — build pySCENIC raw-count input loom for Npp53 tumor cells (no Cyc/Ion)
# =============================================================================
# Mirrors the karthaus_v3clean_pyscenic/00_build_loom.py convention.
# Source cells: 18,162 canonical tumor cells from sc_fitness_plasticity_v10.csv
#   (Ionocyte + Cycling already excluded upstream; also ~79 additional QC drops).
# Expression : npp53_edited.h5ad (20,978 x 33,676), raw counts in layers['counts'];
#              .X is log-normed and MUST NOT be used for SCENIC.
# Filter key : composite (obs['sample_id'], bare_bc) where bare_bc =
#              obs.index.split('_')[1].replace('-1','').  Matches scfp
#              ('sample', 'cellBC') 100% (18162/18162 verified).
# Gene filter: pySCENIC noise filter — drop genes detected in <3 cells.
# Output     : npp53_tumor_no_cyc_ion.loom  (rows=Gene, cols=CellID)
# Env        : pysenic  (loompy + scanpy)
# =============================================================================
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, loompy as lp

H5AD    = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/data/npp53_merge/npp53_edited.h5ad"
SCFP    = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv"
LOOMOUT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/inputs/npp53_tumor_no_cyc_ion.loom"

print("[1] load scfp cellBC list (canonical tumor cells)", flush=True)
scfp = pd.read_csv(SCFP)
scfp_key = (scfp["sample"].astype(str) + "|" + scfp["cellBC"].astype(str)).values
scfp_set = set(scfp_key)
print(f"    scfp n={len(scfp)}  n_unique_keys={len(scfp_set)}", flush=True)

print("[2] load npp53_edited.h5ad in memory + build composite key", flush=True)
adata = sc.read_h5ad(H5AD)
bc  = adata.obs.index.to_series().str.split("_").str[1].str.replace("-1","",regex=False)
sid = adata.obs["sample_id"].astype(str)
adata.obs["_key"] = (sid + "|" + bc).values
print(f"    h5ad shape={adata.shape}  n_matched={adata.obs['_key'].isin(scfp_set).sum()}", flush=True)

print("[3] swap X to raw counts layer (SCENIC needs raw ints)", flush=True)
assert "counts" in adata.layers, "no 'counts' layer in h5ad!"
adata.X = adata.layers["counts"].copy()
sx = adata.X[:200].toarray() if sp.issparse(adata.X) else adata.X[:200]
print(f"    X max={float(sx.max()):.0f}  integer={(sx == sx.astype(int)).all()}", flush=True)

print("[4] subset to scfp tumor cells", flush=True)
adata = adata[adata.obs["_key"].isin(scfp_set)].copy()
print(f"    after cell filter: {adata.shape}", flush=True)
assert adata.n_obs == len(scfp_set), f"filter drop: {len(scfp_set)-adata.n_obs}"

print("[5] pySCENIC gene noise filter — drop genes detected in <3 cells", flush=True)
v0 = adata.n_vars
sc.pp.filter_genes(adata, min_cells=3)
print(f"    genes {v0} -> {adata.n_vars}", flush=True)

print("[6] write pySCENIC loom (v2 spec, rows=Gene cols=CellID)", flush=True)
X  = adata.X
Xd = X.toarray() if sp.issparse(X) else np.asarray(X)
matrix = Xd.T.astype(np.float32)                 # genes x cells
ng = np.asarray((adata.X > 0).sum(axis=1)).flatten().astype(np.int32)
nu = np.asarray(adata.X.sum(axis=1)).flatten().astype(np.float32)
row_attrs = {"Gene":   np.array(adata.var_names, dtype="U")}
col_attrs = {"CellID": np.array(adata.obs_names, dtype="U"), "nGene": ng, "nUMI": nu}
lp.create(LOOMOUT, matrix, row_attrs, col_attrs)
print(f"[OK] {LOOMOUT}  ({matrix.shape[0]} genes x {matrix.shape[1]} cells)", flush=True)

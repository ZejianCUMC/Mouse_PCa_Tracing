#!/usr/bin/env python
"""
00_prep_data.py — Load h5ad + scfp whitelist, filter to 12-state tumor cells.

Emits:
  - _cache/adata_12state.h5ad   (18,162 cells × 33,676 genes, X=log-norm, layers['counts']=raw, obsm['scVI'])
  - _cache/state_mean_lognorm.parquet  (12 × 33,676 mean log-norm expression per state)
  - _cache/state_mean_scvi.parquet     (12 × 30 mean scVI per state)
  - _cache/cell_scvi.parquet           (18,162 × 30 scVI + putative_state per cell)
"""
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

ROOT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
H5AD = f"{ROOT}/data/npp53_merge/npp53_edited.h5ad"
SCFP = f"{ROOT}/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv"
OUT = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/_cache"
os.makedirs(OUT, exist_ok=True)

STATE_ORDER = [
    "Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep", "Luminal-CR",
    "Luminal-Persistent", "Luminal-Squamous",
    "Basal-like",
    "Partial-EMT",
    "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF",
]

def main():
    print("[00] loading h5ad…", flush=True)
    h = ad.read_h5ad(H5AD)
    print("  h5ad shape:", h.shape)
    print("[00] loading scfp whitelist…", flush=True)
    scfp = pd.read_csv(SCFP)
    print("  scfp shape:", scfp.shape)

    # ── build (sample, bare_bc) join key ─────────────────────────────────
    obs = h.obs.copy()
    obs["bare_bc"] = obs.index.to_series().str.split("_").str[1].str.replace("-1", "", regex=False)
    obs["_key"] = obs["sample_id"].astype(str) + "|" + obs["bare_bc"]
    scfp["_key"] = scfp["sample"].astype(str) + "|" + scfp["cellBC"].astype(str)
    whitelist = set(scfp["_key"].values)
    keep = obs["_key"].isin(whitelist).values
    print(f"  h5ad cells matched to scfp whitelist: {keep.sum()}/{len(obs)}  (expect 18,162)")
    adata = h[keep].copy()

    # ── verify state distribution ────────────────────────────────────────
    vc = adata.obs["putative_state"].value_counts()
    print("[00] filtered per-state counts:")
    print(vc)
    missing = [s for s in STATE_ORDER if s not in vc.index]
    if missing:
        print(f"WARN: missing states: {missing}")
    # drop any unused categories
    adata.obs["putative_state"] = adata.obs["putative_state"].astype(str)
    adata = adata[adata.obs["putative_state"].isin(STATE_ORDER)].copy()
    print(f"  final adata shape: {adata.shape}")

    adata.write(f"{OUT}/adata_12state.h5ad")
    print(f"  wrote {OUT}/adata_12state.h5ad")

    # ── state-mean log-norm expression (12 × G) ──────────────────────────
    print("[00] state-mean log-norm expression…", flush=True)
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    df_cells = pd.DataFrame(X, index=adata.obs.index, columns=adata.var_names)
    df_cells["_state"] = adata.obs["putative_state"].values
    state_mean = df_cells.groupby("_state").mean().reindex(STATE_ORDER)
    state_mean.to_parquet(f"{OUT}/state_mean_lognorm.parquet")
    print(f"  wrote state_mean_lognorm.parquet  shape={state_mean.shape}")

    # ── state-mean scVI (12 × 30) ────────────────────────────────────────
    print("[00] state-mean scVI…", flush=True)
    scvi_raw = adata.obsm["scVI"]
    scvi_arr = scvi_raw.to_numpy() if hasattr(scvi_raw, "to_numpy") else np.asarray(scvi_raw)
    cols = [f"scVI_{i:02d}" for i in range(scvi_arr.shape[1])]
    df_scvi = pd.DataFrame(scvi_arr, index=adata.obs.index.astype(str), columns=cols)
    df_scvi["_state"] = adata.obs["putative_state"].values
    state_mean_scvi = df_scvi.groupby("_state").mean().reindex(STATE_ORDER)
    state_mean_scvi.to_parquet(f"{OUT}/state_mean_scvi.parquet")
    print(f"  wrote state_mean_scvi.parquet  shape={state_mean_scvi.shape}")

    # ── cell-level scVI + state (for chooseR / bootstrap) ────────────────
    df_scvi.to_parquet(f"{OUT}/cell_scvi.parquet")
    print(f"  wrote cell_scvi.parquet  shape={df_scvi.shape}")

    print("[00] DONE")

if __name__ == "__main__":
    main()

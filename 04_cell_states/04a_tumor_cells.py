#!/usr/bin/env python
# =============================================================================
# 04a · Major cell types + reporter-positive tumour cells  (scanpy / scVI)
# =============================================================================
# Description : Grafted tumours contain host cells (immune, endothelial,
#               fibroblasts). Tumour cells carry the lineage-recorder transgenes,
#               read out as GFP / mCherry "genes" in the GRCm39_tracing reference.
#               Steps: counts -> scVI latent -> Leiden -> per-cluster compartment
#               call by canonical marker scores (argmax of z-scored means) ->
#               host cluster = host-type markers AND mostly reporter-negative ->
#               tumour = GFP>0 AND mCherry>0 AND not in a host cluster (the Npp53
#               "edited / double-positive, host-removed" definition).
# Input       : ${OUT_QC}/counts.mtx, genes.txt, barcodes.txt, cell_meta.csv
# Output      : ${OUT_ST}/all_cells.h5ad, tumor_cells.h5ad, major_types.csv, plots
# Conda env   : HPC ${ENVS}/scvi   (scanpy 1.11, scvi-tools 1.2)
# Key deps    : scanpy, scvi-tools, leidenalg
# Validated   : 2026-09-27
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import os, argparse
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp
import anndata as ad, scanpy as sc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
DEFAULTS = {
    "n_hvg": 2000,               # HVGs for scVI (seurat_v3 on raw counts)
    "n_latent": 30,              # scVI latent dims (Npp53 cohort setting)
    "n_layers": 2,               # scVI encoder/decoder layers
    "max_epochs": 200,           # scVI training epochs (early stopping on)
    "n_neighbors": 15,           # kNN graph on the scVI latent
    "leiden_res": 0.5,           # coarse clustering for major cell types
    "min_markers": 2,            # a compartment needs >= this many markers present
    "reporter_genes": ["GFP", "mCherry"],   # tumour = all of these > 0 (CellBender counts)
    "exclude_from_hvg": ["GFP", "mCherry", "Cas9", "Cas9_full", "Cas9-full", "Cre", "Cre_full", "Cre-full"],  # transgenes never drive clustering (Seurat turns "_" into "-")
    "host_types": ["Fibroblast", "SmoothMuscle", "Endothelial", "Macrophage_Myeloid", "Tcell_NK", "Bcell"],
    "host_max_reporter_frac": 0.3,   # host cluster = host-type markers AND < 30 % reporter+ cells
                                     # (tumour mesenchyme can look fibroblast-like but is reporter+)
    "seed": int(os.environ.get("SEED", 42)),
}

# Canonical mouse-prostate compartments (prior knowledge, used only to NAME clusters)
MARKERS = {
    "Epithelial_Luminal": ["Epcam", "Krt8", "Krt18", "Krt19", "Cdh1", "Pbsn", "Nkx3-1", "Ar", "Krt7", "Cldn3"],
    "Epithelial_Basal":   ["Krt5", "Krt14", "Trp63", "Krt15"],
    "Mesenchymal_EMT":    ["Vim", "Col3a1", "Col1a1", "Zeb1", "Zeb2", "Fn1", "Igfbp2", "Twist1", "Snai2"],
    "Fibroblast":         ["Dcn", "Lum", "Pdgfra", "Gsn", "Col1a2", "Sfrp2"],
    "SmoothMuscle":       ["Acta2", "Myh11", "Tagln", "Myl9", "Des"],
    "Endothelial":        ["Pecam1", "Cldn5", "Cdh5", "Flt1", "Emcn", "Egfl7"],
    "Macrophage_Myeloid": ["Lyz2", "C1qa", "C1qb", "C1qc", "Adgre1", "Cd68", "Itgam", "Csf1r"],
    "Tcell_NK":           ["Cd3e", "Cd3d", "Cd3g", "Cd8a", "Nkg7", "Gzmb", "Klrb1c"],
    "Bcell":              ["Cd79a", "Cd79b", "Ms4a1", "Igkc"],
    "Neuroendocrine":     ["Chga", "Chgb", "Syp", "Ncam1", "Ascl1"],
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def load_counts(qc_dir):
    m = sio.mmread(os.path.join(qc_dir, "counts.mtx")).tocsr().T.tocsr()      # -> cells x genes
    genes = open(os.path.join(qc_dir, "genes.txt")).read().split("\n")[:m.shape[1]]
    bcs = open(os.path.join(qc_dir, "barcodes.txt")).read().split("\n")[:m.shape[0]]
    meta = pd.read_csv(os.path.join(qc_dir, "cell_meta.csv")).set_index("cell").loc[bcs]
    a = ad.AnnData(X=sp.csr_matrix(m, dtype=np.float32), obs=meta, var=pd.DataFrame(index=genes))
    a.var_names_make_unique()
    a.layers["counts"] = a.X.copy()
    return a

def run_scvi(a, p):
    """HVG (transgenes removed) -> scVI latent in a.obsm['X_scVI'] -> kNN + UMAP."""
    import scvi
    scvi.settings.seed = p["seed"]
    h = a[:, ~a.var_names.isin(p["exclude_from_hvg"])].copy()
    sc.pp.highly_variable_genes(h, n_top_genes=p["n_hvg"], flavor="seurat_v3", layer="counts")
    h = h[:, h.var["highly_variable"]].copy()
    # one 10x run -> no batch key. Never use treatment or mouse as batch (they are the biology).
    scvi.model.SCVI.setup_anndata(h, layer="counts")
    model = scvi.model.SCVI(h, n_latent=p["n_latent"], n_layers=p["n_layers"])
    model.train(max_epochs=p["max_epochs"], early_stopping=True)
    a.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(a, use_rep="X_scVI", n_neighbors=p["n_neighbors"], random_state=p["seed"])
    sc.tl.umap(a, random_state=p["seed"])
    return model

def compartment_scores(a, markers, min_markers):
    """Mean of per-gene z-scored log-normalised expression, per compartment (cells x compartments)."""
    out = {}
    for name, genes in markers.items():
        g = [x for x in genes if x in a.var_names]
        if len(g) < min_markers:
            continue
        X = a[:, g].X
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        sd = X.std(axis=0); sd[sd == 0] = 1
        out[name] = ((X - X.mean(axis=0)) / sd).mean(axis=1)
    return pd.DataFrame(out, index=a.obs_names)

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str, output_dir: str, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    os.makedirs(output_dir, exist_ok=True)
    sc.settings.figdir = output_dir
    a = load_counts(input_path)
    print("cells x genes:", a.shape)

    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)          # a.X = log-normalised
    run_scvi(a, p)
    sc.tl.leiden(a, resolution=p["leiden_res"], key_added="leiden_major",
                 random_state=p["seed"], flavor="igraph", n_iterations=2, directed=False)

    S = compartment_scores(a, MARKERS, p["min_markers"])
    for c in S.columns:
        a.obs[f"score_{c}"] = S[c].values
    per_cl = S.groupby(a.obs["leiden_major"].values).mean()
    cl2type = per_cl.idxmax(axis=1)
    a.obs["major_type"] = a.obs["leiden_major"].map(cl2type).astype("category")
    print("\ncluster -> major type\n", pd.concat([cl2type.rename("major_type"),
          a.obs["leiden_major"].value_counts().rename("n")], axis=1))

    cnt = a.layers["counts"]
    pos = np.ones(a.n_obs, dtype=bool)
    for g in p["reporter_genes"]:
        v = np.asarray(cnt[:, a.var_names.get_loc(g)].todense()).ravel()
        a.obs[f"{g}_counts"] = v
        pos &= v > 0
    a.obs["reporter_pos"] = pos
    frac = a.obs.groupby("leiden_major", observed=True)["reporter_pos"].mean()
    host_cl = [c for c in frac.index if cl2type[c] in p["host_types"] and frac[c] < p["host_max_reporter_frac"]]
    a.obs["host_cluster"] = a.obs["leiden_major"].isin(host_cl)
    a.obs["is_tumor"] = pos & ~a.obs["host_cluster"]
    print("\nper cluster: marker call, reporter+ fraction, host?\n",
          pd.DataFrame({"major_type": cl2type, "reporter_frac": frac.round(3),
                        "host": [c in host_cl for c in cl2type.index]}))
    print("\ntumour cells per mouse\n", a.obs.groupby("sample_id")["is_tumor"].sum())

    for col in ["leiden_major", "major_type", "sample_id", "treatment", "is_tumor"]:
        sc.pl.umap(a, color=col, show=False, title=f"scVI UMAP: {col}")
        plt.savefig(os.path.join(output_dir, f"UMAP_all_{col}.png"), dpi=150, bbox_inches="tight"); plt.close()
    sc.pl.umap(a, color=[f"{g}_counts" for g in p["reporter_genes"]], vmax="p99", show=False)
    plt.savefig(os.path.join(output_dir, "UMAP_all_reporters.png"), dpi=150, bbox_inches="tight"); plt.close()

    a.obs[["sample_id", "treatment", "leiden_major", "major_type", "reporter_pos", "host_cluster", "is_tumor"]]\
        .to_csv(os.path.join(output_dir, "major_types.csv"))
    a.write(os.path.join(output_dir, "all_cells.h5ad"))
    t = a[a.obs["is_tumor"]].copy()
    t.X = t.layers["counts"].copy()                     # 04b starts again from raw counts
    for k in list(t.obsm.keys()): del t.obsm[k]
    t.write(os.path.join(output_dir, "tumor_cells.h5ad"))
    summary = {"n_cells": int(a.n_obs), "n_tumor": int(t.n_obs)}
    print(summary)
    return {"output_path": output_dir, "params_used": p, "summary": summary}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=os.environ.get("OUT_QC"))
    ap.add_argument("--outdir", default=os.environ.get("OUT_ST"))
    ap.add_argument("--leiden_res", type=float, default=DEFAULTS["leiden_res"])
    a = ap.parse_args()
    run_pipeline(a.input, a.outdir, leiden_res=a.leiden_res)

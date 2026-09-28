#!/usr/bin/env python
# =============================================================================
# 04b · Tumour cell states: scVI + Leiden resolution sweep + stability + markers
# =============================================================================
# Description : Re-embeds tumour cells only (new scVI), clusters with Leiden at
#               resolutions 0.1 ... 1.2 and asks, for every resolution:
#                 * k            number of clusters
#                 * ARI_boot     bootstrap reproducibility (80 % subsample, re-cluster, B runs)
#                 * jac_mean/min per-cluster stability = best Jaccard match of each
#                                full-data cluster in the bootstrap clustering
#               A resolution is "stable" when every cluster's median Jaccard >= jac_min.
#               The finest stable resolution is suggested (override with --res).
#               At the chosen resolution: Wilcoxon markers, dotplot, and a
#               signature-score hint for each cluster (Luminal / Basal / EMT /
#               Cycling / Ionocyte ...). The hint is a STARTING POINT for naming —
#               final state names come from reading the marker table.
#               Cycling and Ionocyte clusters are flagged `exclude_from_tree`
#               (the Npp53 recipe removes them before tree building).
#               This is the teaching version. The locked Npp53 procedure is
#               CHOIR-anchored AMI -> top candidates -> chooseR + sc-SHC (HANDSOFF.md §6).
# Input       : ${OUT_ST}/tumor_cells.h5ad (raw counts)
# Output      : ${OUT_ST}/resolution_metrics.csv, cell_states.csv, markers_res*.csv,
#               tumor_states.h5ad, UMAP/dotplot/stability figures
# Conda env   : HPC ${ENVS}/scvi
# Validated   : 2026-09-27
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import os, argparse
import numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
DEFAULTS = {
    "n_hvg": 2000,                       # HVGs for the tumour-only scVI
    "n_latent": 30, "n_layers": 2,       # scVI architecture (Npp53 cohort)
    "max_epochs": 300,                   # scVI epochs (early stopping on)
    "n_neighbors": 15,                   # kNN on scVI latent
    "resolutions": [round(0.1 * i, 1) for i in range(1, 13)],   # Leiden sweep 0.1..1.2
    "n_boot": 20,                        # bootstrap runs per resolution
    "boot_frac": 0.8,                    # fraction of cells per bootstrap
    "jac_min": 0.6,                      # stability threshold (median per-cluster Jaccard)
    "min_cluster_cells": 20,             # smaller clusters are flagged "tiny" in the table
    "exclude_from_hvg": ["GFP", "mCherry", "Cas9", "Cas9-full", "Cre", "Cre-full"],
    "exclude_hints": ["Cycling", "Ionocyte"],   # removed before tree building
    "seed": int(os.environ.get("SEED", 42)),
}

# Candidate tumour-state signatures (mouse prostate; prior knowledge used only for naming hints)
SIGNATURES = {
    "Luminal_AR":     ["Ar", "Nkx3-1", "Pbsn", "Fkbp5", "Tmprss2", "Krt8", "Krt18"],
    "Luminal_ARlow":  ["Krt8", "Krt18", "Krt4", "Tacstd2", "Psca"],
    "Luminal_Squamous": ["Krt13", "Krt4", "Sprr1a", "Krt6a", "Krt17"],
    "Basal":          ["Krt5", "Krt14", "Trp63", "Krt15"],
    "Cycling":        ["Mki67", "Top2a", "Ccnb1", "Cdk1", "Birc5"],
    "Partial_EMT":    ["Vim", "Cd44", "Krt8", "Vcan", "Sparc"],
    "Mesenchymal":    ["Col3a1", "Col1a1", "Zeb1", "Zeb2", "Fn1", "Igfbp2"],
    "Ionocyte":       ["Foxi1", "Ascl3", "Cftr", "Atp6v1g3", "Atp6v0d2"],
    "Neuroendocrine": ["Chga", "Chgb", "Syp", "Ascl1", "Ncam1"],
    "Interferon":     ["Isg15", "Ifit1", "Ifit3", "Irf7", "Stat1"],
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def leiden(a, res, seed, key):
    sc.tl.leiden(a, resolution=res, key_added=key, random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    return a.obs[key].astype(str).values

def per_cluster_jaccard(full, boot):
    """For every full-data cluster: best Jaccard overlap with any bootstrap cluster."""
    out = {}
    for c in np.unique(full):
        A = set(np.where(full == c)[0])
        out[c] = max(len(A & set(np.where(boot == b)[0])) / len(A | set(np.where(boot == b)[0]))
                     for b in np.unique(boot))
    return out

def signature_scores(a, sigs, min_genes=2):
    """Per signature: mean of per-gene z-scored log-normalised expression (cells x signatures).
    Same deterministic scoring as the lab's annotate_npp53.R; rare states (e.g. ionocytes)
    stand out because each gene is scaled across the tumour cells."""
    out = {}
    for name, genes in sigs.items():
        g = [x for x in genes if x in a.var_names]
        if len(g) < min_genes:
            continue
        X = a[:, g].X
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        sd = X.std(axis=0); sd[sd == 0] = 1
        out[name] = ((X - X.mean(axis=0)) / sd).mean(axis=1)
        a.obs[f"sig_{name}"] = out[name]
    return pd.DataFrame(out, index=a.obs_names)

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str, output_dir: str, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    import scvi
    scvi.settings.seed = p["seed"]
    rng = np.random.default_rng(p["seed"])
    os.makedirs(output_dir, exist_ok=True)

    a = sc.read_h5ad(input_path)
    a.layers["counts"] = a.X.copy()
    print("tumour cells:", a.n_obs, "| per mouse:", a.obs["sample_id"].value_counts().to_dict())

    # 1. tumour-only scVI latent
    h = a[:, ~a.var_names.isin(p["exclude_from_hvg"])].copy()
    sc.pp.highly_variable_genes(h, n_top_genes=p["n_hvg"], flavor="seurat_v3", layer="counts")
    h = h[:, h.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(h, layer="counts")
    model = scvi.model.SCVI(h, n_latent=p["n_latent"], n_layers=p["n_layers"])
    model.train(max_epochs=p["max_epochs"], early_stopping=True)
    a.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.neighbors(a, use_rep="X_scVI", n_neighbors=p["n_neighbors"], random_state=p["seed"])
    sc.tl.umap(a, random_state=p["seed"])

    # 2. resolution sweep + bootstrap stability
    rows = []
    for r in p["resolutions"]:
        full = leiden(a, r, p["seed"], f"leiden_{r}")
        sizes = pd.Series(full).value_counts()
        aris, jac = [], {c: [] for c in np.unique(full)}
        for b in range(p["n_boot"]):
            idx = np.sort(rng.choice(a.n_obs, int(p["boot_frac"] * a.n_obs), replace=False))
            s = a[idx].copy()
            sc.pp.neighbors(s, use_rep="X_scVI", n_neighbors=p["n_neighbors"], random_state=b)
            lab = leiden(s, r, b, "boot")
            aris.append(adjusted_rand_score(full[idx], lab))
            for c, j in per_cluster_jaccard(full[idx], lab).items():
                jac[c].append(j)
        med = {c: float(np.median(v)) for c, v in jac.items()}
        rows.append({"resolution": r, "k": len(sizes), "ARI_boot": float(np.mean(aris)),
                     "jac_mean": float(np.mean(list(med.values()))), "jac_min": float(min(med.values())),
                     "n_tiny_clusters": int((sizes < p["min_cluster_cells"]).sum())})
        print(rows[-1], flush=True)
    M = pd.DataFrame(rows)
    M["stable"] = M["jac_min"] >= p["jac_min"]
    M.to_csv(os.path.join(output_dir, "resolution_metrics.csv"), index=False)

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for i, col in enumerate(["k", "ARI_boot", "jac_min"]):
        ax[i].plot(M["resolution"], M[col], "o-"); ax[i].set_xlabel("Leiden resolution"); ax[i].set_ylabel(col)
    ax[2].axhline(p["jac_min"], ls="--", c="grey")
    fig.tight_layout(); fig.savefig(os.path.join(output_dir, "resolution_stability.png"), dpi=150); plt.close(fig)

    # 3. choose resolution
    stable = M[M["stable"] & (M["k"] > 1)]
    suggested = float(stable["resolution"].max()) if len(stable) else float(M.loc[M["jac_min"].idxmax(), "resolution"])
    res = p.get("res") or suggested
    print(f"\nsuggested resolution = {suggested}   | using = {res}")
    key = f"leiden_{res}"
    a.obs["cluster"] = a.obs[key].astype("category")

    # 4. markers + signature hints
    S = signature_scores(a, SIGNATURES)
    hint = S.groupby(a.obs["cluster"].values).mean().idxmax(axis=1)
    sc.tl.rank_genes_groups(a, "cluster", method="wilcoxon", use_raw=False)
    mk = sc.get.rank_genes_groups_df(a, None)
    mk = mk[(mk["pvals_adj"] < 0.05) & (mk["logfoldchanges"] > 0.5)]
    mk.to_csv(os.path.join(output_dir, f"markers_res{res}.csv"), index=False)
    top = mk.groupby("group").head(10).groupby("group")["names"].apply(", ".join)
    tab = pd.DataFrame({"n_cells": a.obs["cluster"].value_counts().sort_index(), "hint": hint, "top_markers": top})
    tab.to_csv(os.path.join(output_dir, f"cluster_summary_res{res}.csv"))
    print(tab.to_string())

    a.obs["state"] = [f"c{c}_{hint[c]}" for c in a.obs["cluster"].astype(str)]
    a.obs["exclude_from_tree"] = [hint[c] in p["exclude_hints"] for c in a.obs["cluster"].astype(str)]
    dot_genes = {k: [g for g in v if g in a.var_names][:4] for k, v in SIGNATURES.items()}
    sc.pl.dotplot(a, dot_genes, groupby="cluster", show=False, standard_scale="var")
    plt.savefig(os.path.join(output_dir, f"dotplot_res{res}.png"), dpi=150, bbox_inches="tight"); plt.close()
    for col in ["cluster", "state", "sample_id", "treatment"]:
        sc.pl.umap(a, color=col, show=False, title=f"scVI UMAP res {res}: {col}")
        plt.savefig(os.path.join(output_dir, f"UMAP_tumor_{col}.png"), dpi=150, bbox_inches="tight"); plt.close()

    # 5. per-cell table used by the tree step (16-bp barcode = join key with Cassiopeia cellBC)
    out = a.obs[["sample_id", "treatment", "cluster", "state", "exclude_from_tree"]].copy()
    out.insert(0, "cellBC", [c.split("-")[0] for c in a.obs_names])
    out.to_csv(os.path.join(output_dir, "cell_states.csv"), index_label="cell")
    a.write(os.path.join(output_dir, "tumor_states.h5ad"))
    return {"output_path": output_dir, "params_used": p,
            "summary": {"resolution": res, "k": int(a.obs["cluster"].nunique())}}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=os.path.join(os.environ.get("OUT_ST", "."), "tumor_cells.h5ad"))
    ap.add_argument("--outdir", default=os.environ.get("OUT_ST"))
    ap.add_argument("--res", type=float, default=None, help="force a Leiden resolution")
    ap.add_argument("--n_boot", type=int, default=DEFAULTS["n_boot"])
    a = ap.parse_args()
    run_pipeline(a.input, a.outdir, res=a.res, n_boot=a.n_boot)

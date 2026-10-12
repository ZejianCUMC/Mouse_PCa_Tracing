#!/usr/bin/env python
"""
Method 4 — Wilcoxon top-50 marker overlap (Jaccard) among the 12 states.
  Then Ward on 1 − Jaccard → dendrogram / heatmap.
"""
import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform

ROOT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
CACHE = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/_cache"
OUT = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/method4_jaccard"
os.makedirs(OUT, exist_ok=True)

STATE_ORDER = [
    "Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep", "Luminal-CR",
    "Luminal-Persistent", "Luminal-Squamous",
    "Basal-like",
    "Partial-EMT",
    "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF",
]
TOP_N = 50

def main():
    print("[m4] loading adata_12state…", flush=True)
    a = ad.read_h5ad(f"{CACHE}/adata_12state.h5ad")
    a.obs["putative_state"] = pd.Categorical(a.obs["putative_state"].astype(str), categories=STATE_ORDER)
    print("  shape:", a.shape)

    # ── Wilcoxon rank-sum markers ────────────────────────────────────────
    print("[m4] scanpy rank_genes_groups (wilcoxon)…", flush=True)
    sc.tl.rank_genes_groups(a, groupby="putative_state", method="wilcoxon", n_genes=TOP_N)

    names = pd.DataFrame(a.uns["rank_genes_groups"]["names"])
    scores = pd.DataFrame(a.uns["rank_genes_groups"]["scores"])
    names.to_csv(f"{OUT}/top{TOP_N}_markers_names.csv", index=False)
    scores.to_csv(f"{OUT}/top{TOP_N}_markers_scores.csv", index=False)

    # ── 12 × 12 Jaccard ──────────────────────────────────────────────────
    top_sets = {s: set(names[s].dropna().tolist()[:TOP_N]) for s in STATE_ORDER}
    n = len(STATE_ORDER)
    J = np.zeros((n, n), dtype=float)
    for i, si in enumerate(STATE_ORDER):
        for j, sj in enumerate(STATE_ORDER):
            u = top_sets[si] | top_sets[sj]
            J[i, j] = len(top_sets[si] & top_sets[sj]) / len(u) if u else np.nan
    J_df = pd.DataFrame(J, index=STATE_ORDER, columns=STATE_ORDER)
    J_df.to_csv(f"{OUT}/jaccard_top{TOP_N}.csv")

    # ── Ward on 1 − Jaccard ──────────────────────────────────────────────
    D = 1.0 - J
    np.fill_diagonal(D, 0.0)
    # symmetrize numerical
    D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), method="ward")
    order = leaves_list(Z)
    J_ord = J_df.iloc[order, order]

    # cluster memberships at k=2..5
    memb = {}
    for k in (2, 3, 4, 5):
        lab = fcluster(Z, t=k, criterion="maxclust")
        memb[k] = pd.Series(lab, index=STATE_ORDER, name=f"k{k}")
    df_memb = pd.concat(memb.values(), axis=1)
    df_memb.to_csv(f"{OUT}/state_membership_k2_k5.csv")
    print("[m4] membership at k=2..5:")
    print(df_memb)

    # ── heatmap ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 6.4), dpi=150)
    im = ax.imshow(J_ord.values, cmap="magma", vmin=0, vmax=max(0.3, np.nanmax(J - np.eye(n))))
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(J_ord.columns, rotation=90, fontsize=8)
    ax.set_yticklabels(J_ord.index, fontsize=8)
    for i in range(n):
        for j in range(n):
            v = J_ord.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color="white" if v < 0.5 else "black", fontsize=6)
    cbar = plt.colorbar(im, ax=ax, shrink=0.7)
    cbar.set_label(f"Jaccard(top-{TOP_N} markers)")
    ax.set_title(f"Method 4 — top-{TOP_N} Wilcoxon marker overlap (Ward-ordered)")
    plt.tight_layout()
    plt.savefig(f"{OUT}/jaccard_heatmap.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT}/jaccard_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # persist linkage
    np.savetxt(f"{OUT}/linkage_matrix.tsv", Z, delimiter="\t", header="a\tb\tdist\tsize", comments="")

    def collapse(k):
        d = df_memb[f"k{k}"]
        return {int(c): sorted(d.index[d == c].tolist()) for c in sorted(d.unique())}
    verdict = {f"k{k}": collapse(k) for k in (2, 3, 4, 5)}
    with open(f"{OUT}/verdict.json", "w") as f:
        json.dump(verdict, f, indent=2)
    print("[m4] verdict.json written")
    print("[m4] DONE")

if __name__ == "__main__":
    main()

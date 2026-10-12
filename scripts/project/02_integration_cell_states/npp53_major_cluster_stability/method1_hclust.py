#!/usr/bin/env python
"""
Method 1 — Hierarchical clustering on per-cluster mean log-norm expression.
  • 12 states × 2000 top-variance HVGs (variance computed across the 12 state means)
  • z-score per gene across the 12 states
  • Ward linkage on Euclidean distances
  • Cut at k = 2, 3, 4, 5 → member tables
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram

ROOT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
CACHE = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/_cache"
OUT = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/method1_hclust"
os.makedirs(OUT, exist_ok=True)

STATE_ORDER = [
    "Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep", "Luminal-CR",
    "Luminal-Persistent", "Luminal-Squamous",
    "Basal-like",
    "Partial-EMT",
    "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF",
]

def main():
    print("[m1] loading state means…", flush=True)
    sm = pd.read_parquet(f"{CACHE}/state_mean_lognorm.parquet").reindex(STATE_ORDER)
    print("  state_mean shape:", sm.shape)

    # ── select 2000 HVGs by across-state variance ────────────────────────
    var_across = sm.var(axis=0)
    top_hvg = var_across.nlargest(2000).index
    M = sm[top_hvg].copy()   # 12 × 2000
    # z-score per gene across the 12 states
    M_z = (M - M.mean(axis=0)) / (M.std(axis=0).replace(0, 1))
    print("  HVG-restricted, z-scored matrix:", M_z.shape)
    M_z.to_parquet(f"{OUT}/state_mean_hvg2000_z.parquet")

    # ── Ward linkage on Euclidean ────────────────────────────────────────
    Z = linkage(M_z.values, method="ward", metric="euclidean")
    np.savetxt(f"{OUT}/linkage_matrix.tsv", Z, delimiter="\t", header="a\tb\tdist\tsize", comments="")

    # cluster memberships at k = 2..5
    memb = {}
    for k in (2, 3, 4, 5):
        lab = fcluster(Z, t=k, criterion="maxclust")
        memb[k] = pd.Series(lab, index=STATE_ORDER, name=f"k{k}")
    df_memb = pd.concat(memb.values(), axis=1)
    df_memb.to_csv(f"{OUT}/state_membership_k2_k5.csv")
    print("[m1] membership at k=2..5:")
    print(df_memb)

    # ── dendrogram with cut heights annotated ────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5), dpi=130)
    ddg = dendrogram(
        Z, labels=STATE_ORDER, leaf_rotation=90, leaf_font_size=9,
        color_threshold=0, above_threshold_color="#333",
    )
    # cut lines for k=2..5 (between consecutive merge heights)
    heights = sorted(Z[:, 2])
    colors = {2: "#d62728", 3: "#ff7f0e", 4: "#2ca02c", 5: "#1f77b4"}
    for k, c in colors.items():
        # cut just below the (n-k)-th merge (0-indexed: heights[len-k])
        idx = len(heights) - k
        if 0 <= idx < len(heights):
            h_cut = heights[idx] - 0.01
            ax.axhline(y=h_cut, color=c, lw=1.2, ls="--", alpha=0.75)
            ax.text(0.02, h_cut, f" k={k}", color=c, fontsize=9, va="bottom")
    ax.set_ylabel("Ward distance")
    ax.set_title("Method 1 — Ward on state-mean HVG expression (12 fine states)")
    plt.tight_layout()
    plt.savefig(f"{OUT}/dendrogram_ward_hvg2k.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT}/dendrogram_ward_hvg2k.pdf", bbox_inches="tight")
    plt.close(fig)

    # ── verdict json ─────────────────────────────────────────────────────
    def collapse(k):
        d = df_memb[f"k{k}"]
        return {int(c): sorted(d.index[d == c].tolist()) for c in sorted(d.unique())}
    verdict = {f"k{k}": collapse(k) for k in (2, 3, 4, 5)}
    with open(f"{OUT}/verdict.json", "w") as f:
        json.dump(verdict, f, indent=2)
    print("[m1] verdict.json written")
    print("[m1] DONE")

if __name__ == "__main__":
    main()

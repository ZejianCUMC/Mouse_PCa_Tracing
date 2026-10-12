#!/usr/bin/env python
"""
Combined 4-panel summary figure and integrated verdict.

Layout:
  top-left  : Method 1 dendrogram (Ward on HVG state means)
  top-right : Method 2 silhouette + gap curves
  bot-left  : Method 3 chooseR-style bootstrap stability curve
  bot-right : Method 4 top-50 marker Jaccard heatmap
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

ROOT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
BASE = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability"
OUT = BASE

STATE_ORDER = [
    "Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep", "Luminal-CR",
    "Luminal-Persistent", "Luminal-Squamous",
    "Basal-like",
    "Partial-EMT",
    "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF",
]

def load_json(p):
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)

def main():
    v1 = load_json(f"{BASE}/method1_hclust/verdict.json")
    v2 = load_json(f"{BASE}/method2_silhouette_gap/verdict.json")
    v3 = load_json(f"{BASE}/method3_chooser/verdict.json")
    v4 = load_json(f"{BASE}/method4_jaccard/verdict.json")
    v5 = load_json(f"{BASE}/method5_canonical/verdict.json")

    # ── figure ───────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(15, 11), dpi=140, constrained_layout=False)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32,
                            left=0.06, right=0.98, top=0.94, bottom=0.06)

    # ── Panel A: Method 1 dendrogram ────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    Zm1 = np.loadtxt(f"{BASE}/method1_hclust/linkage_matrix.tsv", skiprows=1)
    ddg = dendrogram(Zm1, labels=STATE_ORDER, leaf_rotation=90, leaf_font_size=8,
                     color_threshold=0, above_threshold_color="#333", ax=ax1)
    heights = sorted(Zm1[:, 2])
    for k, c in {2: "#d62728", 3: "#ff7f0e", 4: "#2ca02c", 5: "#1f77b4"}.items():
        idx = len(heights) - k
        if 0 <= idx < len(heights):
            h_cut = heights[idx] - 0.01
            ax1.axhline(y=h_cut, color=c, lw=1.2, ls="--", alpha=0.75)
            ax1.text(0.02, h_cut, f" k={k}", color=c, fontsize=8, va="bottom",
                     transform=ax1.get_yaxis_transform())
    ax1.set_ylabel("Ward distance")
    ax1.set_title("A · Method 1 — Ward on state-mean HVG expression", fontsize=11)

    # ── Panel B: Method 2 silhouette + gap ──────────────────────────────
    ax2a = fig.add_subplot(gs[0, 1])
    ks = sorted(int(k) for k in v2["silhouette"].keys())
    sils = [v2["silhouette"][str(k)] for k in ks]
    gap = [v2["gap"][str(k)] for k in ks]
    sd = [v2["gap_sd"][str(k)] for k in ks]
    ax2a.plot(ks, sils, "o-", color="#1f77b4", label="silhouette")
    ax2a.axvline(v2["k_silhouette"], ls="--", color="#1f77b4", alpha=0.6)
    ax2a.set_xlabel("k")
    ax2a.set_ylabel("silhouette", color="#1f77b4")
    ax2a.tick_params(axis="y", labelcolor="#1f77b4")

    ax2b = ax2a.twinx()
    ax2b.errorbar(ks, gap, yerr=sd, fmt="s-", color="#d62728", capsize=3, label="gap")
    ax2b.axvline(v2["k_gap"], ls=":", color="#d62728", alpha=0.7)
    ax2b.set_ylabel("gap (Tibshirani)", color="#d62728")
    ax2b.tick_params(axis="y", labelcolor="#d62728")
    ax2a.set_title(
        f"B · Method 2 — silhouette (peak k={v2['k_silhouette']}) + gap (opt k={v2['k_gap']})",
        fontsize=11,
    )

    # ── Panel C: Method 3 chooseR stability curve ────────────────────────
    ax3 = fig.add_subplot(gs[1, 0])
    if v3 is None or not os.path.exists(f"{BASE}/method3_chooser/stability_curve.csv"):
        ax3.text(0.5, 0.5, "Method 3 chooseR — PENDING / FAILED",
                 ha="center", va="center", fontsize=13, color="grey")
        ax3.set_axis_off()
    else:
        curve = pd.read_csv(f"{BASE}/method3_chooser/stability_curve.csv")
        ax3.errorbar(
            curve["k"], curve["median_cluster_stab"],
            yerr=[curve["median_cluster_stab"] - curve["min_cluster_stab"],
                  curve["max_cluster_stab"] - curve["median_cluster_stab"]],
            fmt="o-", color="#2ca02c", capsize=4,
        )
        k_best = v3["k_best"]
        ax3.axvline(k_best, ls="--", color="#2ca02c", alpha=0.7,
                    label=f"recommended k={k_best}")
        ax3.set_xlabel("k")
        ax3.set_ylabel("bootstrap co-cluster stability (median)")
        ax3.set_title(f"C · Method 3 — chooseR-style scVI bootstrap (recommended k={k_best})",
                      fontsize=11)
        ax3.legend(loc="best", fontsize=9)

    # ── Panel D: Method 4 Jaccard heatmap ────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])
    J = pd.read_csv(f"{BASE}/method4_jaccard/jaccard_top50.csv", index_col=0)
    # apply the Ward reordering used inside method 4
    Zm4 = np.loadtxt(f"{BASE}/method4_jaccard/linkage_matrix.tsv", skiprows=1)
    from scipy.cluster.hierarchy import leaves_list
    order = leaves_list(Zm4)
    J_ord = J.iloc[order, order]
    n = J.shape[0]
    Jvals = J.values.copy()
    np.fill_diagonal(Jvals, np.nan)
    vmax = np.nanmax(Jvals)
    im = ax4.imshow(J_ord.values, cmap="magma", vmin=0, vmax=max(0.3, vmax))
    ax4.set_xticks(range(n))
    ax4.set_yticks(range(n))
    ax4.set_xticklabels(J_ord.columns, rotation=90, fontsize=7)
    ax4.set_yticklabels(J_ord.index, fontsize=7)
    for i in range(n):
        for j in range(n):
            v = J_ord.values[i, j]
            ax4.text(j, i, f"{v:.2f}", ha="center", va="center",
                     color="white" if v < 0.5 else "black", fontsize=5)
    cbar = plt.colorbar(im, ax=ax4, shrink=0.7, fraction=0.05, pad=0.02)
    cbar.set_label("Jaccard(top-50 markers)")
    ax4.set_title("D · Method 4 — top-50 Wilcoxon marker overlap (Ward-ordered)", fontsize=11)

    fig.suptitle(
        "Major-cluster stability of the 12 Npp53 tumor states (5 independent methods)",
        fontsize=13, y=0.985,
    )
    plt.savefig(f"{OUT}/major_cluster_stability_summary.png",
                dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT}/major_cluster_stability_summary.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print("[summary] wrote major_cluster_stability_summary.png/.pdf")


if __name__ == "__main__":
    main()

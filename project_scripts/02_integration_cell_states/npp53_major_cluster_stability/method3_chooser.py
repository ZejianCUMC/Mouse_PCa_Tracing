#!/usr/bin/env python
"""
Method 3 — chooseR-style bootstrap consensus stability on the 30-dim scVI cell
embedding (18,162 cells). Fallback since chooseR isn't in rMAST.

For each k in [2..6]:
  1) K-means on ALL cells at k → reference labels y_ref
  2) N_BOOT subsamples of 80% cells, K-means at k → labels
  3) Realign each bootstrap's labels to y_ref by Hungarian on confusion
  4) Per cell, count fraction of bootstraps in which the cell landed in the
     SAME cluster as its reference (rowwise cocluster consistency)
  5) Silhouette of the resulting 18,162 × k soft co-cluster matrix
     (metric = 1 - Jaccard-of-cluster-indicator vectors), averaged as a
     single stability score (mean per-cluster silhouette).

Also produces the per-state chooseR cluster assignment at the recommended k
by majority-voting cells per fine state.
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.optimize import linear_sum_assignment

ROOT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
CACHE = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/_cache"
OUT = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/method3_chooser"
os.makedirs(OUT, exist_ok=True)

STATE_ORDER = [
    "Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep", "Luminal-CR",
    "Luminal-Persistent", "Luminal-Squamous",
    "Basal-like",
    "Partial-EMT",
    "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF",
]
K_RANGE = list(range(2, 7))
N_BOOT = 100
SUB_FRAC = 0.80
SEED = 42


def align_labels(y_new, y_ref, k):
    """Hungarian-align y_new to y_ref based on their confusion matrix."""
    C = np.zeros((k, k), dtype=int)
    for a, b in zip(y_ref, y_new):
        C[a, b] += 1
    row, col = linear_sum_assignment(-C)  # maximize
    mapping = {old: new for new, old in zip(row, col)}
    return np.asarray([mapping[y] for y in y_new])


def bootstrap_stability(X, k, n_boot=100, sub_frac=0.80, seed=42):
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    ref = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(X)
    y_ref = ref.labels_
    # M[i,c] = fraction of bootstraps in which cell i landed in aligned-cluster c
    M = np.zeros((n, k), dtype=np.float32)
    counts = np.zeros(n, dtype=np.int32)
    for b in range(n_boot):
        idx = rng.choice(n, size=int(sub_frac * n), replace=False)
        Xi = X[idx]
        km = KMeans(n_clusters=k, n_init=5, random_state=seed + b).fit(Xi)
        y_bi = km.labels_
        y_bi_aligned = align_labels(y_bi, y_ref[idx], k)
        for j, c in zip(idx, y_bi_aligned):
            M[j, c] += 1
            counts[j] += 1
    counts_safe = np.where(counts == 0, 1, counts)
    M_freq = M / counts_safe[:, None]
    # per-cell stability = fraction of bootstraps in which cell landed in its ref cluster
    per_cell_stab = M_freq[np.arange(n), y_ref]
    # Mean stability aggregated per reference cluster (chooseR-style)
    df = pd.DataFrame({"cluster": y_ref, "stab": per_cell_stab})
    per_cluster = df.groupby("cluster")["stab"].mean()
    # Overall = median of per-cluster means (chooseR reports median cluster silhouette)
    overall = float(per_cluster.median())
    return overall, per_cluster.values, y_ref, per_cell_stab


def main():
    print("[m3] loading cell scVI…", flush=True)
    df = pd.read_parquet(f"{CACHE}/cell_scvi.parquet")
    X = df.drop(columns=["_state"]).values.astype(np.float32)
    states = df["_state"].values
    print(f"  X shape: {X.shape}   states: {len(set(states))}")

    curve = []
    per_cluster_all = {}
    y_ref_by_k = {}
    per_cell_stab_by_k = {}
    for k in K_RANGE:
        print(f"  bootstrap k={k}  N_BOOT={N_BOOT}…", flush=True)
        overall, per_c, y_ref, per_cell = bootstrap_stability(
            X, k, n_boot=N_BOOT, sub_frac=SUB_FRAC, seed=SEED,
        )
        curve.append({"k": k,
                      "median_cluster_stab": overall,
                      "min_cluster_stab": float(np.min(per_c)),
                      "max_cluster_stab": float(np.max(per_c)),
                      "mean_cluster_stab": float(np.mean(per_c))})
        per_cluster_all[k] = per_c.tolist()
        y_ref_by_k[k] = y_ref
        per_cell_stab_by_k[k] = per_cell
        print(f"    k={k}  median={overall:.3f}  per-cluster={np.round(per_c,3).tolist()}")

    curve_df = pd.DataFrame(curve)
    curve_df.to_csv(f"{OUT}/stability_curve.csv", index=False)

    # chooseR rule: max median-cluster stability (higher = more reproducible)
    k_best = int(curve_df.loc[curve_df["median_cluster_stab"].idxmax(), "k"])
    print(f"[m3] recommended k = {k_best}")

    # ── per-state assignment at k_best (majority vote of cells) ──────────
    y_best = y_ref_by_k[k_best]
    per_cell_stab = per_cell_stab_by_k[k_best]
    cell_asg = pd.DataFrame({
        "state": states,
        "cluster": y_best,
        "per_cell_stability": per_cell_stab,
    })
    cell_asg.to_csv(f"{OUT}/per_cell_assignment_k{k_best}.csv.gz", index=False, compression="gzip")

    # majority-vote per state
    maj = (
        cell_asg.groupby(["state", "cluster"]).size()
        .unstack(fill_value=0)
        .reindex(STATE_ORDER)
    )
    maj.to_csv(f"{OUT}/state_cluster_counts_k{k_best}.csv")
    state_major = maj.idxmax(axis=1)
    state_purity = maj.max(axis=1) / maj.sum(axis=1)
    state_asg = pd.DataFrame({
        "state": STATE_ORDER,
        "chooser_cluster": state_major.values,
        "purity": state_purity.values,
    })
    state_asg.to_csv(f"{OUT}/state_assignment_k{k_best}.csv", index=False)
    print(f"[m3] per-state chooseR assignment at k={k_best}:")
    print(state_asg)

    # ── plot ─────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4.2), dpi=140)
    ax.errorbar(
        curve_df["k"], curve_df["median_cluster_stab"],
        yerr=[curve_df["median_cluster_stab"] - curve_df["min_cluster_stab"],
              curve_df["max_cluster_stab"] - curve_df["median_cluster_stab"]],
        fmt="o-", color="#2ca02c", capsize=4,
    )
    ax.axvline(k_best, ls="--", color="#2ca02c", alpha=0.7, label=f"recommended k={k_best}")
    ax.set_xlabel("k")
    ax.set_ylabel("Bootstrap co-cluster stability (per-cluster median)")
    ax.set_title(f"Method 3 — chooseR-style stability on scVI ({N_BOOT} boots, {int(SUB_FRAC*100)}% subsample)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{OUT}/stability_curve.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT}/stability_curve.pdf", bbox_inches="tight")
    plt.close(fig)

    verdict = {
        "curve": {str(r["k"]): r["median_cluster_stab"] for r in curve},
        "per_cluster_stab": {str(k): v for k, v in per_cluster_all.items()},
        "k_best": k_best,
        "n_boot": N_BOOT,
        "sub_frac": SUB_FRAC,
        "note": "chooseR fallback: bootstrap KMeans + Hungarian alignment on scVI; recommended k = argmax(median per-cluster reproducibility).",
    }
    with open(f"{OUT}/verdict.json", "w") as f:
        json.dump(verdict, f, indent=2)
    print("[m3] DONE")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
Method 2 — Silhouette + Gap statistic (Tibshirani 2001) at k = 2..8
on the 12 × 30 state-mean scVI matrix (K-means).
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

ROOT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
CACHE = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/_cache"
OUT = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/method2_silhouette_gap"
os.makedirs(OUT, exist_ok=True)

STATE_ORDER = [
    "Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep", "Luminal-CR",
    "Luminal-Persistent", "Luminal-Squamous",
    "Basal-like",
    "Partial-EMT",
    "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF",
]

K_RANGE = list(range(2, 9))
N_BOOT = 100
SEED = 0

def kmeans_inertia(X, k, seed):
    km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(X)
    return km.inertia_, km.labels_

def gap_stat(X, k_range, n_boot=100, seed=0):
    """Tibshirani 2001 gap statistic — uniform-box reference."""
    rng = np.random.default_rng(seed)
    mins = X.min(axis=0)
    maxs = X.max(axis=0)
    log_wk = []
    log_wk_ref_mean = []
    log_wk_ref_sd = []
    for k in k_range:
        wk_obs, _ = kmeans_inertia(X, k, seed)
        log_wk.append(np.log(wk_obs))
        refs = []
        for b in range(n_boot):
            X_ref = rng.uniform(mins, maxs, size=X.shape)
            wk_ref, _ = kmeans_inertia(X_ref, k, seed + b)
            refs.append(np.log(wk_ref))
        refs = np.asarray(refs)
        log_wk_ref_mean.append(refs.mean())
        log_wk_ref_sd.append(refs.std() * np.sqrt(1 + 1.0 / n_boot))
    log_wk = np.asarray(log_wk)
    log_wk_ref_mean = np.asarray(log_wk_ref_mean)
    log_wk_ref_sd = np.asarray(log_wk_ref_sd)
    gap = log_wk_ref_mean - log_wk
    return gap, log_wk_ref_sd

def main():
    print("[m2] loading state-mean scVI…", flush=True)
    sm = pd.read_parquet(f"{CACHE}/state_mean_scvi.parquet").reindex(STATE_ORDER)
    X = sm.values
    print("  state_mean_scvi shape:", X.shape)

    sils, inertias = [], []
    labels_by_k = {}
    for k in K_RANGE:
        w, lab = kmeans_inertia(X, k, SEED)
        inertias.append(w)
        # silhouette needs ≥2 clusters and n_samples>k
        if k < X.shape[0]:
            s = silhouette_score(X, lab, metric="euclidean")
        else:
            s = np.nan
        sils.append(s)
        labels_by_k[k] = lab
        print(f"  k={k}  inertia={w:.3f}  silhouette={s:.3f}")

    gap, sd = gap_stat(X, K_RANGE, n_boot=N_BOOT, seed=SEED)
    print("[m2] gap stat:")
    for k, g, s_ in zip(K_RANGE, gap, sd):
        print(f"  k={k}  gap={g:.3f}  sd={s_:.3f}")

    # Tibshirani rule: smallest k such that gap(k) ≥ gap(k+1) − sd(k+1)
    k_gap = None
    for i, k in enumerate(K_RANGE[:-1]):
        if gap[i] >= gap[i + 1] - sd[i + 1]:
            k_gap = k
            break
    if k_gap is None:
        k_gap = int(K_RANGE[int(np.argmax(gap))])
        rule = "argmax(gap) fallback"
    else:
        rule = "Tibshirani smallest-k rule"
    k_sil = int(K_RANGE[int(np.argmax(sils))])
    print(f"[m2] optimal k by silhouette = {k_sil}")
    print(f"[m2] optimal k by gap        = {k_gap}  ({rule})")

    # ── membership tables at k_sil and k_gap ─────────────────────────────
    memb_rows = {}
    for k in sorted({k_sil, k_gap, 2, 3, 4, 5}):
        if k in labels_by_k:
            memb_rows[f"k{k}"] = pd.Series(labels_by_k[k], index=STATE_ORDER)
    df_memb = pd.concat(memb_rows.values(), axis=1)
    df_memb.columns = list(memb_rows.keys())
    df_memb.to_csv(f"{OUT}/state_membership_kmeans.csv")

    # ── plot ─────────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), dpi=140)
    ax1.plot(K_RANGE, sils, "o-", color="#1f77b4")
    ax1.axvline(k_sil, color="#1f77b4", ls="--", alpha=0.6, label=f"peak k={k_sil}")
    ax1.set_xlabel("k")
    ax1.set_ylabel("Silhouette (euclidean)")
    ax1.set_title("K-means silhouette on state-mean scVI (12×30)")
    ax1.legend()

    ax2.errorbar(K_RANGE, gap, yerr=sd, fmt="o-", color="#d62728", capsize=3)
    ax2.axvline(k_gap, color="#d62728", ls="--", alpha=0.6, label=f"opt k={k_gap}")
    ax2.set_xlabel("k")
    ax2.set_ylabel("Gap statistic")
    ax2.set_title(f"Gap (Tibshirani 2001, {N_BOOT} boot)")
    ax2.legend()
    plt.tight_layout()
    plt.savefig(f"{OUT}/silhouette_gap_curves.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT}/silhouette_gap_curves.pdf", bbox_inches="tight")
    plt.close(fig)

    verdict = {
        "silhouette": {str(k): float(s) for k, s in zip(K_RANGE, sils)},
        "gap": {str(k): float(g) for k, g in zip(K_RANGE, gap)},
        "gap_sd": {str(k): float(s) for k, s in zip(K_RANGE, sd)},
        "k_silhouette": int(k_sil),
        "k_gap": int(k_gap),
        "gap_rule": rule,
    }
    with open(f"{OUT}/verdict.json", "w") as f:
        json.dump(verdict, f, indent=2)
    print("[m2] DONE")

if __name__ == "__main__":
    main()

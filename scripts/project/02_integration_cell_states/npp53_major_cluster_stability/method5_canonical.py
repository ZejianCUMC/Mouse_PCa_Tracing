#!/usr/bin/env python
"""
Method 5 — Canonical prostate lineage marker sanity check.

For each candidate macro-grouping (k = 2 Epi/Mes;  k = 3 Lum/Basal/Mes;
k = 4 Lum/Basal/PartialEMT/Mes;  plus whatever the other methods produced),
compute mean z-score (across the 12 states) of each canonical marker per
macro-cluster and check for cleanliness.
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
CACHE = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/_cache"
OUT = f"{ROOT}/analysis/evolution_features_v2/major_cluster_stability/method5_canonical"
os.makedirs(OUT, exist_ok=True)

STATE_ORDER = [
    "Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep", "Luminal-CR",
    "Luminal-Persistent", "Luminal-Squamous",
    "Basal-like",
    "Partial-EMT",
    "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF",
]

MARKERS = {
    "Luminal":       ["Krt8", "Krt18", "Foxa1", "Ar", "Nkx3-1", "Nkx3.1", "Cd24a"],
    "Basal":         ["Krt5", "Krt14", "Trp63", "Krt6a"],
    "PartialEMT":    ["Vim", "Snai1", "Snai2", "Zeb1", "Zeb2"],
    "Mesenchymal":   ["Vim", "Fn1", "S100a4", "Twist1", "Twist2", "Cdh2"],
}

# Candidate groupings — a priori biological hypotheses ± results from other methods.
CANDIDATES = {
    "k2_EpiMes": {
        "Epithelial": ["Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep",
                        "Luminal-CR", "Luminal-Persistent", "Luminal-Squamous",
                        "Basal-like", "Partial-EMT"],
        "Mesenchymal": ["Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF"],
    },
    "k2_EpiMes_PE_in_Mes": {
        "Epithelial": ["Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep",
                        "Luminal-CR", "Luminal-Persistent", "Luminal-Squamous",
                        "Basal-like"],
        "Mesenchymal_wPE": ["Partial-EMT", "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF"],
    },
    "k3_LumBasalMes": {
        "Luminal": ["Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep",
                    "Luminal-CR", "Luminal-Persistent", "Luminal-Squamous"],
        "Basal": ["Basal-like", "Partial-EMT"],
        "Mesenchymal": ["Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF"],
    },
    "k3_LumBasalMes_PE_in_Mes": {
        "Luminal": ["Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep",
                    "Luminal-CR", "Luminal-Persistent", "Luminal-Squamous"],
        "Basal": ["Basal-like"],
        "Mesenchymal_wPE": ["Partial-EMT", "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF"],
    },
    "k4_LumBasalPEMes": {
        "Luminal": ["Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep",
                    "Luminal-CR", "Luminal-Persistent", "Luminal-Squamous"],
        "Basal": ["Basal-like"],
        "PartialEMT": ["Partial-EMT"],
        "Mesenchymal": ["Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF"],
    },
    "k5_LumLumSqBasalPEMes": {
        "Luminal": ["Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep",
                    "Luminal-CR", "Luminal-Persistent"],
        "LumSquamous": ["Luminal-Squamous"],
        "Basal": ["Basal-like"],
        "PartialEMT": ["Partial-EMT"],
        "Mesenchymal": ["Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF"],
    },
}

def main():
    print("[m5] loading state-mean log-norm…", flush=True)
    sm = pd.read_parquet(f"{CACHE}/state_mean_lognorm.parquet").reindex(STATE_ORDER)

    # z-score each gene across the 12 states
    z_all = (sm - sm.mean(axis=0)) / sm.std(axis=0).replace(0, 1)

    # Report each canonical marker's presence in var_names, and per-state z
    marker_summary_rows = []
    for cat, glist in MARKERS.items():
        for g in glist:
            present = g in z_all.columns
            marker_summary_rows.append({
                "category": cat, "gene": g, "present_in_var": present,
                "mean_z_luminal_states": z_all[g][["Luminal-Distal","Luminal-Proximal","Luminal-ARdep","Luminal-CR","Luminal-Persistent","Luminal-Squamous"]].mean() if present else np.nan,
                "z_basal_like": z_all[g]["Basal-like"] if present else np.nan,
                "z_partial_emt": z_all[g]["Partial-EMT"] if present else np.nan,
                "mean_z_mes_states": z_all[g][["Mes-Quiescent","Mes-Collagen","Mes-Stemness","Mes-myCAF"]].mean() if present else np.nan,
            })
    ms = pd.DataFrame(marker_summary_rows)
    ms.to_csv(f"{OUT}/marker_presence_and_state_zscore.csv", index=False)
    print("[m5] marker presence per state z-score written")

    # Build per-candidate table: macro-cluster × canonical-panel mean-z
    panel_summaries = {}
    for cname, grouping in CANDIDATES.items():
        rows = []
        for cat, glist in MARKERS.items():
            present = [g for g in glist if g in z_all.columns]
            if not present:
                continue
            for macro, states in grouping.items():
                sub_z = z_all.loc[states, present].mean(axis=1).mean()
                rows.append({"macro_cluster": macro, "marker_panel": cat, "mean_z": sub_z})
        df = pd.DataFrame(rows).pivot(index="macro_cluster", columns="marker_panel", values="mean_z")
        df = df.reindex(columns=[c for c in MARKERS.keys() if c in df.columns])
        panel_summaries[cname] = df
        df.to_csv(f"{OUT}/panel_zscore__{cname}.csv")
        print(f"\n[m5] {cname}")
        print(df.round(2))

    # cleanliness score per candidate = mean(diagonal-ish specificity)
    # for each canonical panel we designate an ideal high macro-cluster
    IDEAL = {
        "Luminal": ["Luminal", "LumSquamous", "Epithelial", "Epithelial_wPE"],
        "Basal": ["Basal"],
        "PartialEMT": ["PartialEMT", "Mesenchymal_wPE", "Basal"],  # if PE merged, still ok
        "Mesenchymal": ["Mesenchymal", "Mesenchymal_wPE"],
    }
    scores = {}
    for cname, df in panel_summaries.items():
        s = 0.0
        n = 0
        for panel in df.columns:
            ideal = [c for c in IDEAL[panel] if c in df.index]
            not_ideal = [c for c in df.index if c not in ideal]
            if ideal and not_ideal:
                gap = df.loc[ideal, panel].max() - df.loc[not_ideal, panel].max()
                s += gap
                n += 1
        scores[cname] = float(s / n) if n else np.nan
    scores_df = pd.Series(scores, name="cleanliness_gap").sort_values(ascending=False)
    scores_df.to_csv(f"{OUT}/cleanliness_scores.csv")
    print("\n[m5] cleanliness gap (higher = cleaner separation):")
    print(scores_df)

    # heatmap panel — grid of candidate summaries
    fig, axes = plt.subplots(2, 3, figsize=(13, 6.5), dpi=140)
    for ax, (cname, df) in zip(axes.flat, panel_summaries.items()):
        im = ax.imshow(df.values, cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")
        ax.set_xticks(range(df.shape[1]))
        ax.set_xticklabels(df.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(df.shape[0]))
        ax.set_yticklabels(df.index, fontsize=8)
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                v = df.values[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            color="white" if abs(v) > 0.8 else "black", fontsize=7)
        ax.set_title(f"{cname}\ngap={scores.get(cname, np.nan):.2f}", fontsize=9)
    for ax in axes.flat[len(panel_summaries):]:
        ax.axis("off")
    fig.suptitle("Method 5 — canonical-marker mean z-score by candidate grouping", y=1.02)
    plt.tight_layout()
    plt.savefig(f"{OUT}/candidate_panel_summary.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT}/candidate_panel_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    verdict = {
        "cleanliness_gap": scores,
        "best_candidate_by_gap": scores_df.index[0] if len(scores_df) else None,
    }
    with open(f"{OUT}/verdict.json", "w") as f:
        json.dump(verdict, f, indent=2)
    print("[m5] DONE")

if __name__ == "__main__":
    main()

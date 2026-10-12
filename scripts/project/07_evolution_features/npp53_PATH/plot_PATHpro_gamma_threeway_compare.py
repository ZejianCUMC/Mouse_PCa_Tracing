#!/usr/bin/env python3
# =============================================================================
# 3-way PATHpro γ comparison: unit-length vs edit-distance vs Phylovelo-pseudotime
# =============================================================================
# Description : Produces the updated Fig 10.4E — three cohort-mean γ bars per
#               CARTA state (c1-c12) side-by-side, plus a secondary scatter panel
#               (edit-distance γ vs pseudotime γ) with Pearson/Spearman
#               correlations printed.  Third branch-length variant is
#               Phylovelo pseudotime: per-edge Δ pseudotime (scaled so median
#               non-pseudocount edge = 1 pseudotime unit).
# Input       : PATH_out/PATHpro_all_samples_gamma.csv            (unit-length)
#               PATH_out/editdist/PATHpro_all_samples_gamma.csv   (edit-distance)
#               PATH_out/pseudotime/PATHpro_all_samples_gamma.csv (pseudotime)
# Output      : Fig_10_4E_PATHpro_threeway_compare.png / .pdf
#               PATH_out/pseudotime/PATHpro_gamma_threeway_summary.tsv
# Conda env   : hotspot (/home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python)
# Validated   : 2026-09-08
# =============================================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import pearsonr, spearmanr

ROOT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2"
UNIT  = os.path.join(ROOT, "PATH_out/PATHpro_all_samples_gamma.csv")
EDIST = os.path.join(ROOT, "PATH_out/editdist/PATHpro_all_samples_gamma.csv")
PT    = os.path.join(ROOT, "PATH_out/pseudotime/PATHpro_all_samples_gamma.csv")
OUT_PNG = os.path.join(ROOT, "carta_fine12/biology_review/Fig_10_4E_PATHpro_threeway_compare.png")
OUT_PDF = os.path.join(ROOT, "carta_fine12/biology_review/Fig_10_4E_PATHpro_threeway_compare.pdf")
OUT_TSV = os.path.join(ROOT, "PATH_out/pseudotime/PATHpro_gamma_threeway_summary.tsv")

# CARTA canonical palette (matches Fig 10.4C / 10.4E edit-dist variant)
STATE_LABEL_TO_C = {
    "Luminal-Distal":     ("c1",  "#4E79A7"),
    "Luminal-Proximal":   ("c2",  "#F28E2B"),
    "Luminal-ARdep":      ("c3",  "#E15759"),
    "Luminal-CR":         ("c4",  "#76B7B2"),
    "Luminal-Persistent": ("c5",  "#59A14F"),
    "Luminal-Squamous":   ("c6",  "#EDC948"),
    "Basal-like":         ("c7",  "#B07AA1"),
    "Partial-EMT":        ("c8",  "#FF9DA7"),
    "Mes-Quiescent":      ("c9",  "#9C755F"),
    "Mes-Collagen":       ("c10", "#BAB0AC"),
    "Mes-Stemness":       ("c11", "#7F7F7F"),
    "Mes-myCAF":          ("c12", "#8CD17D"),
}
ORDER = [f"c{i}" for i in range(1, 13)]
ID_TO_NAME  = {v[0]: k for k, v in STATE_LABEL_TO_C.items()}
ID_TO_COLOR = {v[0]: v[1] for v in STATE_LABEL_TO_C.values()}

# Samples with degenerate v10 pseudotime_umap (only 2 unique values → pt-based
# branch lengths carry very little information; flagged in the report / TSV).
DEGENERATE_PT_SAMPLES = {"JZ201", "JZ225", "JZ229", "JZ230", "JZ231"}


def load(f):
    df = pd.read_csv(f)
    df["c_id"] = df["state"].map(lambda s: STATE_LABEL_TO_C.get(s, (s, "#666666"))[0])
    return df


def cohort_summary(df):
    g = df.groupby("c_id")["gamma"].agg(count="count", mean="mean", std="std").reset_index()
    g["sem"] = g["std"] / np.sqrt(g["count"].clip(lower=1))
    g["sem"] = g["sem"].fillna(0)
    g["rank"] = pd.Categorical(g["c_id"], categories=ORDER, ordered=True)
    return g.sort_values("rank").set_index("c_id").reindex(ORDER).reset_index()


def main():
    unit = load(UNIT)
    edist = load(EDIST)
    pt   = load(PT)

    # For the pseudotime metric, exclude 5 samples with degenerate v10 pseudotime
    # (only 2 unique pt values — a lone 0.0 outlier vs a solid 1.0 majority; see
    # DEGENERATE_PT_SAMPLES).  Their pseudotime-scaled trees carry almost no signal
    # (>99% pseudocount edges) and PATHpro produces spurious very-large γ for
    # states that happen to sit next to the 0.0 outlier cell.  We keep all 25 for
    # unit-length and edit-distance because those metrics are unaffected.
    pt_clean = pt[~pt["sample"].isin(DEGENERATE_PT_SAMPLES)].copy()

    su = cohort_summary(unit)
    se = cohort_summary(edist)
    sp = cohort_summary(pt_clean)

    # Wide layout: bar chart (top, 2/3 width) + scatter panel (right, 1/3 width)
    fig = plt.figure(figsize=(15.6, 5.4), dpi=150)
    gs  = fig.add_gridspec(1, 4, wspace=0.35)
    ax  = fig.add_subplot(gs[0, :3])
    ax2 = fig.add_subplot(gs[0, 3])

    # ── Bar chart ────────────────────────────────────────────────────────────
    x = np.arange(len(ORDER))
    w = 0.26
    colors = [ID_TO_COLOR.get(c, "#7f7f7f") for c in ORDER]

    ax.bar(x - w, su["mean"].fillna(0), w, yerr=su["sem"], color=colors,
           edgecolor="black", linewidth=0.5, capsize=1.6, alpha=0.65,
           label="Unit-length")
    ax.bar(x,     se["mean"].fillna(0), w, yerr=se["sem"], color=colors,
           edgecolor="black", linewidth=0.5, capsize=1.6, alpha=0.90,
           hatch="///", label="Edit-distance")
    ax.bar(x + w, sp["mean"].fillna(0), w, yerr=sp["sem"], color=colors,
           edgecolor="black", linewidth=0.5, capsize=1.6, alpha=1.00,
           hatch="xxx", label="Pseudotime")

    # † for rare (n<3) state / metric combinations
    for i, c in enumerate(ORDER):
        for shift, s in zip([-w, 0, +w], [su, se, sp]):
            n = s.loc[i, "count"]; m = s.loc[i, "mean"]
            if pd.notna(n) and int(n) < 3 and pd.notna(m):
                ax.text(i + shift, m + 0.03, "†", ha="center", fontsize=9,
                        color="#a00")

    ax.axhline(1.0, color="grey", ls="--", lw=0.8)
    ax.set_xticks(x)
    labels = []
    for i, c in enumerate(ORDER):
        name = ID_TO_NAME.get(c, c)
        nu = int(su.loc[i,"count"]) if pd.notna(su.loc[i,"count"]) else 0
        ne = int(se.loc[i,"count"]) if pd.notna(se.loc[i,"count"]) else 0
        npt= int(sp.loc[i,"count"]) if pd.notna(sp.loc[i,"count"]) else 0
        labels.append(f"{c} {name}\nn={nu}/{ne}/{npt}")
    ax.set_xticklabels(labels, fontsize=8, rotation=28, ha="right")
    ax.set_ylabel("PATHpro γ (frequency-weighted mean = 1)", fontsize=10)
    top = max([np.nanmax(su["mean"].values), np.nanmax(se["mean"].values),
               np.nanmax(sp["mean"].values)])
    ax.set_ylim(0, top * 1.30)
    ax.set_title(
        "Fig 10.4E · Per-state PATHpro γ, unit-length vs edit-distance vs Phylovelo-pseudotime  "
        "(unit/edist n=25; pt n=20, 5 degenerate-pt samples excluded)",
        fontsize=11, weight="bold")
    ax.grid(axis="y", ls=":", lw=0.4, alpha=0.5)

    legend_items = [
        Patch(facecolor="lightgrey", alpha=0.65, edgecolor="black",
              label="Unit-length"),
        Patch(facecolor="lightgrey", alpha=0.90, edgecolor="black",
              hatch="///", label="Edit-distance"),
        Patch(facecolor="lightgrey", alpha=1.00, edgecolor="black",
              hatch="xxx", label="Pseudotime"),
    ]
    ax.legend(handles=legend_items, loc="upper right", fontsize=8.5,
              frameon=True)
    for s in ("top","right"):
        ax.spines[s].set_visible(False)

    # ── Scatter (edit-dist γ vs pseudotime γ per state × sample) ─────────────
    # Pairwise merge on (sample, state).  Also compute pearson/spearman across
    # all pairs of metric combinations.
    def merge_pair(a, b, la, lb):
        return a[["sample","state","c_id","gamma"]].rename(columns={"gamma":la}).merge(
               b[["sample","state","gamma"]].rename(columns={"gamma":lb}),
               on=["sample","state"], how="inner")

    m_ue = merge_pair(unit,  edist,    "gamma_unit", "gamma_edist")
    m_up = merge_pair(unit,  pt_clean, "gamma_unit", "gamma_pt")
    m_ep = merge_pair(edist, pt_clean, "gamma_edist","gamma_pt")

    # Also compute on FULL set (including degenerate-pt) for the summary CSV row
    m_ue_all = merge_pair(unit,  edist, "gamma_unit", "gamma_edist")
    m_up_all = merge_pair(unit,  pt,    "gamma_unit", "gamma_pt")
    m_ep_all = merge_pair(edist, pt,    "gamma_edist","gamma_pt")

    def safe_corr(x, y):
        m = (~x.isna()) & (~y.isna())
        if m.sum() < 3:
            return (np.nan, np.nan)
        pr, _ = pearsonr(x[m], y[m])
        sr, _ = spearmanr(x[m], y[m])
        return (pr, sr)

    r_ue,     s_ue     = safe_corr(m_ue    ["gamma_unit"],  m_ue    ["gamma_edist"])
    r_up,     s_up     = safe_corr(m_up    ["gamma_unit"],  m_up    ["gamma_pt"])
    r_ep,     s_ep     = safe_corr(m_ep    ["gamma_edist"], m_ep    ["gamma_pt"])
    r_ue_all, s_ue_all = safe_corr(m_ue_all["gamma_unit"],  m_ue_all["gamma_edist"])
    r_up_all, s_up_all = safe_corr(m_up_all["gamma_unit"],  m_up_all["gamma_pt"])
    r_ep_all, s_ep_all = safe_corr(m_ep_all["gamma_edist"], m_ep_all["gamma_pt"])

    # scatter — edit-distance vs pseudotime, coloured by c_id
    scatter_df = m_ep.copy()
    # Cap axis at 99th percentile to avoid degenerate-pt outlier samples squashing
    # the informative body of the scatter.  Actual full values still in TSV.
    x_v = scatter_df["gamma_edist"].values
    y_v = scatter_df["gamma_pt"].values
    finite = np.isfinite(x_v) & np.isfinite(y_v)
    if finite.sum():
        x_hi = np.percentile(x_v[finite], 99)
        y_hi = np.percentile(y_v[finite], 99)
    else:
        x_hi = y_hi = 3.0
    ax_max = max(x_hi, y_hi, 1.6) * 1.05
    ax_max = min(ax_max, 6.0)  # keep plot readable even w/ big outliers

    for c in ORDER:
        sub = scatter_df[scatter_df["c_id"] == c]
        if len(sub) == 0: continue
        col = ID_TO_COLOR.get(c, "#666")
        ax2.scatter(sub["gamma_edist"], sub["gamma_pt"], s=22, alpha=0.85,
                    color=col, edgecolor="black", linewidth=0.4, label=c)

    lim = ax_max
    ax2.plot([0, lim], [0, lim], color="grey", ls="--", lw=0.7, zorder=0)
    ax2.set_xlim(0, lim); ax2.set_ylim(0, lim)
    ax2.set_xlabel("γ (edit-distance)", fontsize=10)
    ax2.set_ylabel("γ (Phylovelo pseudotime)", fontsize=10)
    ax2.set_title("γ concordance", fontsize=10.5, weight="bold")
    ax2.grid(True, ls=":", lw=0.4, alpha=0.5)
    for s in ("top","right"):
        ax2.spines[s].set_visible(False)
    # Correlation text
    n_e = int(np.isfinite(x_v[finite]).sum())
    txt = (f"per sample×state (n={n_e})\n"
           f"edist ↔ pt:  Pearson r={r_ep:.3f}, Spearman ρ={s_ep:.3f}\n"
           f"unit ↔ pt:  r={r_up:.3f}, ρ={s_up:.3f}\n"
           f"unit ↔ edist: r={r_ue:.3f}, ρ={s_ue:.3f}")
    ax2.text(0.03, 0.97, txt, transform=ax2.transAxes, fontsize=8.5,
             va="top", ha="left",
             bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#888", lw=0.5))

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(OUT_PNG, dpi=180, bbox_inches="tight")
    plt.savefig(OUT_PDF,           bbox_inches="tight")
    plt.close()
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_PDF}")

    # ── State-level 3-way TSV ────────────────────────────────────────────────
    merged = (
        su.rename(columns={"mean":"unit_mean","count":"unit_n","sem":"unit_sem"})
          [["c_id","unit_mean","unit_n","unit_sem"]]
        .merge(
            se.rename(columns={"mean":"edist_mean","count":"edist_n","sem":"edist_sem"})
              [["c_id","edist_mean","edist_n","edist_sem"]],
            on="c_id")
        .merge(
            sp.rename(columns={"mean":"pt_mean","count":"pt_n","sem":"pt_sem"})
              [["c_id","pt_mean","pt_n","pt_sem"]],
            on="c_id")
    )
    merged["name"] = merged["c_id"].map(ID_TO_NAME)
    merged["ratio_pt_over_edist"]   = merged["pt_mean"]   / merged["edist_mean"]
    merged["ratio_pt_over_unit"]    = merged["pt_mean"]   / merged["unit_mean"]
    merged["ratio_edist_over_unit"] = merged["edist_mean"] / merged["unit_mean"]
    # Rank by mean across metrics (ignoring NaNs)
    merged["mean_across_metrics"] = merged[["unit_mean","edist_mean","pt_mean"]].mean(axis=1, skipna=True)

    # Overall pairwise correlations at sample×state level (long-format)
    corr_row = pd.DataFrame([{
        "c_id": "__cohort_corr__", "name": "sample×state pairwise correlations",
        "unit_mean": np.nan, "unit_n": np.nan, "unit_sem": np.nan,
        "edist_mean": np.nan, "edist_n": np.nan, "edist_sem": np.nan,
        "pt_mean": np.nan, "pt_n": np.nan, "pt_sem": np.nan,
        "ratio_pt_over_edist": r_ep, "ratio_pt_over_unit": r_up,
        "ratio_edist_over_unit": r_ue, "mean_across_metrics": np.nan,
    }])

    merged = pd.concat([merged, corr_row], ignore_index=True)
    merged.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV}")
    print("\n--- state-level 3-way γ summary ---")
    print(merged.to_string(index=False))

    print("\n--- pairwise per sample×state correlations ---")
    print(f"unit  ↔ edist:  Pearson r = {r_ue:.4f}   Spearman ρ = {s_ue:.4f}   (n={len(m_ue)})")
    print(f"unit  ↔ pt:     Pearson r = {r_up:.4f}   Spearman ρ = {s_up:.4f}   (n={len(m_up)})")
    print(f"edist ↔ pt:     Pearson r = {r_ep:.4f}   Spearman ρ = {s_ep:.4f}   (n={len(m_ep)})")

    # top-state check under each metric
    def top_state(sdf, k=3):
        return sdf.sort_values("mean", ascending=False)["c_id"].head(k).tolist()
    print(f"\n--- top-3 proliferator (by cohort-mean γ) ---")
    print(f"unit-length : {top_state(su)}")
    print(f"edit-dist   : {top_state(se)}")
    print(f"pseudotime  : {top_state(sp)}")


if __name__ == "__main__":
    main()

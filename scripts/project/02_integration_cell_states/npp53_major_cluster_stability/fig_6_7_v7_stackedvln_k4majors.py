#!/usr/bin/env python
"""
Fig 6.7 v7 — Seurat-style stacked violin plot.

  Rows    : 12 fine tumor states (c1..c12), Cycling + Ionocyte excluded.
  Cols    : combined marker set —
              (a) top-5 rank_genes_groups markers per k=4 major cell type
                  (data-driven block on the LEFT)
              (b) canonical prostate-lineage markers per fine state
                  (canonical block on the RIGHT, separated by a thin divider)
  Row bg  : very-light tint of that state's k=4 major color.
  Fill    : k=4 major color per state.
  Left    : color square (k=4 major) + `c[id]_<state>` label per row.
  X-axis  : gene names rotated + colored underline bar per gene (k=4 color).
  NO title, NO suptitle (per Zejian).

Output (do NOT overwrite existing fig_6_7_k4_major_markers.*):
  fig_6_7_v7_stackedvln_k4majors.{png,pdf}
"""
from pathlib import Path
import re
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

# ── Paths & constants ────────────────────────────────────────────────────────
BASE = Path(
    "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/"
    "1-mouse_based/Npp53_tracing"
)
OUT = BASE / "analysis/evolution_features_v2/major_cluster_stability"
CACHE = OUT / "_cache/adata_12state.h5ad"

FINE_TO_K4 = {
    "Luminal-Distal": "DistalAR_Lum", "Luminal-ARdep": "DistalAR_Lum", "Luminal-CR": "DistalAR_Lum",
    "Luminal-Proximal": "ProxBasalSq", "Luminal-Persistent": "ProxBasalSq",
    "Luminal-Squamous": "ProxBasalSq", "Basal-like": "ProxBasalSq",
    "Mes-Collagen": "CoreMes", "Mes-Stemness": "CoreMes", "Mes-myCAF": "CoreMes",
    "Partial-EMT": "EMT_Bridge", "Mes-Quiescent": "EMT_Bridge",
}
K4_ORDER = ["DistalAR_Lum", "ProxBasalSq", "EMT_Bridge", "CoreMes"]
K4_COLORS = {
    "DistalAR_Lum": "#1f5fa8", "ProxBasalSq": "#00a19b",
    "EMT_Bridge": "#e69f00", "CoreMes": "#d13845",
}

# Canonical 1-based order c1..c12 (per Zejian's spec)
FINE_ORDER = [
    "Luminal-Distal", "Luminal-Proximal", "Luminal-ARdep", "Luminal-CR",
    "Luminal-Persistent", "Luminal-Squamous", "Basal-like",
    "Partial-EMT", "Mes-Quiescent", "Mes-Collagen", "Mes-Stemness", "Mes-myCAF",
]
FINE_LABEL = {fs: f"c{i+1}_{fs}" for i, fs in enumerate(FINE_ORDER)}

CANONICAL = {
    "Luminal-Distal":    ["Tgm4", "Pbsn", "Krt8"],
    "Luminal-Proximal":  ["Msmb", "Sbp", "Aqp3"],
    "Luminal-ARdep":      ["Ar", "Fkbp5", "Nkx3-1", "Klk3"],
    "Luminal-CR":        ["Nectin4", "Elf3", "Klk4"],
    "Luminal-Persistent": ["Krt15", "Krt6a", "Ppp1r1b"],
    "Luminal-Squamous":  ["Krt5", "Sfn", "Krt6a"],
    "Basal-like":        ["Krt14", "Trp63", "Snai2"],
    "Partial-EMT":       ["Vim", "Zeb1", "Fn1", "Nes"],
    "Mes-Quiescent":     ["Igfbp2", "Tagln", "Cnn1"],
    "Mes-Collagen":      ["Col1a1", "Col1a2", "Col3a1"],
    "Mes-Stemness":      ["Sox2", "Cd44"],
    "Mes-myCAF":         ["Fap", "Postn", "Acta2", "Pdgfra"],
}

EXCL = re.compile(
    r"^(mt-|Rp[sl]|Malat1$|Neat1$|Gm[0-9]|Xist$|.*Cas9.*|.*mCherry.*|.*GFP.*|.*Cre$|.*Tomato.*)"
)


def hex_with_alpha(hex_color, alpha):
    """Return matplotlib RGBA tuple from hex + alpha (0..1)."""
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
    return (r, g, b, alpha)


def main():
    print("[fig6.7v7] Loading cached 12-state adata...")
    a = sc.read_h5ad(CACHE)
    print(f"          shape={a.shape}")

    # ── Order categoricals ────────────────────────────────────────────────
    a.obs["putative_state"] = (
        a.obs["putative_state"].astype("category").cat.set_categories(FINE_ORDER, ordered=True)
    )
    a.obs["k4_major"] = (
        a.obs["putative_state"].map(FINE_TO_K4).astype("category").cat.set_categories(K4_ORDER, ordered=True)
    )

    # ── Log-norm from counts (fresh, ensures rank_genes_groups + violin match) ──
    print("[fig6.7v7] Log-normalising counts layer...")
    a_m = sc.AnnData(X=a.layers["counts"].copy(), obs=a.obs.copy(), var=a.var.copy())
    sc.pp.normalize_total(a_m, target_sum=1e4)
    sc.pp.log1p(a_m)

    # ── Data-driven top-5 markers per k4 major ────────────────────────────
    print("[fig6.7v7] rank_genes_groups per k4_major (Wilcoxon)...")
    sc.tl.rank_genes_groups(a_m, groupby="k4_major", method="wilcoxon", n_genes=200)
    TOP_N = 5
    top_by_k4 = {}
    for grp in K4_ORDER:
        names = a_m.uns["rank_genes_groups"]["names"][grp]
        picked = []
        for g in names:
            if EXCL.match(g):
                continue
            if g in picked:
                continue
            picked.append(g)
            if len(picked) >= TOP_N:
                break
        top_by_k4[grp] = picked
        print(f"          {grp}: {picked}")

    # ── Assemble gene x-axis: data-driven → canonical (dedup, keep first) ──
    dd_genes, dd_group = [], []
    seen = set()
    for k4 in K4_ORDER:
        for g in top_by_k4[k4]:
            if g in seen or g not in a_m.var_names:
                continue
            dd_genes.append(g); dd_group.append(k4); seen.add(g)

    can_genes, can_group, can_fine = [], [], []
    dropped_canonical = []
    for fs in FINE_ORDER:
        for g in CANONICAL.get(fs, []):
            if g not in a_m.var_names:
                dropped_canonical.append((fs, g))
                continue
            if g in seen:
                continue  # already used in data-driven block or earlier canonical
            can_genes.append(g); can_group.append(FINE_TO_K4[fs]); can_fine.append(fs); seen.add(g)

    all_genes = dd_genes + can_genes
    all_group = dd_group + can_group
    n_dd, n_can, n_g = len(dd_genes), len(can_genes), len(dd_genes) + len(can_genes)
    print(f"[fig6.7v7] gene axis: {n_g} total  (data-driven={n_dd}, canonical={n_can})")
    print(f"          dropped canonical (absent in var_names): {dropped_canonical}")

    # ── Extract expression matrix (log-norm) ──────────────────────────────
    print("[fig6.7v7] extracting expression matrix...")
    X = a_m[:, all_genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    expr = pd.DataFrame(X, columns=all_genes, index=a_m.obs.index)
    expr["_state"] = a_m.obs["putative_state"].values

    # per-row (state) max used for y-lim (99.5% to trim wild outliers)
    per_state_max = expr.groupby("_state", observed=False)[all_genes].quantile(0.995).reindex(FINE_ORDER)
    per_state_max["_row"] = per_state_max.max(axis=1)
    print("          per-state row-max (log-norm):")
    print(per_state_max["_row"].round(2).to_string())

    # ── Figure ────────────────────────────────────────────────────────────
    NROW = len(FINE_ORDER)
    # Layout (all in inches, measured from bottom of figure):
    #   [0.00 .. 0.05]                 : bottom pad
    #   [0.05 .. 0.35]                 : k4 legend (LEGEND_H_IN=0.30")
    #   [0.35 .. 1.20]                 : rotated gene tick labels (TICK_H_IN=0.85")
    #   [1.20 .. 1.30]                 : colored underline strip (UL_H_IN=0.10")
    #   [1.30 .. 1.32]                 : UL_GAP_IN=0.02"
    #   [1.32 .. 1.32 + PANEL_H]       : 12 stacked violin panels
    #   [top .. top + TOP_MAR_IN=0.15"]: top pad
    LEFT_MAR_IN  = 2.40   # for k4 color-square + c[id]_ label (longest: c11_Mes-Stemness)
    PER_GENE_IN  = 0.26
    RIGHT_MAR_IN = 0.35   # for sec.axis-style top-value tick
    TOP_MAR_IN   = 0.15   # NO title
    ROW_H_IN     = 0.60   # per fine-state row
    UL_H_IN      = 0.10   # underline bar strip height
    UL_GAP_IN    = 0.02   # tiny gap between panels and underline strip
    TICK_H_IN    = 0.85   # rotated gene names (~10 chars @ fontsize 8, 90°)
    LEGEND_H_IN  = 0.30   # k4 legend at bottom
    LEG_PAD_IN   = 0.05   # bottom pad below legend
    PANEL_H      = NROW * ROW_H_IN
    BOT_MAR_IN   = LEG_PAD_IN + LEGEND_H_IN + TICK_H_IN + UL_H_IN + UL_GAP_IN

    fig_w = LEFT_MAR_IN + n_g * PER_GENE_IN + RIGHT_MAR_IN
    fig_h = TOP_MAR_IN + PANEL_H + BOT_MAR_IN
    print(f"[fig6.7v7] figure {fig_w:.2f}\" x {fig_h:.2f}\"  ({NROW}rows x {n_g}cols)")

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Fractional coordinates
    l_frac    = LEFT_MAR_IN / fig_w
    r_frac    = (fig_w - RIGHT_MAR_IN) / fig_w
    ul_b_in   = LEG_PAD_IN + LEGEND_H_IN + TICK_H_IN
    ul_t_in   = ul_b_in + UL_H_IN
    b_in      = ul_t_in + UL_GAP_IN
    t_in      = b_in + PANEL_H
    ul_b_frac = ul_b_in / fig_h
    ul_t_frac = ul_t_in / fig_h
    b_frac    = b_in / fig_h
    t_frac    = t_in / fig_h

    gs = fig.add_gridspec(
        NROW, 1,
        left=l_frac, right=r_frac, bottom=b_frac, top=t_frac,
        hspace=0,
    )
    axes = [fig.add_subplot(gs[i, 0]) for i in range(NROW)]

    # ── Left color-square + label + violin per row ────────────────────────
    for i, fs in enumerate(FINE_ORDER):
        ax = axes[i]
        k4 = FINE_TO_K4[fs]
        color = K4_COLORS[k4]
        sub = expr[expr["_state"] == fs]
        n_cells = len(sub)

        # Draw violins
        for j, g in enumerate(all_genes):
            vals = sub[g].values
            if len(vals) < 3 or vals.max() == 0 or np.unique(vals).size < 2:
                continue  # KDE undefined
            vp = ax.violinplot(
                vals, positions=[j], widths=0.85,
                showmeans=False, showmedians=False, showextrema=False,
            )
            for pc in vp["bodies"]:
                pc.set_facecolor(color)
                pc.set_edgecolor(color)
                pc.set_alpha(0.85)
                pc.set_linewidth(0.35)

        # Per-row y-limit — use 99.5% row max with a floor
        row_max = float(per_state_max.loc[fs, "_row"])
        if not np.isfinite(row_max) or row_max < 0.3:
            row_max = 0.5
        ax.set_ylim(-row_max * 0.03, row_max * 1.04)
        ax.set_xlim(-0.7, n_g - 0.3)

        # Very-light row-background tint of k4 color
        ax.set_facecolor(hex_with_alpha(color, 0.055))

        # Row separator lines — keep all four spines but thin them.
        # With hspace=0, adjacent rows' bottom-top spines overlap into a single
        # separator line; keeping them both visible gives clean edges on the
        # top of row 0 and the bottom of row 11 as well.
        for sp_name in ("top", "bottom", "left", "right"):
            ax.spines[sp_name].set_visible(True)
            ax.spines[sp_name].set_linewidth(0.6)
            ax.spines[sp_name].set_color("#333")

        # y-axis: no ticks on left; sec-axis tick on right showing "row_max"
        ax.set_yticks([])
        # Right side "sec.axis" style: draw a small tick + text at the top
        rlabel = f"{row_max:.1f}" if row_max >= 1 else f"{row_max:.2f}"
        ax.text(
            n_g - 0.35, row_max, rlabel,
            ha="right", va="top", fontsize=6.5, color="#555",
        )

        # No x-ticks on ANY panel row — gene labels are drawn on the underline
        # strip axis so they never overlap with the underline bars.
        ax.set_xticks([])

        # thin dashed divider between data-driven and canonical blocks
        ax.axvline(n_dd - 0.5, color="#333", lw=0.7, ls="--", alpha=0.65)

        # ── Left-of-panel: color square + c[id]_state label ────────────────
        # Use axis-coord annotation so it scales with the row.
        # Square first (a small filled Rectangle in axes-fraction coords)
        # positioned in the figure margin area.
        pos = ax.get_position()
        square_w_fig = 0.014
        square_x0_fig = l_frac - square_w_fig - 0.010
        square_h_fig = (pos.y1 - pos.y0) * 0.55
        square_y0_fig = pos.y0 + (pos.y1 - pos.y0 - square_h_fig) / 2.0
        rect = mpatches.Rectangle(
            (square_x0_fig, square_y0_fig), square_w_fig, square_h_fig,
            transform=fig.transFigure, facecolor=color, edgecolor="#222",
            linewidth=0.5, clip_on=False, zorder=5,
        )
        fig.patches.append(rect)
        # Label text: `c[id]_state`  (right-aligned so it stops before the square)
        label = FINE_LABEL[fs]
        text_x_fig = square_x0_fig - 0.006
        text_y_fig = (pos.y0 + pos.y1) / 2.0
        fig.text(
            text_x_fig, text_y_fig, label,
            ha="right", va="center", fontsize=9, color="#0f2747",
        )

    # ── Bottom underline strip axis: colored bar per gene + gene tick labels ──
    ax_ul = fig.add_axes(
        [l_frac, ul_b_frac, r_frac - l_frac, ul_t_frac - ul_b_frac]
    )
    ax_ul.set_xlim(-0.7, n_g - 0.3)
    ax_ul.set_ylim(0, 1)
    for sp in ax_ul.spines.values():
        sp.set_visible(False)
    ax_ul.set_yticks([])
    for j, g in enumerate(all_genes):
        color = K4_COLORS[all_group[j]]
        rect = mpatches.Rectangle(
            (j - 0.42, 0.10), 0.84, 0.80,
            facecolor=color, edgecolor="none",
        )
        ax_ul.add_patch(rect)
    # divider on the underline bar too
    ax_ul.axvline(n_dd - 0.5, color="#333", lw=0.8, ls="--", alpha=0.7)
    # Gene tick labels — drawn BELOW this axis (extends into the tick-area)
    ax_ul.set_xticks(range(n_g))
    ax_ul.set_xticklabels(all_genes, rotation=90, fontsize=8, fontstyle="italic")
    ax_ul.tick_params(axis="x", length=0, pad=2)

    # ── k4 legend at the very bottom ───────────────────────────────────────
    handles = [
        mpatches.Patch(
            facecolor=K4_COLORS[k], edgecolor="#222", linewidth=0.5,
            label=f"{i+1}. {k}",
        )
        for i, k in enumerate(K4_ORDER)
    ]
    # legend y anchor in fig coords: centre it inside the legend band
    legend_y_frac = (LEG_PAD_IN + LEGEND_H_IN / 2.0) / fig_h
    fig.legend(
        handles=handles, loc="center", ncol=4, frameon=False,
        fontsize=9, bbox_to_anchor=(0.5, legend_y_frac),
        handlelength=1.6, handleheight=1.0, columnspacing=1.6,
        borderaxespad=0,
    )

    # NO SUPTITLE — per Zejian.

    for fmt in ("png", "pdf"):
        out = OUT / f"fig_6_7_v7_stackedvln_k4majors.{fmt}"
        fig.savefig(out, dpi=200, bbox_inches=None)  # keep exact layout coordinates
        print(f"[fig6.7v7] saved: {out}")
    plt.close(fig)

    # ── Marker manifest ────────────────────────────────────────────────────
    marker_tab = pd.DataFrame({
        "col_index": range(n_g),
        "gene": all_genes,
        "block": ["data_driven"] * n_dd + ["canonical"] * n_can,
        "k4_group": all_group,
        "fine_state": [""] * n_dd + can_fine,
    })
    marker_tab.to_csv(OUT / "fig_6_7_v7_markers.csv", index=False)
    print(f"[fig6.7v7] saved: {OUT/'fig_6_7_v7_markers.csv'}")

    # dropped canonical list
    if dropped_canonical:
        pd.DataFrame(dropped_canonical, columns=["fine_state", "gene_missing"]).to_csv(
            OUT / "fig_6_7_v7_canonical_dropped.csv", index=False
        )
        print(f"[fig6.7v7] dropped canonical genes -> fig_6_7_v7_canonical_dropped.csv")

    print("[fig6.7v7] DONE.")


if __name__ == "__main__":
    main()

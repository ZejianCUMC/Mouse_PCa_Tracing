#!/usr/bin/env python3
"""
=============================================================================
PATHpro γ — per-state proliferation-rate aggregation + plots
=============================================================================
Description : Aggregate per-sample PATHpro γ (from PATHpro_all_samples_gamma.csv)
              into cohort-level per-state summaries + per-arm heatmap.
Input       : PATH_out/PATHpro_all_samples_gamma.csv (long: sample × state × γ)
              PATH_out/PATHpro_run_log.csv           (per-sample run status)
Output      : biology_review/npp53_PATHpro_gamma_per_state.{png,pdf,tsv}
              biology_review/npp53_PATHpro_gamma_per_arm_heatmap.{png,pdf,tsv}
              (also synced to Dropbox share/Npp53_tracing/evolution_features/)
Conda env   : rMAST (matplotlib + pandas + numpy)
=============================================================================
"""
from __future__ import annotations
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ── PATHS ──────────────────────────────────────────────────────────────
PROJ       = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing")
PATH_OUT   = PROJ / "analysis" / "evolution_features_v2" / "PATH_out"
GAMMA_CSV  = PATH_OUT / "PATHpro_all_samples_gamma.csv"
LOG_CSV    = PATH_OUT / "PATHpro_run_log.csv"
OUT_DIR    = PROJ / "analysis" / "evolution_features_v2" / "carta_fine12" / "biology_review"
DBX_DIR    = Path("/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/evolution_features")

# ── STATE PALETTE ──────────────────────────────────────────────────────
STATE_COLORS = {
    'Luminal-Distal':'#1f5fa8','Luminal-Proximal':'#9467bd','Luminal-ARdep':'#d23030',
    'Luminal-CR':'#f1903d','Luminal-Persistent':'#117c70','Luminal-Squamous':'#f78ec7',
    'Basal-like':'#8c564b','Partial-EMT':'#bcbd22','Mes-Quiescent':'#c49c94',
    'Mes-Collagen':'#5b5b5b','Mes-Stemness':'#8c1a1a','Mes-myCAF':'#4daf4a',
}
STATE_LABEL_SHORT = {
    'Luminal-Distal':'c1 Lum-Distal','Luminal-Proximal':'c2 Lum-Proximal',
    'Luminal-ARdep':'c3 Lum-ARdep','Luminal-CR':'c4 Lum-CR',
    'Luminal-Persistent':'c5 Lum-Persistent','Luminal-Squamous':'c6 Lum-Squamous',
    'Basal-like':'c7 Basal-like','Partial-EMT':'c8 Partial-EMT',
    'Mes-Quiescent':'c9 Mes-Quiescent','Mes-Collagen':'c10 Mes-Collagen',
    'Mes-Stemness':'c11 Mes-Stemness','Mes-myCAF':'c12 Mes-myCAF',
}
SHORT_INT = {s:i+1 for i,s in enumerate([
    'Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
    'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
    'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF'])}
ARM_ORDER  = ["Intact","Castration","Castration + Enzalutamide"]
ARM_LABEL  = {"Intact":"Intact","Castration":"Cas","Castration + Enzalutamide":"Cas+Enz"}

# ── Rendering knobs (tuned to avoid clipping / overlap) ────────────────
DPI       = 200
BAR_W_IN  = 9.5
BAR_H_IN  = 5.6
HM_W_IN   = 7.8
HM_H_IN   = 5.0

# ══════════════════════════════════════════════════════════════════════
# 1) LOAD + AGGREGATE
# ══════════════════════════════════════════════════════════════════════
def load_gamma() -> pd.DataFrame:
    if not GAMMA_CSV.exists():
        raise SystemExit(f"missing {GAMMA_CSV}; run npp53_PATHpro_run.R first")
    df = pd.read_csv(GAMMA_CSV)
    # normalise treatment column into the 3-arm buckets used elsewhere
    df["arm"] = df["treatment"].where(df["treatment"].isin(ARM_ORDER),
                                       df["treatment"])
    return df

def per_state_summary(gamma_df: pd.DataFrame) -> pd.DataFrame:
    """mean / median / SEM of γ across samples, per state (+ per-arm mean)."""
    agg = (gamma_df.groupby("state")["gamma"]
                   .agg(n_samples="count", mean="mean", median="median",
                        sd="std", min="min", max="max")
                   .reset_index())
    agg["sem"] = agg["sd"] / np.sqrt(agg["n_samples"].clip(lower=1))
    # per-arm mean, wide
    per_arm = (gamma_df.groupby(["state","arm"])["gamma"].mean()
                       .unstack("arm"))
    for arm in ARM_ORDER:
        col = f"mean_{ARM_LABEL[arm]}"
        agg[col] = agg["state"].map(per_arm.get(arm, pd.Series(dtype=float)))
    # short label + int ordering
    agg["short"]     = agg["state"].map(STATE_LABEL_SHORT)
    agg["state_int"] = agg["state"].map(SHORT_INT)
    return agg.sort_values("mean", ascending=False).reset_index(drop=True)

def per_arm_matrix(gamma_df: pd.DataFrame) -> pd.DataFrame:
    """states (canonical order) × arms → mean γ."""
    mat = (gamma_df.groupby(["state","arm"])["gamma"].mean()
                   .unstack("arm").reindex(columns=ARM_ORDER))
    # order rows canonically 1..12
    ordered_states = sorted(mat.index, key=lambda s: SHORT_INT.get(s, 99))
    return mat.reindex(ordered_states)

# ══════════════════════════════════════════════════════════════════════
# 2) PLOTS
# ══════════════════════════════════════════════════════════════════════
def plot_bar(summary: pd.DataFrame, out_stem: Path) -> None:
    """Horizontal bar chart: per-state mean γ, sorted descending.

    n<3 states get no error bar (SEM undefined) and are visually flagged.
    """
    df = summary.copy()  # already sorted
    df["sem_plot"] = df["sem"].fillna(0.0)              # avoid matplotlib warnings
    n  = len(df)
    fig, ax = plt.subplots(figsize=(BAR_W_IN, BAR_H_IN), dpi=DPI)
    y = np.arange(n)[::-1]                              # top = highest γ
    colors = [STATE_COLORS.get(s, "#999999") for s in df["state"]]
    edges  = ["#222" if r["n_samples"] >= 3 else "#c11" for _,r in df.iterrows()]
    lws    = [0.6 if r["n_samples"] >= 3 else 1.3      for _,r in df.iterrows()]
    ax.barh(y, df["mean"], xerr=df["sem_plot"],
            color=colors, edgecolor=edges, linewidth=lws,
            error_kw=dict(elinewidth=0.9, capsize=3, ecolor="#333"))
    # y ticks with short label + n
    labels = [f"{r['short']}  (n={int(r['n_samples'])})" for _,r in df.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("PATHpro γ (per-state proliferation rate, total_prolif=1 norm)",
                  fontsize=11)
    # single-line title (was 2 lines and collided with the γ=1 caption)
    ax.set_title("Fig 10.4C · Per-state mean γ across 25-sample Npp53 cohort",
                 fontsize=11, loc="left", pad=14)
    # normalization reference line
    ax.axvline(1.0, color="#666", lw=0.7, ls="--", alpha=0.7)
    xmax = float(np.nanmax(df["mean"] + df["sem_plot"]))
    if not np.isfinite(xmax) or xmax <= 0:
        xmax = float(df["mean"].max() or 1.0)
    # annotate mean value at bar's right edge (with n-tag when n<3)
    for yi, (_, r) in zip(y, df.iterrows()):
        v      = float(r["mean"])
        sem    = float(r["sem_plot"])
        tag    = f"{v:.2f}" if r["n_samples"] >= 3 else f"{v:.2f}†"
        color  = "#222"     if r["n_samples"] >= 3 else "#a00"
        ax.text(v + sem + 0.015*xmax, yi, tag,
                fontsize=9, va="center", ha="left", color=color)
    ax.set_xlim(0, xmax * 1.20)
    ax.set_ylim(-0.6, n - 0.4)
    ax.grid(axis="x", ls=":", color="#bbb", alpha=0.6)
    ax.set_axisbelow(True)
    for spine in ("top","right"): ax.spines[spine].set_visible(False)
    # "γ = 1" caption anchored TOP-RIGHT (axes coords), well past the left-anchored
    # title's endpoint, and with a leader offset so it hugs the corner
    ax.annotate("γ = 1  (freq-weighted normalization)",
                xy=(1.0, 1.0), xycoords="axes fraction",
                xytext=(-2, 4), textcoords="offset points",
                fontsize=8.5, color="#555", ha="right", va="bottom")
    # combined subtitle + † legend as a single figure-level footer BELOW the axes
    fig.text(0.5, 0.005,
             "error bars = SEM   ·   red-edged bars / † values: n<3 samples "
             "(unreliable — single-sample for n=1, 2-point SEM for n=2)",
             fontsize=8.5, color="#444", ha="center", va="bottom")
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(f"{out_stem}.png", dpi=DPI, bbox_inches="tight")
    plt.savefig(f"{out_stem}.pdf",              bbox_inches="tight")
    plt.close(fig)

def plot_heatmap(mat: pd.DataFrame, out_stem: Path) -> None:
    """per-arm γ heatmap: 3 arms (cols) × 12 states (rows), z-scored PER STATE."""
    m = mat.copy()
    # z-score per state (row-wise) — highlights arm shifts, not overall γ level
    row_mean = m.mean(axis=1)
    row_sd   = m.std(axis=1).replace(0, np.nan)
    z = (m.sub(row_mean, axis=0)).div(row_sd, axis=0)
    z = z.fillna(0.0)
    # keep only states present (rows with all-NaN mean drop out)
    keep = m.notna().any(axis=1)
    z = z.loc[keep]; m = m.loc[keep]
    n_rows = len(z); n_cols = len(z.columns)
    fig, ax = plt.subplots(figsize=(HM_W_IN, HM_H_IN), dpi=DPI)
    # diverging blue↔red
    cmap = LinearSegmentedColormap.from_list("bwr_soft",
             ["#1f5fa8","#eef2f7","#d23030"], N=256)
    vmax = max(abs(z.values.min()), abs(z.values.max()))
    vmax = max(vmax, 1e-6)
    im = ax.imshow(z.values, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    # ticks
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([ARM_LABEL[c] for c in z.columns], fontsize=11)
    row_labels = [STATE_LABEL_SHORT.get(s, s) for s in z.index]
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=10)
    # cell annotations: raw γ value ("—" if state not observed in that arm at ≥20 cells).
    # Use %.2g for values below 0.01 so numerically-zero (~1e-9, ~4e-3) do NOT all collapse to "0.00".
    def _fmt(v: float) -> str:
        av = abs(v)
        if av == 0:            return "0"
        if av < 0.01:          return f"{v:.2g}"      # 2 sig figs (matches title text)
        return f"{v:.2f}"
    for i in range(n_rows):
        for j in range(n_cols):
            v = m.iloc[i,j]
            zv = z.iloc[i,j]
            if pd.isna(v):
                ax.text(j, i, "—", ha="center", va="center",
                        fontsize=10, color="#777")
                continue
            tc = "#111" if abs(zv) < 0.7*vmax else "#f8f8f8"
            ax.text(j, i, _fmt(v), ha="center", va="center",
                    fontsize=8.8, color=tc)
    ax.set_title("Fig 10.4D · Per-arm mean γ heatmap\n"
                 "(cells show raw γ, %.2g when |γ|<0.01; color = row z-score; "
                 "'—' = state not observed at ≥20 cells in that arm)",
                 fontsize=11, loc="left")
    ax.set_xlabel("Treatment arm", fontsize=10)
    ax.set_ylabel("Cell state (CARTA c1-c12)", fontsize=10)
    # colorbar
    cb = fig.colorbar(im, ax=ax, fraction=0.036, pad=0.03)
    cb.set_label("Row z-score of γ", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    # thin grid between cells
    ax.set_xticks(np.arange(-.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="#ffffff", linestyle="-", linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)
    plt.tight_layout()
    plt.savefig(f"{out_stem}.png", dpi=DPI, bbox_inches="tight")
    plt.savefig(f"{out_stem}.pdf",              bbox_inches="tight")
    plt.close(fig)

# ══════════════════════════════════════════════════════════════════════
# 3) DELIVERY
# ══════════════════════════════════════════════════════════════════════
def sync_to_dropbox(files: list[Path]) -> None:
    DBX_DIR.mkdir(parents=True, exist_ok=True)
    for src in files:
        dst = DBX_DIR / src.name
        shutil.copy2(src, dst)
        print(f"[sync] {src.name} → {dst}")

def main():
    gamma_df = load_gamma()
    print(f"[load] gamma table: {len(gamma_df)} rows "
          f"({gamma_df['sample'].nunique()} samples, "
          f"{gamma_df['state'].nunique()} states)")
    summary = per_state_summary(gamma_df)
    per_arm = per_arm_matrix(gamma_df)
    # TSV outputs
    summary_tsv = OUT_DIR / "npp53_PATHpro_gamma_per_state.tsv"
    per_arm_tsv = OUT_DIR / "npp53_PATHpro_gamma_per_arm_heatmap.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)
    per_arm.to_csv(per_arm_tsv, sep="\t")
    # Plots
    bar_stem     = OUT_DIR / "npp53_PATHpro_gamma_per_state"
    heatmap_stem = OUT_DIR / "npp53_PATHpro_gamma_per_arm_heatmap"
    plot_bar(summary, bar_stem)
    plot_heatmap(per_arm, heatmap_stem)
    # Report top 3 (n>=3 samples only — otherwise SEM is undefined)
    print("\n=== TOP-3 states by cohort mean γ (n>=3 samples) ===")
    robust = summary[summary["n_samples"] >= 3].reset_index(drop=True)
    for _,r in robust.head(3).iterrows():
        print(f"  {r['state']:20s}  mean={r['mean']:.3f}  SEM={r['sem']:.3f}  "
              f"(n={int(r['n_samples'])})")
    print("\n=== ARM-SPECIFIC γ shifts (per-state Intact vs Cas vs Cas+Enz) ===")
    for _,r in robust.iterrows():
        vals = [r.get(f"mean_{ARM_LABEL[a]}", np.nan) for a in ARM_ORDER]
        show = "  ".join(
            (f"{ARM_LABEL[a]}={v:.2f}" if pd.notna(v) else f"{ARM_LABEL[a]}=NA")
            for a, v in zip(ARM_ORDER, vals))
        v_arr = np.array([v for v in vals if pd.notna(v)], dtype=float)
        if v_arr.size >= 2:
            spread = float(v_arr.max() - v_arr.min())
        else:
            spread = float("nan")
        print(f"  {r['state']:20s}  {show}  Δmax-min={spread:.2f}")
    # Sync deliverables
    to_sync = [
        Path(f"{bar_stem}.png"),     Path(f"{bar_stem}.pdf"),     summary_tsv,
        Path(f"{heatmap_stem}.png"), Path(f"{heatmap_stem}.pdf"), per_arm_tsv,
    ]
    sync_to_dropbox(to_sync)
    print("\nDone.")

if __name__ == "__main__":
    main()

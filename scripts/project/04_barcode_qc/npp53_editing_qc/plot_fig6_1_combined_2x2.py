#!/usr/bin/env python3
# =============================================================================
# Fig 6.1 — Combined 2×2 QC Panel
# =============================================================================
# Description : 2×2 figure combining:
#               (A) Per-sample editing rate (detected only) — 3 stacked panels
#               (B) Cas9 edits per cell — violin by treatment arm
#               (C) Allele uniqueness — grouped bar
#               (D) Cross-sample allele sharing — grouped bar
# Input       : per_library allele tables, per_cell_qc.csv,
#               all_samples_allele_table.csv, npp53_clean_49619_metadata.csv
# Output      : fig6_1_combined_qc.png/.pdf
# Conda env   : scvi
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy import stats
import os, shutil, time
from collections import defaultdict

# ── PATHS ────────────────────────────────────────────────────────────────────
PROJECT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
AT_DIR  = os.path.join(PROJECT, "lineage_tracing/allele_tables/per_library")
QC_DIR  = os.path.join(PROJECT, "QC")
AT_QC   = os.path.join(PROJECT, "lineage_tracing/allele_table_QC")
DROPBOX = "/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/Cassiopeia_QC"
STEM    = "fig6_1_combined_qc"

# ── CONSTANTS ────────────────────────────────────────────────────────────────
LIBRARIES = ["MJZ008", "MJZ019", "MJZ021", "MJZ022", "MJZ023", "MJZ024", "MJZ025"]

ARM_ORDER  = ["Intact", "Castration", "Castration + Enzalutamide"]
ARM_SHORT  = {"Intact": "Intact", "Castration": "Cas", "Castration + Enzalutamide": "Cas+Enz"}
ARM_COLORS = {"Intact": "#4878CF", "Castration": "#D65F5F", "Castration + Enzalutamide": "#6ACC65"}


# ═══════════════════════════════════════════════════════════════════════════════
# DATA — Panel A: per-sample editing rate (detected only)
# ═══════════════════════════════════════════════════════════════════════════════
print("Loading metadata …")
meta = pd.read_csv(os.path.join(QC_DIR, "npp53_clean_49619_metadata.csv"))
meta["library"]  = meta["cell"].str.split("_").str[0]
meta["bare_bc"]  = meta["cell"].str.split("_").str[1].str.replace("-1", "")

neo = meta[~meta["celltype"].isin(["Macrophage_Myeloid", "Endothelial"])]
bc_lookup = {}
for lib in LIBRARIES:
    m = neo[neo["library"] == lib]
    bc_lookup[lib] = {r["bare_bc"]: (r["sample_id"], r["treatment"])
                      for _, r in m.iterrows()}

print("Parsing allele tables for Panel A …")
t0 = time.time()
rows_data = []
for lib in LIBRARIES:
    bmap = bc_lookup.get(lib, {})
    if not bmap:
        continue
    atf = os.path.join(AT_DIR, f"{lib}_allele_table.csv")
    if not os.path.exists(atf):
        continue
    print(f"  {lib} …", end="", flush=True)
    at = pd.read_csv(atf)
    at_m = at[at["cellBC"].isin(set(bmap.keys()))]
    print(f" {len(at_m)} rows")
    for _, row in at_m.iterrows():
        sid, trt = bmap[row["cellBC"]]
        d = {"sample": sid, "treatment": trt}
        for site in ["r1", "r2", "r3"]:
            val = row[site]
            if pd.isna(val):
                d[f"{site}_state"] = "undetected"
            elif "[None]" in str(val):
                d[f"{site}_state"] = "uncut"
            else:
                d[f"{site}_state"] = "edited"
        rows_data.append(d)

df_edit = pd.DataFrame(rows_data)
print(f"  {len(df_edit)} rows in {time.time()-t0:.1f}s")

# Per-sample detected-only rate
edit_records = []
for (sid, trt), grp in df_edit.groupby(["sample", "treatment"]):
    rec = {"sample": sid, "treatment": trt}
    for site in ["r1", "r2", "r3"]:
        edited = (grp[f"{site}_state"] == "edited").sum()
        uncut  = (grp[f"{site}_state"] == "uncut").sum()
        detected = edited + uncut
        rec[f"{site}_rate"] = edited / detected * 100 if detected > 0 else np.nan
    edit_records.append(rec)

ss = pd.DataFrame(edit_records)
ss["arm_idx"] = ss["treatment"].map({a: i for i, a in enumerate(ARM_ORDER)})
ss = ss.sort_values(["arm_idx", "sample"]).reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════════════════════
# DATA — Panels B, C, D
# ═══════════════════════════════════════════════════════════════════════════════
print("\nLoading QC data for panels B-D …")
cell_qc = pd.read_csv(os.path.join(AT_QC, "per_cell_qc.csv"), low_memory=False)
cell_qc_sample = cell_qc[cell_qc["scope"] == "sample"].copy()
sample_treatment = dict(cell_qc_sample[["sample", "treatment"]].drop_duplicates().values)

allele_df = pd.read_csv(os.path.join(PROJECT, "lineage_tracing/allele_tables/all_samples_allele_table.csv"),
                         low_memory=False)
allele_sample = allele_df[allele_df["sample_id"].str.startswith("JZ")].copy()
unique_samples = sorted(allele_sample["sample_id"].unique())


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE — 2×2 layout
# ═══════════════════════════════════════════════════════════════════════════════
print("\nPlotting 2×2 …")

fig = plt.figure(figsize=(16, 12))
outer_gs = GridSpec(2, 2, hspace=0.32, wspace=0.28,
                    left=0.06, right=0.92, top=0.93, bottom=0.06)

BAR_WIDTH = 0.65

# ── PANEL A: Per-sample editing rate — 3 touching sub-panels ─────────────
inner_gs_a = GridSpecFromSubplotSpec(3, 1, subplot_spec=outer_gs[0, 0], hspace=0)

n_samples = len(ss)
x_a = np.arange(n_samples)

sites_info = [("r1_rate", "Barcode-1"),
              ("r2_rate", "Barcode-2"),
              ("r3_rate", "Barcode-3")]

ax_a_list = []
for si, (col, subtitle) in enumerate(sites_info):
    if si == 0:
        ax = fig.add_subplot(inner_gs_a[si])
    else:
        ax = fig.add_subplot(inner_gs_a[si], sharex=ax_a_list[0])
    ax_a_list.append(ax)

    colors = [ARM_COLORS[t] for t in ss["treatment"]]
    vals = ss[col].values
    ax.bar(x_a, vals, width=0.7, color=colors,
           edgecolor="white", linewidth=0.5, zorder=3)

    # y-axis: ticks at 2, 50, 98 labeled as 0, 50, 100
    ax.set_yticks([2, 50, 98])
    ax.set_yticklabels(["0", "50", "100"], fontsize=7)
    ax.set_ylim(-2, 107)

    # Right-side subtitle
    ax.text(1.02, 0.5, subtitle, transform=ax.transAxes,
            fontsize=8, fontweight="bold", va="center", ha="left",
            rotation=-90)

    # Frame: keep all 4 spines
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)

    # Arm separators
    for i in range(1, n_samples):
        if ss.iloc[i]["treatment"] != ss.iloc[i-1]["treatment"]:
            ax.axvline(i - 0.5, color="grey", linestyle=":", linewidth=0.8, alpha=0.5)

    # Hide x-tick labels for top two panels
    if si < 2:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis='x', length=0)

# Shared y-axis title (centered on panel A)
pos_top = ax_a_list[0].get_position()
pos_bot = ax_a_list[-1].get_position()
mid_y_a = (pos_top.y1 + pos_bot.y0) / 2
fig.text(0.02, mid_y_a, "Editing rate (%)", va="center", ha="center",
         rotation="vertical", fontsize=9, fontweight="bold")

# X labels on bottom sub-panel
ax_a_list[-1].set_xticks(x_a)
ax_a_list[-1].set_xticklabels(ss["sample"], rotation=45, ha="right", fontsize=6)

# Legend at top of Panel A
handles_a = [mpatches.Patch(facecolor=ARM_COLORS[a], label=ARM_SHORT[a]) for a in ARM_ORDER]
ax_a_list[0].legend(handles=handles_a, fontsize=6, loc="upper right", ncol=3, framealpha=0.9)

# Panel label "A"
ax_a_list[0].text(-0.08, 1.10, "A", transform=ax_a_list[0].transAxes,
                  fontsize=16, fontweight="bold", va="top")


# ── PANEL B: Violin — Cas9 edits per cell ──────────────────────────────────
ax_b = fig.add_subplot(outer_gs[0, 1])

violin_data = []
for arm in ARM_ORDER:
    vals = cell_qc_sample.loc[cell_qc_sample["treatment"] == arm, "n_edited_sites"].values
    violin_data.append(vals)

parts = ax_b.violinplot(violin_data, positions=range(len(ARM_ORDER)),
                        showmedians=False, showextrema=False, widths=0.7)
for pc, arm in zip(parts["bodies"], ARM_ORDER):
    pc.set_facecolor(ARM_COLORS[arm])
    pc.set_edgecolor("black")
    pc.set_linewidth(0.6)
    pc.set_alpha(0.75)

for i, vals in enumerate(violin_data):
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    ax_b.vlines(i, q1, q3, color="black", linewidth=2.5, zorder=3)
    ax_b.scatter(i, med, color="white", edgecolor="black", s=30, zorder=4, linewidth=1)

for i in range(len(ARM_ORDER)):
    ax_b.text(i, -1, f"{len(violin_data[i]):,}", ha="center", va="top", fontsize=6)

pairs = [(0, 1), (0, 2), (1, 2)]
y_top = max(v.max() for v in violin_data)
for pidx, (i, j) in enumerate(pairs):
    _, p = stats.mannwhitneyu(violin_data[i], violin_data[j], alternative="two-sided")
    y_bracket = y_top + 4 + pidx * 6
    ax_b.plot([i, i, j, j], [y_bracket - 1, y_bracket, y_bracket, y_bracket - 1],
              color="black", linewidth=0.8)
    p_str = "p < 1e-16" if p < 1e-16 else f"p = {p:.2e}"
    ax_b.text((i + j) / 2, y_bracket + 0.5, p_str, ha="center", va="bottom", fontsize=5)

ax_b.set_xticks(range(len(ARM_ORDER)))
ax_b.set_xticklabels([ARM_SHORT[a] for a in ARM_ORDER], fontsize=8)
ax_b.set_ylabel("Cas9 edits per cell", fontsize=9)
ax_b.set_ylim(bottom=-5)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)
ax_b.text(-0.08, 1.05, "B", transform=ax_b.transAxes,
          fontsize=16, fontweight="bold", va="top")


# ── PANEL C: Allele uniqueness ──────────────────────────────────────────────
ax_c = fig.add_subplot(outer_gs[1, 0])

CELL_BINS = [1, 2, 3, 4, "≥5"]
CELL_BIN_LABELS = ["1", "2", "3", "4", "≥5"]

sample_fracs_g = {}
for sid in unique_samples:
    sub = allele_sample[allele_sample["sample_id"] == sid]
    counts = sub.groupby(["intBC", "allele"])["cellBC"].nunique()
    total = len(counts)
    if total == 0:
        sample_fracs_g[sid] = [0.0] * len(CELL_BINS)
        continue
    fracs = []
    for b in CELL_BINS:
        n = (counts >= 5).sum() if b == "≥5" else (counts == b).sum()
        fracs.append(n / total)
    sample_fracs_g[sid] = fracs

n_bins = len(CELL_BINS)
n_arms = len(ARM_ORDER)
bar_width = 0.22
x_pos = np.arange(n_bins)

for arm_idx, arm in enumerate(ARM_ORDER):
    arm_samples = [s for s in unique_samples if sample_treatment.get(s) == arm]
    n_s = len(arm_samples)
    if n_s == 0:
        continue
    mat = np.array([sample_fracs_g[s] for s in arm_samples])
    means = mat.mean(axis=0)
    sems = mat.std(axis=0, ddof=1) / np.sqrt(n_s) if n_s > 1 else np.zeros(n_bins)
    off = (arm_idx - (n_arms - 1) / 2) * bar_width
    ax_c.bar(x_pos + off, means, width=bar_width,
             color=ARM_COLORS[arm], edgecolor="white", linewidth=0.5,
             label=f"{ARM_SHORT[arm]} (n={n_s})", zorder=3)
    ax_c.errorbar(x_pos + off, means, yerr=sems,
                  fmt="none", ecolor="black", capsize=2, linewidth=0.8, zorder=4)

ax_c.set_xticks(x_pos)
ax_c.set_xticklabels(CELL_BIN_LABELS, fontsize=9)
ax_c.set_xlabel("Cells per unique allele", fontsize=9)
ax_c.set_ylabel("Proportion", fontsize=9)
ax_c.set_ylim(0, 1.0)
ax_c.legend(fontsize=6, loc="upper right", title="Arm", title_fontsize=7,
            framealpha=0.9, edgecolor="grey")
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
ax_c.text(-0.08, 1.05, "C", transform=ax_c.transAxes,
          fontsize=16, fontweight="bold", va="top")


# ── PANEL D: Cross-sample allele sharing ────────────────────────────────────
ax_d = fig.add_subplot(outer_gs[1, 1])

SAMPLE_BINS = [1, 2, 3, 4, "≥5"]
SAMPLE_BIN_LABELS = ["1", "2", "3", "4", "≥5"]

allele_to_samples = defaultdict(set)
for _, row in allele_sample.iterrows():
    allele_to_samples[(row["intBC"], row["allele"])].add(row["sample_id"])

sample_fracs_h = {}
for sid in unique_samples:
    sub = allele_sample[allele_sample["sample_id"] == sid]
    sample_alleles = sub.groupby(["intBC", "allele"]).size().index.tolist()
    total = len(sample_alleles)
    if total == 0:
        sample_fracs_h[sid] = [0.0] * len(SAMPLE_BINS)
        continue
    sharing_counts = [len(allele_to_samples[k]) for k in sample_alleles]
    fracs = []
    for b in SAMPLE_BINS:
        n = sum(1 for c in sharing_counts if c >= 5) if b == "≥5" else sum(1 for c in sharing_counts if c == b)
        fracs.append(n / total)
    sample_fracs_h[sid] = fracs

n_bins_h = len(SAMPLE_BINS)

for arm_idx, arm in enumerate(ARM_ORDER):
    arm_samples = [s for s in unique_samples if sample_treatment.get(s) == arm]
    n_s = len(arm_samples)
    if n_s == 0:
        continue
    mat = np.array([sample_fracs_h[s] for s in arm_samples])
    means = mat.mean(axis=0)
    sems = mat.std(axis=0, ddof=1) / np.sqrt(n_s) if n_s > 1 else np.zeros(n_bins_h)
    off = (arm_idx - (n_arms - 1) / 2) * bar_width
    ax_d.bar(x_pos + off, means, width=bar_width,
             color=ARM_COLORS[arm], edgecolor="white", linewidth=0.5,
             label=f"{ARM_SHORT[arm]} (n={n_s})", zorder=3)
    ax_d.errorbar(x_pos + off, means, yerr=sems,
                  fmt="none", ecolor="black", capsize=2, linewidth=0.8, zorder=4)

ax_d.set_xticks(x_pos)
ax_d.set_xticklabels(SAMPLE_BIN_LABELS, fontsize=9)
ax_d.set_xlabel("Samples sharing allele", fontsize=9)
ax_d.set_ylabel("Proportion", fontsize=9)
ax_d.set_ylim(0, 1.0)
ax_d.legend(fontsize=6, loc="upper right", title="Arm", title_fontsize=7,
            framealpha=0.9, edgecolor="grey")
ax_d.spines["top"].set_visible(False)
ax_d.spines["right"].set_visible(False)
ax_d.text(-0.08, 1.05, "D", transform=ax_d.transAxes,
          fontsize=16, fontweight="bold", va="top")


# ── SAVE ─────────────────────────────────────────────────────────────────────
os.makedirs(QC_DIR, exist_ok=True)
out_png = os.path.join(QC_DIR, f"{STEM}.png")
out_pdf = os.path.join(QC_DIR, f"{STEM}.pdf")

fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
print(f"\nSaved: {out_png}")
print(f"Saved: {out_pdf}")

for ext in [".png", ".pdf"]:
    src = os.path.join(QC_DIR, f"{STEM}{ext}")
    dst = os.path.join(DROPBOX, f"{STEM}{ext}")
    shutil.copy2(src, dst)
    print(f"Synced: {dst}")

plt.close()
print("Done.")

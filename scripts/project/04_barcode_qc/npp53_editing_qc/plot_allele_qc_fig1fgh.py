#!/usr/bin/env python3
# =============================================================================
# Allele-table QC — Lu et al. 2024 Nature Fig 1f–h analogs
# =============================================================================
# Description : Three-panel QC figure for Npp53 Cas9 lineage tracing:
#               (f) Cas9 edits per cell — violin by treatment arm
#               (g) Allele uniqueness — grouped bar: x = #cells sharing allele,
#                   bars = treatment arm, y = proportion (mean ± SEM)
#               (h) Cross-sample sharing — grouped bar: x = #samples sharing
#                   allele, bars = treatment arm, y = proportion (mean ± SEM)
# Input       : per_cell_qc.csv, all_samples_allele_table.csv
# Output      : basic_qc_like_LuNat_Fig1f_to_1h.png/.pdf
# Conda env   : scvi
# Style ref   : Lu et al. 2024 Nature 636:233–240, Fig 1f–h
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from collections import defaultdict
from scipy import stats
import shutil, os

# ── PATHS ────────────────────────────────────────────────────────────────────
PROJECT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
QC_DIR  = os.path.join(PROJECT, "lineage_tracing/allele_table_QC")
OUT_DIR = os.path.join(PROJECT, "QC")
DROPBOX = "/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/Cassiopeia_QC"

STEM = "basic_qc_like_LuNat_Fig1f_to_1h"

# ── ARM ORDERING & COLORS ───────────────────────────────────────────────────
ARM_ORDER = ["Intact", "Castration", "Castration + Enzalutamide"]
ARM_SHORT = {"Intact": "Intact", "Castration": "Cas", "Castration + Enzalutamide": "Cas+Enz"}
ARM_COLORS = {"Intact": "#4878CF", "Castration": "#D65F5F", "Castration + Enzalutamide": "#6ACC65"}

# ── LOAD DATA ────────────────────────────────────────────────────────────────
print("Loading per_cell_qc …")
cell_qc = pd.read_csv(os.path.join(QC_DIR, "per_cell_qc.csv"), low_memory=False)
cell_qc_sample = cell_qc[cell_qc["scope"] == "sample"].copy()
print(f"  {len(cell_qc_sample)} sample-level cells")

# Build sample → treatment map
sample_treatment = cell_qc_sample[["sample", "treatment"]].drop_duplicates()
sample_treatment = dict(zip(sample_treatment["sample"], sample_treatment["treatment"]))

print("Loading allele table …")
allele_df = pd.read_csv(os.path.join(PROJECT, "lineage_tracing/allele_tables/all_samples_allele_table.csv"),
                        low_memory=False)
# Keep only sample-level entries (JZ*)
allele_sample = allele_df[allele_df["sample_id"].str.startswith("JZ")].copy()
print(f"  {len(allele_sample)} sample-level allele rows")

# ── FIGURE SETUP ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5.0), gridspec_kw={"width_ratios": [1.0, 1.2, 1.2]})
fig.suptitle("Npp53 Lineage Tracing — Allele Table QC  (cf. Lu et al. 2024 Fig 1f–h)",
             fontsize=11, fontweight="bold", y=0.98)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL f: Cas9 edits per cell — violin by treatment arm
# Lu style: violin + count annotation + median label + P-values
# ══════════════════════════════════════════════════════════════════════════════
ax = axes[0]
ax.set_title("", fontsize=13, fontweight="bold", loc="left", pad=8)

violin_data = []
for arm in ARM_ORDER:
    vals = cell_qc_sample.loc[cell_qc_sample["treatment"] == arm, "n_edited_sites"].values
    violin_data.append(vals)

parts = ax.violinplot(violin_data, positions=range(len(ARM_ORDER)),
                      showmedians=False, showextrema=False, widths=0.7)
for pc, arm in zip(parts["bodies"], ARM_ORDER):
    pc.set_facecolor(ARM_COLORS[arm])
    pc.set_edgecolor("black")
    pc.set_linewidth(0.6)
    pc.set_alpha(0.75)

# Overlay median + IQR box
for i, vals in enumerate(violin_data):
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    ax.vlines(i, q1, q3, color="black", linewidth=2.5, zorder=3)
    ax.scatter(i, med, color="white", edgecolor="black", s=30, zorder=4, linewidth=1)

# Sample counts at bottom (like Lu's count annotations)
ymin = 0
for i, arm in enumerate(ARM_ORDER):
    n = len(violin_data[i])
    ax.text(i, ymin - 1, f"{n:,}", ha="center", va="top", fontsize=7, color="black")

# Pairwise Wilcoxon rank-sum P-values (like Lu)
pairs = [(0, 1), (0, 2), (1, 2)]
y_top = max(v.max() for v in violin_data)
for pidx, (i, j) in enumerate(pairs):
    _, p = stats.mannwhitneyu(violin_data[i], violin_data[j], alternative="two-sided")
    # Bracket
    y_bracket = y_top + 4 + pidx * 6
    ax.plot([i, i, j, j], [y_bracket - 1, y_bracket, y_bracket, y_bracket - 1],
            color="black", linewidth=0.8)
    p_str = "p < 1e-16" if p < 1e-16 else f"p = {p:.2e}"
    ax.text((i + j) / 2, y_bracket + 0.5, p_str, ha="center", va="bottom", fontsize=6)

ax.set_xticks(range(len(ARM_ORDER)))
ax.set_xticklabels([ARM_SHORT[a] for a in ARM_ORDER], fontsize=9)
ax.set_ylabel("Number of Cas9 edits per cell", fontsize=10)
ax.set_ylim(bottom=ymin - 5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL g: Allele uniqueness — grouped bar chart
# Lu style: x = number of cells with unique barcode (1, 2, 3, 4, ≥5)
#           bars grouped by sample type (= treatment arm)
#           y = proportion, mean ± SEM across samples within each arm
# ══════════════════════════════════════════════════════════════════════════════
ax = axes[1]
ax.set_title("", fontsize=13, fontweight="bold", loc="left", pad=8)

CELL_BINS = [1, 2, 3, 4, "≥5"]
CELL_BIN_LABELS = ["1", "2", "3", "4", "≥5"]

print("Computing allele uniqueness per sample (grouped bar) …")

# Per sample: compute proportion of alleles in each cell-count bin
sample_fracs_g = {}  # sample_id -> list of fracs for each bin
unique_samples = sorted(allele_sample["sample_id"].unique())
for sid in unique_samples:
    sub = allele_sample[allele_sample["sample_id"] == sid]
    counts = sub.groupby(["intBC", "allele"])["cellBC"].nunique()
    total = len(counts)
    if total == 0:
        sample_fracs_g[sid] = [0.0] * len(CELL_BINS)
        continue
    fracs = []
    for b in CELL_BINS:
        if b == "≥5":
            n = (counts >= 5).sum()
        else:
            n = (counts == b).sum()
        fracs.append(n / total)
    sample_fracs_g[sid] = fracs

# Group by arm, compute mean ± SEM
n_bins = len(CELL_BINS)
n_arms = len(ARM_ORDER)
bar_width = 0.22
x_positions = np.arange(n_bins)

for arm_idx, arm in enumerate(ARM_ORDER):
    arm_samples = [s for s in unique_samples if sample_treatment.get(s) == arm]
    n_s = len(arm_samples)
    if n_s == 0:
        continue
    mat = np.array([sample_fracs_g[s] for s in arm_samples])
    means = mat.mean(axis=0)
    sems = mat.std(axis=0, ddof=1) / np.sqrt(n_s) if n_s > 1 else np.zeros(n_bins)

    offset = (arm_idx - (n_arms - 1) / 2) * bar_width
    bars = ax.bar(x_positions + offset, means, width=bar_width,
                  color=ARM_COLORS[arm], edgecolor="white", linewidth=0.5,
                  label=f"{ARM_SHORT[arm]} (n={n_s})", zorder=3)
    ax.errorbar(x_positions + offset, means, yerr=sems,
                fmt="none", ecolor="black", capsize=2, linewidth=0.8, zorder=4)

ax.set_xticks(x_positions)
ax.set_xticklabels(CELL_BIN_LABELS, fontsize=10)
ax.set_xlabel("Number of cells with\nunique allele", fontsize=10)
ax.set_ylabel("Proportion", fontsize=10)
ax.set_ylim(0, 1.0)
ax.legend(fontsize=7, loc="upper right", title="Sample type", title_fontsize=8,
          framealpha=0.9, edgecolor="grey")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL h: Cross-sample allele sharing — grouped bar chart
# Lu style: x = number of samples with unique barcode (1, 2, 3, 4, ≥5)
#           bars grouped by sample type (= treatment arm)
#           y = proportion, mean ± SEM across samples within each arm
# ══════════════════════════════════════════════════════════════════════════════
ax = axes[2]
ax.set_title("", fontsize=13, fontweight="bold", loc="left", pad=8)

SAMPLE_BINS = [1, 2, 3, 4, "≥5"]
SAMPLE_BIN_LABELS = ["1", "2", "3", "4", "≥5"]

print("Computing cross-sample allele sharing (grouped bar) …")

# First, build global map: (intBC, allele) -> set of sample_ids
allele_to_samples = defaultdict(set)
for _, row in allele_sample.iterrows():
    allele_to_samples[(row["intBC"], row["allele"])].add(row["sample_id"])

# Per sample: for each of its alleles, how many samples share that allele?
# Then bin the proportions.
sample_fracs_h = {}
for sid in unique_samples:
    sub = allele_sample[allele_sample["sample_id"] == sid]
    # Unique alleles in this sample
    sample_alleles = sub.groupby(["intBC", "allele"]).size().index.tolist()
    total = len(sample_alleles)
    if total == 0:
        sample_fracs_h[sid] = [0.0] * len(SAMPLE_BINS)
        continue
    # For each allele, how many samples share it?
    sharing_counts = [len(allele_to_samples[k]) for k in sample_alleles]
    fracs = []
    for b in SAMPLE_BINS:
        if b == "≥5":
            n = sum(1 for c in sharing_counts if c >= 5)
        else:
            n = sum(1 for c in sharing_counts if c == b)
        fracs.append(n / total)
    sample_fracs_h[sid] = fracs

# Group by arm, compute mean ± SEM
n_bins_h = len(SAMPLE_BINS)

for arm_idx, arm in enumerate(ARM_ORDER):
    arm_samples = [s for s in unique_samples if sample_treatment.get(s) == arm]
    n_s = len(arm_samples)
    if n_s == 0:
        continue
    mat = np.array([sample_fracs_h[s] for s in arm_samples])
    means = mat.mean(axis=0)
    sems = mat.std(axis=0, ddof=1) / np.sqrt(n_s) if n_s > 1 else np.zeros(n_bins_h)

    offset = (arm_idx - (n_arms - 1) / 2) * bar_width
    bars = ax.bar(x_positions + offset, means, width=bar_width,
                  color=ARM_COLORS[arm], edgecolor="white", linewidth=0.5,
                  label=f"{ARM_SHORT[arm]} (n={n_s})", zorder=3)
    ax.errorbar(x_positions + offset, means, yerr=sems,
                fmt="none", ecolor="black", capsize=2, linewidth=0.8, zorder=4)

ax.set_xticks(x_positions)
ax.set_xticklabels(SAMPLE_BIN_LABELS, fontsize=10)
ax.set_xlabel("Number of samples with\nunique allele", fontsize=10)
ax.set_ylabel("Proportion", fontsize=10)
ax.set_ylim(0, 1.0)
ax.legend(fontsize=7, loc="upper right", title="Sample type", title_fontsize=8,
          framealpha=0.9, edgecolor="grey")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── SAVE ─────────────────────────────────────────────────────────────────────
plt.tight_layout(rect=[0, 0.02, 1, 0.95])

out_png = os.path.join(OUT_DIR, f"{STEM}.png")
out_pdf = os.path.join(OUT_DIR, f"{STEM}.pdf")

fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")

# Copy to Dropbox
for ext in [".png", ".pdf"]:
    src = os.path.join(OUT_DIR, f"{STEM}{ext}")
    dst = os.path.join(DROPBOX, f"{STEM}{ext}")
    shutil.copy2(src, dst)
    print(f"Synced: {dst}")

plt.close()
print("Done.")

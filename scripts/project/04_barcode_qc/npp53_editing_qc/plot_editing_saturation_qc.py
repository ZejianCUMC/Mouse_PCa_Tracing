#!/usr/bin/env python3
# =============================================================================
# Cas9 Editing Saturation QC
# =============================================================================
# Description : Tests whether Cas9 editing has reached saturation across
#               the 3 barcode target sites (r1, r2, r3).
#               Panel A: Per-site editing rate by treatment arm
#               Panel B: Per-cell unedited site profile (0/1/2/3 sites edited)
#               Panel C: Editing rate by arm (box + strip) — per sample
#               Panel D: Allele diversity vs editing rate
# Input       : per_library allele tables, npp53_clean_49619_metadata.csv
# Output      : editing_saturation_qc.png/.pdf
# Conda env   : scvi
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
import re, os, shutil, time

# ── PATHS ────────────────────────────────────────────────────────────────────
PROJECT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
AT_DIR  = os.path.join(PROJECT, "lineage_tracing/allele_tables/per_library")
QC_DIR  = os.path.join(PROJECT, "QC")
DROPBOX = "/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/Cassiopeia_QC"
STEM    = "editing_saturation_qc"

LIBRARIES = ["MJZ019", "MJZ021", "MJZ022", "MJZ023", "MJZ024", "MJZ025"]

ARM_ORDER  = ["Intact", "Castration", "Castration + Enzalutamide"]
ARM_SHORT  = {"Intact": "Intact", "Castration": "Cas", "Castration + Enzalutamide": "Cas+Enz"}
ARM_COLORS = {"Intact": "#4878CF", "Castration": "#D65F5F", "Castration + Enzalutamide": "#6ACC65"}

# ── LOAD METADATA ────────────────────────────────────────────────────────────
print("Loading metadata …")
meta = pd.read_csv(os.path.join(QC_DIR, "npp53_clean_49619_metadata.csv"))
meta["library"]  = meta["cell"].str.split("_").str[0]
meta["bare_bc"]  = meta["cell"].str.split("_").str[1].str.replace("-1", "")

# Only neoplasm cells (tumor cells that have the viral construct)
neoplasm_mask = ~meta["celltype"].isin(["Macrophage_Myeloid", "Endothelial"])
meta_neo = meta[neoplasm_mask].copy()

# Build bc_lookup: {library: {bare_bc: (sample_id, treatment)}}
bc_lookup = {}
for lib in LIBRARIES:
    m = meta_neo[meta_neo["library"] == lib]
    bc_lookup[lib] = {row["bare_bc"]: (row["sample_id"], row["treatment"])
                      for _, row in m.iterrows()}

print(f"Neoplasm barcodes: {sum(len(v) for v in bc_lookup.values())}")


# ── PARSE ALLELE TABLES ─────────────────────────────────────────────────────
print("\nParsing allele tables …")
t0 = time.time()

# Collect per-row data: sample, treatment, cellBC, site edited (T/F), allele string
rows_data = []  # list of dicts

for lib in LIBRARIES:
    bmap = bc_lookup.get(lib, {})
    if not bmap:
        continue
    atf = os.path.join(AT_DIR, f"{lib}_allele_table.csv")
    if not os.path.exists(atf):
        continue

    print(f"  Loading {lib} …", end="", flush=True)
    at = pd.read_csv(atf)
    bc_set = set(bmap.keys())
    at_match = at[at["cellBC"].isin(bc_set)].copy()
    if at_match.empty:
        print(" skip"); continue

    print(f" {len(at_match)} rows, {at_match['cellBC'].nunique()} cells")

    for _, row in at_match.iterrows():
        cbc = row["cellBC"]
        sid, trt = bmap[cbc]

        # For each site, check if edited (not "[None]" and not NaN/empty)
        r1_edited = "[None]" not in str(row["r1"]) and str(row["r1"]) not in ("nan", "")
        r2_edited = "[None]" not in str(row["r2"]) and str(row["r2"]) not in ("nan", "")
        r3_edited = "[None]" not in str(row["r3"]) and str(row["r3"]) not in ("nan", "")

        rows_data.append({
            "sample": sid,
            "treatment": trt,
            "cellBC": cbc,
            "intBC": row["intBC"],
            "r1_edited": r1_edited,
            "r2_edited": r2_edited,
            "r3_edited": r3_edited,
            "n_sites_edited": int(r1_edited) + int(r2_edited) + int(r3_edited),
            "allele": str(row.get("allele", "")),
        })

df = pd.DataFrame(rows_data)
elapsed = time.time() - t0
print(f"  Parsed {len(df)} rows in {elapsed:.1f}s")
print(f"  Unique cells: {df['cellBC'].nunique()}")
print(f"  Unique samples: {df['sample'].nunique()}")


# ═══════════════════════════════════════════════════════════════════════════════
# COMPUTE METRICS
# ═══════════════════════════════════════════════════════════════════════════════

# ── 1. Per-site editing rate (across all intBCs) by treatment arm ─────────
print("\n=== Per-site editing rate by arm ===")
site_rates = {}
for arm in ARM_ORDER:
    sub = df[df["treatment"] == arm]
    n = len(sub)
    r1 = sub["r1_edited"].sum() / n * 100
    r2 = sub["r2_edited"].sum() / n * 100
    r3 = sub["r3_edited"].sum() / n * 100
    site_rates[arm] = {"r1": r1, "r2": r2, "r3": r3, "n_rows": n}
    print(f"  {ARM_SHORT[arm]:8s}: r1={r1:.1f}% r2={r2:.1f}% r3={r3:.1f}% (n={n:,} intBC rows)")


# ── 2. Per-cell: fraction of intBCs with each site edited ────────────────
print("\n=== Per-cell site editing profile ===")
cell_site = df.groupby(["cellBC", "sample", "treatment"]).agg(
    n_intBC=("intBC", "count"),
    r1_edited_count=("r1_edited", "sum"),
    r2_edited_count=("r2_edited", "sum"),
    r3_edited_count=("r3_edited", "sum"),
    n_unique_alleles=("allele", "nunique"),
).reset_index()

cell_site["r1_frac"] = cell_site["r1_edited_count"] / cell_site["n_intBC"]
cell_site["r2_frac"] = cell_site["r2_edited_count"] / cell_site["n_intBC"]
cell_site["r3_frac"] = cell_site["r3_edited_count"] / cell_site["n_intBC"]
cell_site["mean_site_frac"] = (cell_site["r1_frac"] + cell_site["r2_frac"] + cell_site["r3_frac"]) / 3

# How many cells have ALL intBCs edited at a given site?
for site in ["r1", "r2", "r3"]:
    fully_edited = (cell_site[f"{site}_frac"] == 1.0).sum()
    pct = fully_edited / len(cell_site) * 100
    print(f"  {site}: {fully_edited:,} / {len(cell_site):,} cells ({pct:.1f}%) have ALL intBCs edited")

# Cells with at least 1 unedited intBC at any site
any_unedited = ((cell_site["r1_frac"] < 1) | (cell_site["r2_frac"] < 1) | (cell_site["r3_frac"] < 1)).sum()
pct_any = any_unedited / len(cell_site) * 100
print(f"  Cells with ≥1 unedited intBC at any site: {any_unedited:,} / {len(cell_site):,} ({pct_any:.1f}%)")

# Cells where ALL sites on ALL intBCs are edited (fully saturated)
fully_saturated = ((cell_site["r1_frac"] == 1) & (cell_site["r2_frac"] == 1) & (cell_site["r3_frac"] == 1)).sum()
pct_sat = fully_saturated / len(cell_site) * 100
print(f"  Fully saturated cells (all 3 sites × all intBCs edited): {fully_saturated:,} ({pct_sat:.1f}%)")


# ── 3. Per-sample editing rate for each site ─────────────────────────────
sample_site = df.groupby(["sample", "treatment"]).agg(
    n_rows=("intBC", "count"),
    r1_sum=("r1_edited", "sum"),
    r2_sum=("r2_edited", "sum"),
    r3_sum=("r3_edited", "sum"),
    n_unique_alleles=("allele", "nunique"),
    n_cells=("cellBC", "nunique"),
).reset_index()
sample_site["r1_rate"] = sample_site["r1_sum"] / sample_site["n_rows"] * 100
sample_site["r2_rate"] = sample_site["r2_sum"] / sample_site["n_rows"] * 100
sample_site["r3_rate"] = sample_site["r3_sum"] / sample_site["n_rows"] * 100
sample_site["mean_rate"] = (sample_site["r1_rate"] + sample_site["r2_rate"] + sample_site["r3_rate"]) / 3
sample_site["diversity"] = sample_site["n_unique_alleles"] / sample_site["n_cells"]

print("\n=== Per-sample editing rates (mean ± sd) ===")
for arm in ARM_ORDER:
    sub = sample_site[sample_site["treatment"] == arm]
    for site in ["r1", "r2", "r3"]:
        vals = sub[f"{site}_rate"]
        print(f"  {ARM_SHORT[arm]:8s} {site}: {vals.mean():.1f} ± {vals.std():.1f}%")


# ── 4. Per-row: number of sites edited (0, 1, 2, 3) ─────────────────────
print("\n=== Distribution of sites edited per intBC row ===")
for n in range(4):
    ct = (df["n_sites_edited"] == n).sum()
    pct = ct / len(df) * 100
    print(f"  {n} sites edited: {ct:,} ({pct:.1f}%)")


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT — 2×2 figure
# ═══════════════════════════════════════════════════════════════════════════════
print("\nPlotting …")

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
fig.suptitle("Cas9 Editing Saturation QC — Npp53 Lineage Tracing",
             fontsize=13, fontweight="bold", y=0.98)

# ── Panel A: Per-site editing rate by treatment arm (grouped bar) ────────
ax = axes[0, 0]
sites = ["r1", "r2", "r3"]
site_labels = ["Barcode-1", "Barcode-2", "Barcode-3"]
x = np.arange(len(sites))
bar_w = 0.22

for i, arm in enumerate(ARM_ORDER):
    rates = [site_rates[arm][s] for s in sites]
    offset = (i - 1) * bar_w
    bars = ax.bar(x + offset, rates, width=bar_w, color=ARM_COLORS[arm],
                  edgecolor="white", linewidth=0.5,
                  label=ARM_SHORT[arm], zorder=3)
    # Add value labels
    for b, v in zip(bars, rates):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.5,
                f"{v:.1f}%", ha="center", va="bottom", fontsize=6, fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(site_labels, fontsize=9)
ax.set_ylabel("Editing rate (%)", fontsize=10)
ax.set_ylim(0, 105)
ax.axhline(100, color="red", linestyle="--", linewidth=0.8, alpha=0.5, zorder=1)
ax.text(2.5, 101, "Saturation", color="red", fontsize=8, alpha=0.7, ha="right")
ax.legend(fontsize=8, loc="lower right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.text(-0.08, 1.05, "A", transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")
ax.set_title("Per-site editing rate (all intBC rows)", fontsize=10, pad=10)


# ── Panel B: Per-cell editing profile — fraction of intBCs edited per site
ax = axes[0, 1]

for arm in ARM_ORDER:
    sub = cell_site[cell_site["treatment"] == arm]
    for si, site in enumerate(["r1", "r2", "r3"]):
        vals = sub[f"{site}_frac"] * 100
        # Violin at position
        xpos = si * 4 + ARM_ORDER.index(arm)
        vp = ax.violinplot([vals.values], positions=[xpos], showmedians=False,
                           showextrema=False, widths=0.8)
        for pc in vp["bodies"]:
            pc.set_facecolor(ARM_COLORS[arm])
            pc.set_alpha(0.6)
        med = np.median(vals)
        ax.scatter(xpos, med, color="white", edgecolor="black", s=20, zorder=5, linewidth=0.8)

# x-axis labels
ax.set_xticks([1, 5, 9])
ax.set_xticklabels(site_labels, fontsize=9)
ax.set_ylabel("% of intBCs edited per cell", fontsize=10)
ax.set_ylim(-5, 110)
ax.axhline(100, color="red", linestyle="--", linewidth=0.8, alpha=0.5, zorder=1)

# Legend
handles = [mpatches.Patch(facecolor=ARM_COLORS[arm], alpha=0.6, label=ARM_SHORT[arm])
           for arm in ARM_ORDER]
ax.legend(handles=handles, fontsize=7, loc="lower left")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.text(-0.08, 1.05, "B", transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")
ax.set_title("Per-cell editing fraction (within intBCs)", fontsize=10, pad=10)


# ── Panel C: Per-intBC-row — sites edited (0, 1, 2, 3) stacked bar ──────
ax = axes[1, 0]

site_count_colors = {0: "#d9d9d9", 1: "#fdae6b", 2: "#6baed6", 3: "#2171b5"}
site_count_labels = {0: "0 sites", 1: "1 site", 2: "2 sites", 3: "3 sites (all)"}

x = np.arange(len(ARM_ORDER))
bottoms = np.zeros(len(ARM_ORDER))

for n_sites in [0, 1, 2, 3]:
    heights = []
    for arm in ARM_ORDER:
        sub = df[df["treatment"] == arm]
        pct = (sub["n_sites_edited"] == n_sites).sum() / len(sub) * 100
        heights.append(pct)
    ax.bar(x, heights, width=0.6, bottom=bottoms,
           color=site_count_colors[n_sites], edgecolor="white", linewidth=0.5,
           label=site_count_labels[n_sites], zorder=3)
    # Add % label inside bar if >3%
    for xi, h in enumerate(heights):
        if h > 3:
            ax.text(xi, bottoms[xi] + h/2, f"{h:.1f}%",
                    ha="center", va="center", fontsize=7, fontweight="bold",
                    color="white" if n_sites == 3 else "black")
    bottoms += heights

ax.set_xticks(x)
ax.set_xticklabels([ARM_SHORT[a] for a in ARM_ORDER], fontsize=9)
ax.set_ylabel("% of intBC rows", fontsize=10)
ax.set_ylim(0, 105)
ax.legend(fontsize=7, loc="upper right", title="Sites edited per intBC", title_fontsize=8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.text(-0.08, 1.05, "C", transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")
ax.set_title("Sites edited per intBC row", fontsize=10, pad=10)


# ── Panel D: Allele diversity vs mean editing rate (per sample) ──────────
ax = axes[1, 1]

for arm in ARM_ORDER:
    sub = sample_site[sample_site["treatment"] == arm]
    ax.scatter(sub["mean_rate"], sub["diversity"],
               c=ARM_COLORS[arm], s=50, alpha=0.8, edgecolors="black",
               linewidth=0.5, label=ARM_SHORT[arm], zorder=3)

# Correlation
from scipy import stats as sp_stats
r, p = sp_stats.pearsonr(sample_site["mean_rate"], sample_site["diversity"])
ax.text(0.05, 0.95, f"r = {r:.2f}, p = {p:.2e}",
        transform=ax.transAxes, fontsize=8, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

ax.set_xlabel("Mean editing rate across 3 sites (%)", fontsize=10)
ax.set_ylabel("Allele diversity\n(unique alleles / cells per sample)", fontsize=10)
ax.legend(fontsize=8, loc="lower right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.text(-0.08, 1.05, "D", transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")
ax.set_title("Allele diversity vs editing rate", fontsize=10, pad=10)


# ── SAVE ─────────────────────────────────────────────────────────────────────
plt.tight_layout(rect=[0, 0, 1, 0.95])

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

# ── SUMMARY ──────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("SATURATION ASSESSMENT SUMMARY")
print("="*70)
overall_r1 = df["r1_edited"].mean() * 100
overall_r2 = df["r2_edited"].mean() * 100
overall_r3 = df["r3_edited"].mean() * 100
print(f"Overall per-site editing rate: r1={overall_r1:.1f}%, r2={overall_r2:.1f}%, r3={overall_r3:.1f}%")
pct_3sites = (df["n_sites_edited"] == 3).mean() * 100
pct_0sites = (df["n_sites_edited"] == 0).mean() * 100
print(f"intBC rows with all 3 sites edited: {pct_3sites:.1f}%")
print(f"intBC rows with 0 sites edited (uncut): {pct_0sites:.1f}%")
print(f"Fully saturated cells (all sites × all intBCs): {fully_saturated:,} / {len(cell_site):,} ({pct_sat:.1f}%)")
if pct_3sites > 90:
    print("→ WARNING: Editing appears NEAR SATURATION — most intBCs fully edited")
elif pct_3sites > 70:
    print("→ CAUTION: Editing is HIGH but not fully saturated — some diversity headroom remains")
else:
    print("→ GOOD: Editing is NOT saturated — substantial diversity headroom for lineage resolution")
print("="*70)
print("Done.")

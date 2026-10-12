#!/usr/bin/env python3
# =============================================================================
# Per-sample editing rate bar plot — split by r1/r2/r3
# =============================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os, shutil, time

PROJECT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
AT_DIR  = os.path.join(PROJECT, "lineage_tracing/allele_tables/per_library")
QC_DIR  = os.path.join(PROJECT, "QC")
DROPBOX = "/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/Cassiopeia_QC"

LIBRARIES = ["MJZ019", "MJZ021", "MJZ022", "MJZ023", "MJZ024", "MJZ025"]
ARM_ORDER  = ["Intact", "Castration", "Castration + Enzalutamide"]
ARM_COLORS = {"Intact": "#4878CF", "Castration": "#D65F5F", "Castration + Enzalutamide": "#6ACC65"}

# ── LOAD METADATA ────────────────────────────────────────────────────────────
print("Loading metadata …")
meta = pd.read_csv(os.path.join(QC_DIR, "npp53_clean_49619_metadata.csv"))
meta["library"]  = meta["cell"].str.split("_").str[0]
meta["bare_bc"]  = meta["cell"].str.split("_").str[1].str.replace("-1", "")

neoplasm_mask = ~meta["celltype"].isin(["Macrophage_Myeloid", "Endothelial"])
meta_neo = meta[neoplasm_mask].copy()

bc_lookup = {}
for lib in LIBRARIES:
    m = meta_neo[meta_neo["library"] == lib]
    bc_lookup[lib] = {row["bare_bc"]: (row["sample_id"], row["treatment"])
                      for _, row in m.iterrows()}

# ── PARSE ALLELE TABLES ─────────────────────────────────────────────────────
print("Parsing allele tables …")
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
    at_match = at[at["cellBC"].isin(set(bmap.keys()))].copy()
    if at_match.empty:
        print(" skip"); continue
    print(f" {len(at_match)} rows")

    for _, row in at_match.iterrows():
        cbc = row["cellBC"]
        sid, trt = bmap[cbc]
        r1_ed = "[None]" not in str(row["r1"]) and str(row["r1"]) not in ("nan", "")
        r2_ed = "[None]" not in str(row["r2"]) and str(row["r2"]) not in ("nan", "")
        r3_ed = "[None]" not in str(row["r3"]) and str(row["r3"]) not in ("nan", "")
        rows_data.append({"sample": sid, "treatment": trt,
                          "r1_edited": r1_ed, "r2_edited": r2_ed, "r3_edited": r3_ed})

df = pd.DataFrame(rows_data)
print(f"  {len(df)} total rows, {elapsed:.1f}s" if False else f"  {len(df)} total rows in {time.time()-t0:.1f}s")

# ── COMPUTE PER-SAMPLE EDITING RATE ─────────────────────────────────────────
sample_stats = df.groupby(["sample", "treatment"]).agg(
    n_rows=("r1_edited", "count"),
    r1_sum=("r1_edited", "sum"),
    r2_sum=("r2_edited", "sum"),
    r3_sum=("r3_edited", "sum"),
).reset_index()
sample_stats["r1_rate"] = sample_stats["r1_sum"] / sample_stats["n_rows"] * 100
sample_stats["r2_rate"] = sample_stats["r2_sum"] / sample_stats["n_rows"] * 100
sample_stats["r3_rate"] = sample_stats["r3_sum"] / sample_stats["n_rows"] * 100

# Sort by arm order then sample name
sample_stats["arm_idx"] = sample_stats["treatment"].map({a: i for i, a in enumerate(ARM_ORDER)})
sample_stats = sample_stats.sort_values(["arm_idx", "sample"]).reset_index(drop=True)

# Print the calculation detail for each sample
print("\n" + "="*90)
print(f"{'Sample':<15} {'Arm':<12} {'intBC rows':>10}  {'r1 edited':>10} {'r1 rate':>8}  "
      f"{'r2 edited':>10} {'r2 rate':>8}  {'r3 edited':>10} {'r3 rate':>8}")
print("="*90)
for _, r in sample_stats.iterrows():
    print(f"{r['sample']:<15} {r['treatment'][:11]:<12} {int(r['n_rows']):>10,}  "
          f"{int(r['r1_sum']):>10,} {r['r1_rate']:>7.1f}%  "
          f"{int(r['r2_sum']):>10,} {r['r2_rate']:>7.1f}%  "
          f"{int(r['r3_sum']):>10,} {r['r3_rate']:>7.1f}%")
print("="*90)

# ── PLOT: 3-panel (r1 / r2 / r3), x=sample, bar colored by treatment ────────
fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
fig.suptitle("Per-sample Cas9 editing rate at each barcode site",
             fontsize=13, fontweight="bold", y=0.98)

sites = [("r1_rate", "Barcode-1 (r1)"), ("r2_rate", "Barcode-2 (r2)"), ("r3_rate", "Barcode-3 (r3)")]
n_samples = len(sample_stats)
x = np.arange(n_samples)

for ax_i, (col, title) in enumerate(sites):
    ax = axes[ax_i]
    colors = [ARM_COLORS[t] for t in sample_stats["treatment"]]
    bars = ax.bar(x, sample_stats[col], width=0.7, color=colors,
                  edgecolor="white", linewidth=0.5, zorder=3)

    # Add value labels on top of bars
    for xi, val in enumerate(sample_stats[col]):
        ax.text(xi, val + 0.8, f"{val:.0f}", ha="center", va="bottom",
                fontsize=6, fontweight="bold", rotation=0)

    ax.set_ylabel("Editing rate (%)", fontsize=10)
    ax.set_ylim(0, 105)
    ax.axhline(100, color="red", linestyle="--", linewidth=0.6, alpha=0.4)

    # Arm mean lines
    for arm in ARM_ORDER:
        mask = sample_stats["treatment"] == arm
        if mask.sum() > 0:
            mean_val = sample_stats.loc[mask, col].mean()
            idx = np.where(mask)[0]
            ax.hlines(mean_val, idx[0] - 0.45, idx[-1] + 0.45,
                      color=ARM_COLORS[arm], linestyle="-", linewidth=1.5,
                      alpha=0.7, zorder=4)

    ax.set_title(title, fontsize=11, fontweight="bold", loc="left", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Add arm separators
    arm_changes = []
    for i in range(1, n_samples):
        if sample_stats.iloc[i]["treatment"] != sample_stats.iloc[i-1]["treatment"]:
            arm_changes.append(i - 0.5)
    for xc in arm_changes:
        ax.axvline(xc, color="grey", linestyle=":", linewidth=0.8, alpha=0.5)

# X labels on bottom panel only
axes[-1].set_xticks(x)
axes[-1].set_xticklabels(sample_stats["sample"], rotation=45, ha="right", fontsize=8)

# Legend
handles = [mpatches.Patch(facecolor=ARM_COLORS[a], label=a) for a in ARM_ORDER]
axes[0].legend(handles=handles, fontsize=8, loc="lower right", ncol=3,
               bbox_to_anchor=(1.0, 1.02), framealpha=0.9)

plt.tight_layout(rect=[0, 0, 1, 0.95])

stem = "editing_rate_per_sample_by_site"
out_png = os.path.join(QC_DIR, f"{stem}.png")
out_pdf = os.path.join(QC_DIR, f"{stem}.pdf")
fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
print(f"\nSaved: {out_png}")
print(f"Saved: {out_pdf}")

for ext in [".png", ".pdf"]:
    shutil.copy2(os.path.join(QC_DIR, f"{stem}{ext}"),
                 os.path.join(DROPBOX, f"{stem}{ext}"))
    print(f"Synced to Dropbox")

plt.close()
print("Done.")

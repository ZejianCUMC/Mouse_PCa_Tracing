#!/usr/bin/env python3
# =============================================================================
# Manhattan-style Mutation Frequency Plot — Lu et al. 2024 Fig 1d analog
# =============================================================================
# Description : Two-panel Manhattan plot of per-site editing frequency
#               along barcode coordinates (~100–230 bp).
#               Upper panel: Normal (Macrophage_Myeloid + Endothelial)
#               Lower panel: Neoplasm (GFP+/mCherry+ tumor cells)
#               Each dot = one mouse sample at one edit position.
# Input       : per_library allele tables, npp53_clean_49619_metadata.csv
# Output      : manhattan_mutation_freq_like_LuNat_Fig1d.png/.pdf
# Conda env   : scvi
# Style ref   : Lu et al. 2024 Nature 636:233–240, Fig 1d
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re, os, shutil
from collections import defaultdict

# ── PATHS ────────────────────────────────────────────────────────────────────
PROJECT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
AT_DIR  = os.path.join(PROJECT, "lineage_tracing/allele_tables/per_library")
QC_DIR  = os.path.join(PROJECT, "QC")
DROPBOX = "/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/Cassiopeia_QC"
STEM    = "manhattan_mutation_freq_like_LuNat_Fig1d"

# ── BARCODE REGION DEFINITIONS ───────────────────────────────────────────────
# Real barcode coordinates from allele encoding analysis
REGIONS = {
    "r1": {"start": 89, "end": 125, "center": 113, "color": "#4878CF", "label": "Barcode-1"},
    "r2": {"start": 140, "end": 179, "center": 159, "color": "#D65F5F", "label": "Barcode-2"},
    "r3": {"start": 200, "end": 233, "center": 221, "color": "#6ACC65", "label": "Barcode-3"},
}

LIBRARIES = ["MJZ019", "MJZ021", "MJZ022", "MJZ023", "MJZ024", "MJZ025"]
EDIT_PAT  = re.compile(r'\[(\d+):(\d+)([DIS])\]')

# ── LOAD METADATA ────────────────────────────────────────────────────────────
print("Loading metadata …")
meta = pd.read_csv(os.path.join(QC_DIR, "npp53_clean_49619_metadata.csv"))

# Parse library + bare barcode from cell name
meta["library"]  = meta["cell"].str.split("_").str[0]
meta["bare_bc"]  = meta["cell"].str.split("_").str[1].str.replace("-1", "")

# Define groups
normal_mask = meta["celltype"].isin(["Macrophage_Myeloid", "Endothelial"])
neoplasm_mask = ~normal_mask

normal_bcs  = {}  # library -> {bare_bc: sample_id}
neoplasm_bcs = {}
for lib in LIBRARIES:
    m_n = meta[(meta["library"] == lib) & normal_mask]
    normal_bcs[lib] = dict(zip(m_n["bare_bc"], m_n["sample_id"]))
    m_t = meta[(meta["library"] == lib) & neoplasm_mask]
    neoplasm_bcs[lib] = dict(zip(m_t["bare_bc"], m_t["sample_id"]))

print(f"Normal barcodes:   {sum(len(v) for v in normal_bcs.values())} across libraries")
print(f"Neoplasm barcodes: {sum(len(v) for v in neoplasm_bcs.values())} across libraries")

# ── PARSE ALLELE TABLES — PER-SAMPLE, PER-POSITION MUTATION FREQUENCY ───────
def _affected_positions(edit_str, site):
    """Return frozenset of barcode-region positions affected by an edit string.

    Deletions are expanded to all covered positions within the region.
    """
    val = str(edit_str)
    if "[None]" in val:
        return frozenset()
    matches = EDIT_PAT.findall(val)
    if not matches:
        return frozenset()
    reg_s, reg_e = REGIONS[site]["start"], REGIONS[site]["end"]
    positions = set()
    for pos_str, length_str, etype in matches:
        pos = int(pos_str)
        length = int(length_str)
        if etype in ("D", "S"):
            # Deletion/substitution spans pos to pos+length-1
            for p in range(max(pos, reg_s), min(pos + length, reg_e + 1)):
                positions.add(p)
        else:  # Insertion — affects one position
            if reg_s <= pos <= reg_e:
                positions.add(pos)
    return frozenset(positions)


def parse_per_sample_position(bc_lookup, label=""):
    """Per-sample, per-position mutation frequency.

    Each dot = one (sample, position) combination.
    Frequency = fraction of that sample's cell×intBC rows where an edit
    covers that position (deletions expanded to all spanned positions).
    """
    # Pre-accumulate per (sample, site, position) → affected_count
    # and per (sample, site) → total_rows
    sp_affected = defaultdict(lambda: defaultdict(int))  # (sid,site) -> {pos: count}
    sp_total    = defaultdict(lambda: defaultdict(int))  # sid -> {site: total}
    sample_cells = defaultdict(set)

    for lib in LIBRARIES:
        bc_map = bc_lookup.get(lib, {})
        if not bc_map:
            continue
        atf = os.path.join(AT_DIR, f"{lib}_allele_table.csv")
        if not os.path.exists(atf):
            continue

        bc_set = set(bc_map.keys())
        at = pd.read_csv(atf)
        at_match = at[at["cellBC"].isin(bc_set)].copy()
        if at_match.empty:
            continue

        # Map cellBC → sample_id
        at_match["sid"] = at_match["cellBC"].map(bc_map)

        for site in ["r1", "r2", "r3"]:
            # Parse unique edit strings once (much faster than per-row)
            unique_edits = at_match[site].unique()
            edit_to_pos = {e: _affected_positions(e, site) for e in unique_edits}
            at_match[f"_pos_{site}"] = at_match[site].map(edit_to_pos)

        print(f"    {lib}: {len(at_match)} matched rows ({label})")

        # Accumulate counts
        for _, row in at_match.iterrows():
            sid = row["sid"]
            sample_cells[sid].add(row["cellBC"])
            for site in ["r1", "r2", "r3"]:
                sp_total[sid][site] += 1
                for p in row[f"_pos_{site}"]:
                    sp_affected[(sid, site)][p] += 1

    # Build dots: one per (sample, position)
    dots = []
    for (sid, site), pos_counts in sp_affected.items():
        total = sp_total[sid][site]
        n_cells = len(sample_cells[sid])
        for pos, count in pos_counts.items():
            freq = count / total if total > 0 else 0
            dots.append((pos, freq, sid, site, n_cells))

    stats = {sid: len(cells) for sid, cells in sample_cells.items()}
    return dots, stats


print("Parsing Normal per-sample per-position …")
normal_dots_raw, normal_stats = parse_per_sample_position(normal_bcs, "Normal")
print(f"  {len(normal_stats)} samples, {sum(normal_stats.values())} cells, "
      f"{len(normal_dots_raw)} dots")

print("Parsing Neoplasm per-sample per-position …")
neoplasm_dots_raw, neoplasm_stats = parse_per_sample_position(neoplasm_bcs, "Neoplasm")
print(f"  {len(neoplasm_stats)} samples, {sum(neoplasm_stats.values())} cells, "
      f"{len(neoplasm_dots_raw)} dots")

# ── PLOT ─────────────────────────────────────────────────────────────────────
fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(12, 7), sharex=True,
    gridspec_kw={"height_ratios": [1, 1], "hspace": 0.08}
)

fig.suptitle("Npp53 Lineage Tracing — Mutation Frequency Along Barcode Coordinates\n"
             "(cf. Lu et al. 2024 Nature Fig 1d)",
             fontsize=11, fontweight="bold", y=0.98)

# Sample → treatment for color reference
sample_treatment = dict(zip(meta["sample_id"], meta["treatment"]))

def plot_panel(ax, dots_raw, stats, title, invert=False):
    """Plot one Manhattan panel. Each dot = one (sample, site).

    dots_raw: list of (sample_id, site_name, editing_rate, n_cells)
    x = cut site center (jittered within region for visual spread)
    y = per-site editing rate for that sample
    """
    if not dots_raw:
        ax.set_ylabel("Mutation frequency")
        ax.text(0.5, 0.5, f"No data ({title})", transform=ax.transAxes,
                ha="center", va="center", fontsize=12, color="grey")
        return

    rng = np.random.default_rng(42)
    xs, ys, colors, sizes = [], [], [], []

    for sid, site, rate, n_cells in dots_raw:
        rinfo = REGIONS[site]
        # Jitter x within region for visual spread
        x = rinfo["center"] + rng.uniform(-8, 8)
        xs.append(x)
        ys.append(rate)
        colors.append(rinfo["color"])
        # Scale dot size by cell count (larger = more cells = more reliable)
        sizes.append(max(15, min(80, n_cells * 0.8)))

    xs = np.array(xs)
    ys = np.array(ys)

    if invert:
        ys_plot = -ys
    else:
        ys_plot = ys

    ax.scatter(xs, ys_plot, c=colors, s=sizes, alpha=0.65,
               edgecolors="black", linewidths=0.4, zorder=3)

    ax.set_ylabel("Mutation frequency", fontsize=10)

    n_samp = len(stats)
    n_cells_total = sum(stats.values())

    if invert:
        ax.set_ylim(-1.05, 0.05)
        ax.set_yticks([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0])
        ax.set_yticklabels(["1.0", "0.8", "0.6", "0.4", "0.2", "0.0"])
        ax.text(0.02, 0.05, f"{title}\n({n_cells_total:,} cells, {n_samp} samples)",
                transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    else:
        ax.set_ylim(-0.05, 1.05)
        ax.text(0.02, 0.95, f"{title}\n({n_cells_total:,} cells, {n_samp} samples)",
                transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# Upper: Normal
plot_panel(ax_top, normal_dots_raw, normal_stats,
           "Normal (GFP⁻/mCherry⁻)", invert=False)

# Lower: Neoplasm (inverted y like Lu)
plot_panel(ax_bot, neoplasm_dots_raw, neoplasm_stats,
           "Neoplasm (GFP⁺/mCherry⁺)", invert=True)

# ── BARCODE REGION BANDS (at the x-axis junction) ───────────────────────────
x_min = 82   # just before r1 start (89)
x_max = 240  # just after r3 end (233)
ax_top.set_xlim(x_min, x_max)
ax_bot.set_xlim(x_min, x_max)

# Draw colored bands at top panel bottom and bot panel top
band_height_top = 0.03  # fraction of y range
band_height_bot = 0.03

for site_name, rinfo in REGIONS.items():
    w = rinfo["end"] - rinfo["start"]
    # Band on top panel (below y=0)
    rect_top = mpatches.Rectangle(
        (rinfo["start"], -0.06), w, 0.06,
        facecolor=rinfo["color"], alpha=0.8,
        edgecolor="black", linewidth=0.5, zorder=5
    )
    ax_top.add_patch(rect_top)
    ax_top.text((rinfo["start"] + rinfo["end"]) / 2, -0.03,
                rinfo["label"], ha="center", va="center", fontsize=6,
                fontweight="bold", color="white", zorder=6)

    # Band on bottom panel (above y=0)
    rect_bot = mpatches.Rectangle(
        (rinfo["start"], 0.0), w, 0.06,
        facecolor=rinfo["color"], alpha=0.8,
        edgecolor="black", linewidth=0.5, zorder=5
    )
    ax_bot.add_patch(rect_bot)
    ax_bot.text((rinfo["start"] + rinfo["end"]) / 2, 0.03,
                rinfo["label"], ha="center", va="center", fontsize=6,
                fontweight="bold", color="white", zorder=6)

ax_bot.set_xlabel("Position along lineage barcode construct (bp)", fontsize=10)

# ── SAVE ─────────────────────────────────────────────────────────────────────
os.makedirs(QC_DIR, exist_ok=True)
out_png = os.path.join(QC_DIR, f"{STEM}.png")
out_pdf = os.path.join(QC_DIR, f"{STEM}.pdf")

fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")

# Copy to Dropbox
for ext in [".png", ".pdf"]:
    src = os.path.join(QC_DIR, f"{STEM}{ext}")
    dst = os.path.join(DROPBOX, f"{STEM}{ext}")
    shutil.copy2(src, dst)
    print(f"Synced: {dst}")

plt.close()
print("Done.")

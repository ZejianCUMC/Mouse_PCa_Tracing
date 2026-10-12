#!/usr/bin/env python3
# =============================================================================
# Manhattan Stacked-Bar Mutation Spectrum — Lu et al. 2024 Fig 1d analog
# =============================================================================
# Description : Per-position stacked bar plot of mean edited-cell counts
#               along barcode coordinates. Bars stacked by edit type
#               (Deletion / Insertion / Substitution / Complex).
#               Upper panel: Normal (broken y-axis), Lower: Neoplasm (mirrored).
# Input       : per_library allele tables, npp53_clean_49619_metadata.csv
# Output      : manhattan_stacked_mutation_spectrum.png/.pdf
# Conda env   : scvi
# Style ref   : Lu et al. 2024 Nature 636:233–240, Fig 1d
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import re, os, shutil
from collections import defaultdict
import time

# ── PATHS ────────────────────────────────────────────────────────────────────
PROJECT = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
AT_DIR  = os.path.join(PROJECT, "lineage_tracing/allele_tables/per_library")
QC_DIR  = os.path.join(PROJECT, "QC")
DROPBOX = "/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/Cassiopeia_QC"
STEM    = "manhattan_stacked_mutation_spectrum"

# ── BARCODE REGION DEFINITIONS ───────────────────────────────────────────────
REGIONS = {
    "r1": {"start": 89, "end": 125, "center": 113, "color": "#4878CF", "label": "Barcode-1"},
    "r2": {"start": 140, "end": 179, "center": 159, "color": "#D65F5F", "label": "Barcode-2"},
    "r3": {"start": 200, "end": 233, "center": 221, "color": "#6ACC65", "label": "Barcode-3"},
}

# All integer positions within barcode regions
ALL_POSITIONS = []
POS_TO_SITE = {}
for site, rinfo in REGIONS.items():
    for p in range(rinfo["start"], rinfo["end"] + 1):
        ALL_POSITIONS.append(p)
        POS_TO_SITE[p] = site
ALL_POSITIONS = sorted(ALL_POSITIONS)

LIBRARIES = ["MJZ019", "MJZ021", "MJZ022", "MJZ023", "MJZ024", "MJZ025"]
EDIT_PAT  = re.compile(r'\[(\d+):(\d+)([DIS])\]')

# Edit type colors
ETYPE_COLORS = {
    "D": "#3274A1",   # steel blue — deletion
    "I": "#E1812C",   # orange — insertion
    "S": "#3A923A",   # green — substitution
    "X": "#9467BD",   # purple — complex (multiple types)
}
ETYPE_LABELS = {"D": "Deletion", "I": "Insertion", "S": "Substitution", "X": "Complex"}

# ── LOAD METADATA ────────────────────────────────────────────────────────────
print("Loading metadata …")
meta = pd.read_csv(os.path.join(QC_DIR, "npp53_clean_49619_metadata.csv"))
meta["library"]  = meta["cell"].str.split("_").str[0]
meta["bare_bc"]  = meta["cell"].str.split("_").str[1].str.replace("-1", "")

normal_mask = meta["celltype"].isin(["Macrophage_Myeloid", "Endothelial"])

# Build bc_lookup: {library: {bare_bc: sample_id}}
normal_bcs, neoplasm_bcs = {}, {}
for lib in LIBRARIES:
    m_n = meta[(meta["library"] == lib) & normal_mask]
    normal_bcs[lib] = dict(zip(m_n["bare_bc"], m_n["sample_id"]))
    m_t = meta[(meta["library"] == lib) & ~normal_mask]
    neoplasm_bcs[lib] = dict(zip(m_t["bare_bc"], m_t["sample_id"]))

print(f"Normal barcodes:   {sum(len(v) for v in normal_bcs.values())}")
print(f"Neoplasm barcodes: {sum(len(v) for v in neoplasm_bcs.values())}")


# ── HELPER: parse edit string → {position: edit_type} ───────────────────────
def _parse_edit_positions(edit_str, site):
    """Return dict {position: edit_type} for positions within the barcode region.

    Deletions are expanded to all covered positions.
    If multiple edits in one string cover the same position, mark as 'X' (complex).
    """
    val = str(edit_str)
    if "[None]" in val:
        return {}
    matches = EDIT_PAT.findall(val)
    if not matches:
        return {}

    reg_s, reg_e = REGIONS[site]["start"], REGIONS[site]["end"]
    pos_types = {}  # position -> set of types

    for pos_str, length_str, etype in matches:
        pos = int(pos_str)
        length = int(length_str)
        if etype == "D":
            for p in range(max(pos, reg_s), min(pos + length, reg_e + 1)):
                pos_types.setdefault(p, set()).add("D")
        elif etype == "S":
            for p in range(max(pos, reg_s), min(pos + length, reg_e + 1)):
                pos_types.setdefault(p, set()).add("S")
        else:  # I
            if reg_s <= pos <= reg_e:
                pos_types.setdefault(pos, set()).add("I")

    # Collapse to single type or complex
    result = {}
    for p, types in pos_types.items():
        if len(types) == 1:
            result[p] = types.pop()
        else:
            result[p] = "X"
    return result


# ── MAIN COMPUTATION ────────────────────────────────────────────────────────
def compute_per_sample_position_counts(bc_lookup, label=""):
    """For each sample, at each barcode position, count UNIQUE CELLS
    with an edit covering that position, classified by edit type.

    Returns:
        sample_pos_type_count: {sample_id: {position: {etype: count}}}
        sample_total_cells: {sample_id: total_unique_cells}
    """
    # Per (sample, cell) → {position: set_of_types_across_intBCs}
    cell_pos_types = defaultdict(lambda: defaultdict(set))
    sample_cells = defaultdict(set)

    t0 = time.time()

    for lib in LIBRARIES:
        bc_map = bc_lookup.get(lib, {})
        if not bc_map:
            continue
        atf = os.path.join(AT_DIR, f"{lib}_allele_table.csv")
        if not os.path.exists(atf):
            continue

        print(f"  Loading {lib} …", end="", flush=True)
        at = pd.read_csv(atf)
        bc_set = set(bc_map.keys())
        at_match = at[at["cellBC"].isin(bc_set)].copy()
        if at_match.empty:
            print(" 0 rows, skip")
            continue

        at_match["sid"] = at_match["cellBC"].map(bc_map)
        print(f" {len(at_match)} rows, {at_match['cellBC'].nunique()} cells")

        # Pre-parse unique edit strings for each site (big speedup)
        site_edit_maps = {}
        for site in ["r1", "r2", "r3"]:
            unique_edits = at_match[site].unique()
            site_edit_maps[site] = {e: _parse_edit_positions(e, site) for e in unique_edits}

        # Process each row: accumulate per-cell position types
        for _, row in at_match.iterrows():
            sid = row["sid"]
            cbc = row["cellBC"]
            sample_cells[sid].add(cbc)
            cell_key = (sid, cbc)

            for site in ["r1", "r2", "r3"]:
                pos_types = site_edit_maps[site].get(row[site], {})
                for p, etype in pos_types.items():
                    cell_pos_types[cell_key][p].add(etype)

    elapsed = time.time() - t0
    print(f"  Row processing done in {elapsed:.1f}s")

    # Collapse per-cell to final type and count per sample
    print("  Collapsing per-cell types …")
    sample_pos_type_count = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    sample_total_cells = {sid: len(cells) for sid, cells in sample_cells.items()}

    for (sid, cbc), pos_type_sets in cell_pos_types.items():
        for p, types in pos_type_sets.items():
            if len(types) == 1:
                final_type = types.pop()
            else:
                final_type = "X"
            sample_pos_type_count[sid][p][final_type] += 1

    n_samples = len(sample_total_cells)
    n_cells = sum(sample_total_cells.values())
    print(f"  {label}: {n_samples} samples, {n_cells} cells")

    return sample_pos_type_count, sample_total_cells


print("\n=== Normal ===")
normal_counts, normal_cells = compute_per_sample_position_counts(normal_bcs, "Normal")

print("\n=== Neoplasm ===")
neo_counts, neo_cells = compute_per_sample_position_counts(neoplasm_bcs, "Neoplasm")


# ── COMPUTE MEAN ± SEM ACROSS SAMPLES ───────────────────────────────────────
def compute_mean_sem(sample_pos_type_count, sample_total_cells):
    """For each position and edit type, compute mean and SEM of cell counts
    across samples.

    Returns:
        mean_dict: {position: {etype: mean_count}}
        sem_dict:  {position: {etype: sem_count}}
    """
    samples = sorted(sample_total_cells.keys())
    n = len(samples)
    etypes = ["D", "I", "S", "X"]

    mean_dict = defaultdict(lambda: defaultdict(float))
    sem_dict  = defaultdict(lambda: defaultdict(float))

    for pos in ALL_POSITIONS:
        for et in etypes:
            vals = [sample_pos_type_count[sid][pos].get(et, 0) for sid in samples]
            vals = np.array(vals, dtype=float)
            mean_dict[pos][et] = vals.mean()
            sem_dict[pos][et] = vals.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0

    return mean_dict, sem_dict


print("\nComputing means …")
normal_mean, normal_sem = compute_mean_sem(normal_counts, normal_cells)
neo_mean, neo_sem = compute_mean_sem(neo_counts, neo_cells)

# Print summary stats
for label, mean_d in [("Normal", normal_mean), ("Neoplasm", neo_mean)]:
    peak_pos = max(ALL_POSITIONS, key=lambda p: sum(mean_d[p].values()))
    peak_val = sum(mean_d[peak_pos].values())
    print(f"  {label} peak: pos {peak_pos}, mean {peak_val:.1f} cells")


# ── PLOT ─────────────────────────────────────────────────────────────────────
print("\nPlotting …")

# Determine y-axis limits from data
neo_peak = max(sum(neo_mean[p].values()) for p in ALL_POSITIONS)
normal_peak = max(sum(normal_mean[p].values()) for p in ALL_POSITIONS)
print(f"  Neoplasm peak mean: {neo_peak:.1f}")
print(f"  Normal peak mean: {normal_peak:.1f}")

# Y-axis scale: Neoplasm ymax = 2000
ymax = 2000
STEP = 500  # tick step

# Normal ylim: just above the Normal data, rounded up
normal_ylim = int(np.ceil(normal_peak / 5) * 5) + 5  # e.g. 20

print(f"  ymax: {ymax}")
print(f"  Normal ylim: {normal_ylim}")

# ── BUILD CONTIGUOUS X-AXIS (remove inter-barcode gaps) ─────────────────────
# Map real positions to contiguous display positions with small separator
SEPARATOR = 3  # empty gap between regions for visual separation
display_positions = []   # contiguous x values for plotting
real_positions = []      # actual bp positions (for tick labels)
pos_to_display = {}      # real_pos -> display_x

offset = 0
region_display_ranges = {}  # site -> (display_start, display_end)
for site in ["r1", "r2", "r3"]:
    rinfo = REGIONS[site]
    start_disp = offset
    for p in range(rinfo["start"], rinfo["end"] + 1):
        pos_to_display[p] = offset
        display_positions.append(offset)
        real_positions.append(p)
        offset += 1
    region_display_ranges[site] = (start_disp, offset - 1)
    offset += SEPARATOR  # gap between regions

DISPLAY_POS = display_positions  # for bar plotting

# ── LAYOUT ───────────────────────────────────────────────────────────────────
# Two panels: Normal (top, bars up) and Neoplasm (bottom, bars down/mirrored)
# No compressed panel — removed per user request

fig = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 1, height_ratios=[1, 7], hspace=0.08)
ax_norm = fig.add_subplot(gs[0])     # Normal
ax_neo  = fig.add_subplot(gs[1], sharex=ax_norm)  # Neoplasm

fig.suptitle("Npp53 Lineage Tracing — Per-Position Mutation Spectrum\n"
             "(cf. Lu et al. 2024 Nature Fig 1d)",
             fontsize=12, fontweight="bold", y=0.98)


BAR_WIDTH = 0.65  # < 1.0 for small gaps between bars

def plot_stacked_bars(ax, mean_dict, real_pos_list, disp_pos_list, bottom_up=True):
    """Draw stacked bars using display (contiguous) x positions."""
    etypes = ["D", "I", "S", "X"]
    x = np.array(disp_pos_list, dtype=float)
    cumulative = np.zeros(len(disp_pos_list))

    for et in etypes:
        heights = np.array([mean_dict[p].get(et, 0) for p in real_pos_list])
        if bottom_up:
            ax.bar(x, heights, width=BAR_WIDTH, bottom=cumulative,
                   color=ETYPE_COLORS[et], edgecolor="none",
                   label=ETYPE_LABELS[et], linewidth=0, zorder=3)
        else:
            ax.bar(x, -heights, width=BAR_WIDTH, bottom=-cumulative,
                   color=ETYPE_COLORS[et], edgecolor="none",
                   label=ETYPE_LABELS[et], linewidth=0, zorder=3)
        cumulative += heights


# ── Normal panel (bars going UP) ────────────────────────────────────────────
plot_stacked_bars(ax_norm, normal_mean, real_positions, DISPLAY_POS, bottom_up=True)

ax_norm.set_ylim(0, normal_ylim)
step_nz = max(1, normal_ylim // 4)
ax_norm.set_yticks(range(0, normal_ylim + 1, step_nz))
ax_norm.spines["top"].set_visible(False)
ax_norm.spines["right"].set_visible(False)
ax_norm.spines["bottom"].set_visible(False)  # hide bottom spine (shared boundary)
ax_norm.tick_params(bottom=False, labelbottom=False)


# ── Neoplasm panel (mirrored, bars going DOWN) ──────────────────────────────
plot_stacked_bars(ax_neo, neo_mean, real_positions, DISPLAY_POS, bottom_up=False)

ax_neo.set_ylim(-ymax, 0)
neo_ticks = list(range(0, ymax + 1, STEP))
ax_neo.set_yticks([-t for t in neo_ticks])
ax_neo.set_yticklabels([str(t) for t in neo_ticks])

# X-axis line at TOP of Neoplasm panel (y=0 boundary), not bottom
ax_neo.spines["top"].set_visible(True)
ax_neo.spines["bottom"].set_visible(False)
ax_neo.spines["right"].set_visible(False)
ax_neo.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)


# ── Y-axis labels: "Normal" and "Neoplasm" subtitles + shared ylabel ────────
fig.text(0.015, 0.50, "Mean edited cells per sample", fontsize=10,
         ha="center", va="center", rotation=90)
# "Normal" subtitle near the Normal panel
ax_norm.set_ylabel("Normal", fontsize=10, fontweight="bold", labelpad=10)
# "Neoplasm" subtitle near the Neoplasm panel
ax_neo.set_ylabel("Neoplasm", fontsize=10, fontweight="bold", labelpad=10)


# ── Barcode region shading (light background spans on both panels) ──────────
for site_name in ["r1", "r2", "r3"]:
    rinfo = REGIONS[site_name]
    ds, de = region_display_ranges[site_name]
    for ax in [ax_norm, ax_neo]:
        ax.axvspan(ds - 0.5, de + 0.5, color=rinfo["color"], alpha=0.08, zorder=0)

# Region labels below x-axis (on Neoplasm panel)
from matplotlib.transforms import blended_transform_factory
trans_label = blended_transform_factory(ax_neo.transData, ax_neo.transAxes)
for site_name in ["r1", "r2", "r3"]:
    rinfo = REGIONS[site_name]
    ds, de = region_display_ranges[site_name]
    mid = (ds + de) / 2
    ax_neo.text(mid, -0.06, rinfo["label"],
                transform=trans_label, ha="center", va="top",
                fontsize=9, fontweight="bold", color=rinfo["color"])


# ── X-axis (relative positions within each barcode, starting at 1) ──────────
x_min_disp = -1
x_max_disp = max(DISPLAY_POS) + 1
ax_norm.set_xlim(x_min_disp, x_max_disp)

# Custom x-ticks: show relative bp positions (1-based within each barcode)
xtick_positions = []
xtick_labels_list = []
for site in ["r1", "r2", "r3"]:
    rinfo = REGIONS[site]
    region_len = rinfo["end"] - rinfo["start"] + 1
    # Show ticks at start (1), center, and end of each barcode
    key_relative = [1, (region_len + 1) // 2, region_len]
    for rel_pos in key_relative:
        real_bp = rinfo["start"] + rel_pos - 1
        if real_bp in pos_to_display:
            xtick_positions.append(pos_to_display[real_bp])
            xtick_labels_list.append(str(rel_pos))

ax_neo.set_xticks(xtick_positions)
ax_neo.set_xticklabels(xtick_labels_list, fontsize=8)

# "Position along lineage barcode construct (bp)" centered under Barcode-2
ds2, de2 = region_display_ranges["r2"]
mid2 = (ds2 + de2) / 2
ax_neo.text(mid2, -0.10, "Position along lineage barcode construct (bp)",
            transform=trans_label, ha="center", va="top",
            fontsize=10, color="black")


# ── Legend — outside plot area on the right ─────────────────────────────────
handles = [mpatches.Patch(facecolor=ETYPE_COLORS[et], label=ETYPE_LABELS[et])
           for et in ["D", "I", "S", "X"]]
ax_neo.legend(handles=handles, loc="center left", fontsize=8,
              title="Edit type", title_fontsize=9, framealpha=0.9,
              bbox_to_anchor=(1.01, 0.5))


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

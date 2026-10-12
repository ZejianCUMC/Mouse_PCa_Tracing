# =============================================================================
# Per-sample barcode heatmap (P3 poster style) Pipeline
# =============================================================================
# Description : Reproduce the ancestorOT poster panel P3_JZ138_barcode_heatmap for every Npp53 and
#               Np53Rb1 mouse: rectangular Greedy tree (leaf x = PhyloVelo pseudotime) + cell-state
#               strip + barcode-edit heatmap (intBC groups sorted by completeness, then diversity;
#               white = uncut, grey = missing). Plotting code/geometry copied from
#               23-ancestorOT/report/gen_jz138_barcode_heatmap.py; helpers rect_layout /
#               allele_colormap / cm_to_rgb imported read-only from gen_barcode_joint.py.
# Input       : Npp53 v8 Greedy trees + per-sample allele tables (matrix rebuilt with the v8 call),
#               Np53Rb1 Greedy trees + the saved character matrices they were solved from,
#               cell states, PhyloVelo pseudotime.
# Output      : <project>/figures/barcode_heatmaps_20260925/*_barcode_heatmap.{png,pdf}, legends,
#               JOINS.csv, rebuilt Npp53 matrices (_rebuilt_character_matrices/)
# Conda env   : stage "npp53_matrices": /home/jrchristin/miniconda3/envs/Cassiopeia_v2/bin/python
#               stage "plot":           /home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python
# Key deps    : cassiopeia 2.0 (matrices), ete3, matplotlib, pandas
# Validated   : 2026-09-25
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import argparse, json, re, sys, warnings
from collections import OrderedDict, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
MB = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
N, R = MB / "Npp53_tracing", MB / "Npp53Rb1_Tracing"
RL = R / "lineage_tracing"
DEFAULTS = {
    "ref_dir": "/media/jrchristin/Backup10TB1/agent/1-project/23-ancestorOT/report",   # read-only helpers
    "cut_sites": ["r1", "r2", "r3"],
    "allele_rep_thresh": 0.99,                    # same as tree building (v8 + Rb1 runs)
    "exclude_states_npp53": ["Ionocyte", "Cycling"],   # v8 tree-build exclusion
    # Npp53
    "npp53_samples": ["JZ136"] + [f"JZ{i}" for i in list(range(201, 215)) + list(range(224, 233)) + [234]],
    "npp53_trees": N / "lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy",
    "npp53_allele": N / "lineage_tracing/allele_tables/per_sample",
    "npp53_obs": N / "data/npp53_merge/obs_metadata.csv",
    "npp53_pt": N / "analysis/evolution_features_v2/carta_fine12/biology_review/evopath_trajectory_features/combined_pseudotime_percell.csv",
    "npp53_pt_fallback": N / "lineage_tracing/cassiopeia_trees_v8_rebuilt/phylovelo_out",
    "npp53_viewer": N / "reports/Npp53_tree_viewer_interactive.html",     # DATA.palettes.state
    "npp53_out": N / "figures/barcode_heatmaps_20260925",
    # Np53Rb1: (sample, cassiopeia run dir, tree-viewer run dir); JZ135 only via its n3000 subsample tree
    "rb1_runs": [(s, RL / "cassiopeia_eight_20260923_seed42", RL / "tree_viewer_20260923_seed42")
                 for s in ["JZ137", "JZ138", "JZ256", "JZ259", "JZ261", "JZ262", "JZ263"]]
                + [("JZ135_n3000", RL / "tree_viewer_20260923_seed42/subsample_cassiopeia", RL / "tree_viewer_20260923_seed42")]
                + [(s, RL / "cassiopeia_mjz026_028_20260925_seed42", RL / "tree_viewer_mjz026_028_20260925_seed42")
                   for s in ["JZ265", "JZ266", "JZ267", "JZ268", "JZ269", "JZ270", "JZ271", "JZ272", "JZ273",
                             "JZ274", "JZ276", "JZ277", "JZ279", "JZ280"]],
    "rb1_palette": R / "analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/exclude_current_c17_c19_seed42/SPLIT_PLOTS_COMPLETE.json",
    "rb1_out": R / "figures/barcode_heatmaps_20260925",
    # plot (reference geometry)
    "figsize": (16, 11), "width_ratios": [0.9, 0.10, 5.4], "dpi": 500,
    "header_font_max": 6.5, "site_font_max": 5.0,
    "header_gap_px": 6,        # min free output-px between neighbouring intBC headers (width-based header rule)
    "min_title_gap_px": 15,    # min output-px between "Cell states" title and heatmap site numbers/headers
    "npp53_matrix_dir": N / "figures/barcode_heatmaps_20260925/_rebuilt_character_matrices",
    "rb1_display_trim": {"JZ137": 0.10, "JZ138": 0.10},   # Zejian 2026-09-25: display-only hide of intBCs with <10 % cell coverage (legacy panels)
}
TREAT = {"Intact": "Intact", "Castration": "Castration", "Castration+ENZA": "CastrationENZA",
         "Castration + Enzalutamide": "CastrationENZA"}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def column_keys(at, cut_sites, thresh):
    """Replicate cassiopeia.pp.convert_alleletable_to_character_matrix column order: '<intBC><site>'
    keys in first-appearance order over cells (insertion order), minus low-diversity keys."""
    filtered = defaultdict(OrderedDict)
    for cell, intbc, *vals in at[["cellBC", "intBC", *cut_sites]].itertuples(index=False):
        for c, v in zip(cut_sites, vals):
            filtered[cell].setdefault(f"{intbc}{c}", []).append(v)
    dist = defaultdict(list)
    for cell in filtered:
        for key in filtered[cell]:
            dist[key].extend(list(set(filtered[cell][key])))
    keep = []
    for key in dist:
        props = np.unique([str(x) for x in dist[key]], return_counts=True)[1] / len(dist[key])
        if not np.any(props > thresh): keep.append(key)
    return keep, filtered

def verify_keys(cm, keys, filtered):
    """Each column: cells lacking that intBC-site key must be exactly the -1 cells."""
    assert len(keys) == cm.shape[1], (len(keys), cm.shape)
    isna = lambda v: not isinstance(v, str) and pd.isna(v)
    for j, k in enumerate(keys):
        miss = set(cm.index[cm.iloc[:, j].astype(str) == "-1"])
        for c in cm.index:
            vals = filtered.get(c, {}).get(k)
            if vals is None: assert c in miss, (k, c, "absent key not -1")                 # intBC absent -> -1
            elif all(isna(v) for v in vals): assert c in miss, (k, c, "blank not -1")      # un-detected site -> -1
            else: assert c not in miss, (k, c, "observed site coded -1")

def to_numeric_cm(cm):
    """Ints stay; ambiguous tuple states '(a, b)' -> NaN (drawn as missing/grey)."""
    n_amb = int(cm.astype(str).apply(lambda s: s.str.startswith("(")).values.sum())
    return cm.apply(pd.to_numeric, errors="coerce"), n_amb

def npp53_matrices(p):
    """Stage 1 (Cassiopeia_v2 env): rebuild each Npp53 v8 matrix with the exact v8 call; save + column keys."""
    import cassiopeia as cas
    from ete3 import Tree
    wd = p["npp53_out"] / "_rebuilt_character_matrices"; wd.mkdir(parents=True, exist_ok=True)
    for s in p["npp53_samples"]:
        f = wd / f"{s}_character_matrix.csv"
        if f.exists(): print("exists", f.name); continue
        at = pd.read_csv(p["npp53_allele"] / f"{s}_allele_table.csv")          # default read, as v8
        at = at[~at["putative_state"].isin(set(p["exclude_states_npp53"]))]
        cm, _, _ = cas.pp.convert_alleletable_to_character_matrix(at, allele_rep_thresh=p["allele_rep_thresh"])
        leaves = set(Tree(open(p["npp53_trees"] / f"{s}_tree_greedy.nwk").read(), format=1).get_leaf_names())
        assert set(cm.index) == leaves, f"{s}: matrix rows != v8 tree leaves"
        keys, filt = column_keys(at, p["cut_sites"], p["allele_rep_thresh"])
        verify_keys(cm, keys, filt)
        cm.to_csv(f); (wd / f"{s}_columns.json").write_text(json.dumps(keys))
        print("built", s, cm.shape, flush=True)

def group_columns(cm, keys):
    """Group columns by intBC (sites in order), sort groups by completeness then diversity (reference rule)."""
    groups = OrderedDict()
    for j, k in enumerate(keys):
        groups.setdefault(k[:-2], []).append((k[-2:], j))
    stats = []
    for g, cols in groups.items():
        sub = cm.iloc[:, [j for _, j in sorted(cols)]]
        v = sub.values
        stats.append((g, int((sub > 0).sum().sum()), int(pd.unique(v[v > 0]).size)))
    order = [g for g, _, _ in sorted(stats, key=lambda s: (-s[1], -s[2]))]
    cols, sites, bounds = [], [], []
    for g in order:
        c = sorted(groups[g]); bounds.append((len(cols), len(cols) + len(c)))
        cols += [j for _, j in c]; sites += [int(s[1]) for s, _ in c]
    return cm.iloc[:, cols], sites, bounds

def plot_panel(tr, states, cm, keys, pt, scolor, title, out_stem, p, note=None, trim_cov=None):
    """Reference plotting code (gen_jz138_barcode_heatmap.py) with intBC groups of 1-3 kept sites."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from gen_barcode_joint import allele_colormap, cm_to_rgb, rect_layout
    tr.ladderize()
    leaves, y, depth = rect_layout(tr, pt)
    order = [lf.name for lf in leaves]
    cm = cm.reindex(order)
    cm, sites, bounds = group_columns(cm, keys)
    n_all, hidden = len(bounds), 0
    if trim_cov is not None:   # DISPLAY-ONLY trim: hide intBC groups detected (>=0, i.e. not -1/missing) in < trim_cov of cells
        keep_cols, keep_sites, keep_bounds = [], [], []
        for a, b in bounds:
            cov = float((cm.iloc[:, a:b] >= 0).any(axis=1).mean())
            if cov < trim_cov: hidden += 1; continue
            keep_bounds.append((len(keep_cols), len(keep_cols) + (b - a)))
            keep_cols += list(range(a, b)); keep_sites += sites[a:b]
        cm, sites, bounds = cm.iloc[:, keep_cols], keep_sites, keep_bounds
        tnote = f"Display: {hidden} of {n_all} intBCs with <{trim_cov:.0%} cell coverage hidden (tree built from all {n_all})"
        note = tnote if not note else f"{note} · {tnote}"
    amap = allele_colormap(cm); rgb = cm_to_rgb(cm, amap); Nc = len(order)
    fig = plt.figure(figsize=p["figsize"])
    gs = fig.add_gridspec(1, 3, width_ratios=p["width_ratios"], wspace=0.03)
    axT = fig.add_subplot(gs[0, 0]); segs = []
    for nd in tr.traverse():
        if nd.is_root(): continue
        segs.append([[depth[nd.up], y[nd]], [depth[nd], y[nd]]])
    for nd in tr.traverse():
        if nd.is_leaf(): continue
        ys = [y[c] for c in nd.children]; segs.append([[depth[nd], min(ys)], [depth[nd], max(ys)]])
    axT.add_collection(LineCollection(segs, colors="0.15", linewidths=0.30))
    axT.set_xlim(0, 1.02); axT.set_ylim(-1, Nc); axT.invert_yaxis()
    axT.set_yticks([]); axT.set_xticks([0, 0.5, 1.0]); axT.tick_params(labelsize=9)
    axT.set_xlabel("PhyloVelo pseudotime", fontsize=9)
    for s in ["top", "right", "left"]: axT.spines[s].set_visible(False)
    axT.set_title("Phylogenetic\ntree", fontsize=11)
    axS = fig.add_subplot(gs[0, 1], sharey=axT)
    strip = np.array([[matplotlib.colors.to_rgb(scolor(s)) for s in [states.get(c, "NA") for c in order]]]).transpose(1, 0, 2)
    axS.imshow(strip, aspect="auto", extent=[0, 1, Nc - 0.5, -0.5]); axS.set_xticks([]); axS.set_yticks([])
    axS.set_title("Cell\nstates", fontsize=11)
    axH = fig.add_subplot(gs[0, 2], sharey=axT)
    axH.imshow(rgb, aspect="auto", extent=[0, cm.shape[1], Nc - 0.5, -0.5], interpolation="none")
    for a, _ in bounds[1:]: axH.axvline(a, color="0.75", lw=0.5)
    nbc = len(bounds); ncol = cm.shape[1]
    hfs = max(3.5, min(p["header_font_max"], p["header_font_max"] * 26 / nbc))
    sfs = max(3.0, min(p["site_font_max"], p["site_font_max"] * 78 / ncol))
    axH.set_xticks([j + 0.5 for j in range(ncol)]); axH.set_xticklabels([str(s) for s in sites], fontsize=sfs)
    axH.xaxis.set_ticks_position("top"); axH.tick_params(top=False, labeltop=True, labelbottom=False, pad=1)
    tb = axH.get_xaxis_transform()
    axH.set_yticks([]); axH.set_xticks([], minor=True); axH.set_xlabel("")
    axH.text(0.5, 1.075, title, transform=axH.transAxes, ha="center", va="bottom", fontsize=16, fontweight="bold")
    if note: axH.text(1.0, -0.012, note, transform=axH.transAxes, ha="right", va="top", fontsize=7, color="0.35")
    for s in ["top", "right", "left"]: axH.spines[s].set_visible(False)
    fig.tight_layout()
    # Header labels chosen by AVAILABLE PIXEL WIDTH per intBC group (same rule on every panel):
    # "intBC\nN" if it fits, else bare "N", else bare "N" rotated 90°.
    fig.canvas.draw(); rdr = fig.canvas.get_renderer()
    hbox = axH.get_window_extent(rdr); px_per_col = hbox.width / ncol; gap = p["header_gap_px"] * fig.dpi / p["dpi"]
    heads = []
    for i, (a, b) in enumerate(bounds):
        avail = (b - a) * px_per_col - gap
        for txt, rot in [(f"intBC\n{i+1}", 0), (f"{i+1}", 0), (f"{i+1}", 90)]:
            t = axH.text((a + b) / 2, 1.03, txt, transform=tb, ha="center", va="bottom", fontsize=hfs, fontweight="bold", rotation=rot)
            if t.get_window_extent(rdr).width <= avail or rot == 90: break
            t.remove()
        heads.append(t)
    # "Cell states" title: keep >= min_title_gap_px (output px) from the site numbers / headers to its right.
    fig.canvas.draw(); ttl = axS.title; tbb = ttl.get_window_extent(rdr)
    others = [lab.get_window_extent(rdr) for lab in axH.get_xticklabels() if lab.get_text()] + [h.get_window_extent(rdr) for h in heads]
    right_neigh = [o for o in others if o.y1 > tbb.y0 and o.y0 < tbb.y1 and o.x0 >= tbb.x0]
    need = p["min_title_gap_px"] * fig.dpi / p["dpi"]
    if right_neigh:
        limit = min(o.x0 for o in right_neigh) - need
        if tbb.x1 > limit:
            ttl.set_ha("right"); ttl.set_x(axS.transAxes.inverted().transform((limit, 0))[0])
    for ext in ["png", "pdf"]:
        f = Path(f"{out_stem}.{ext}"); assert not f.exists(), f"refusing to overwrite {f}"
        fig.savefig(f, dpi=p["dpi"], bbox_inches="tight")
    plt.close(fig)
    return nbc, ncol, n_all, hidden

def legend(pal, present, title, out_stem):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    keys = [k for k in pal if k in present]
    fig = plt.figure(figsize=(3.6, 0.32 * len(keys) + 0.8))
    h = [Line2D([0], [0], marker="s", color="w", markerfacecolor=pal[k], markersize=11, label=k) for k in keys]
    fig.legend(handles=h, title=title, loc="center", frameon=False, fontsize=9, title_fontsize=10)
    for ext in ["png", "pdf"]:
        f = Path(f"{out_stem}.{ext}"); assert not f.exists(), f"refusing to overwrite {f}"
        fig.savefig(f, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    sys.path.insert(0, p["ref_dir"])
    from ete3 import Tree
    joins = []
    # ---------- Npp53 ----------
    out = p["npp53_out"]; out.mkdir(parents=True, exist_ok=True); wd = p["npp53_matrix_dir"]
    html = p["npp53_viewer"].read_text()
    pal53 = json.loads(re.search(r"const DATA = (\{.*?\});\n", html, re.S).group(1))["palettes"]["state"]
    norm = lambda s: re.sub(r"^c\d+_", "", s).replace("Lum-", "Luminal-")
    name53 = {norm(k): k for k in pal53}                     # 'Luminal-Distal' -> 'c1_Lum-Distal'
    obs = pd.read_csv(p["npp53_obs"], index_col=0)
    pt_all = pd.read_csv(p["npp53_pt"])
    present53 = set()
    for s in p["npp53_samples"]:
        tr = Tree(open(p["npp53_trees"] / f"{s}_tree_greedy.nwk").read(), format=1)
        leaves = set(tr.get_leaf_names())
        cm, n_amb = to_numeric_cm(pd.read_csv(wd / f"{s}_character_matrix.csv", index_col=0))
        keys = json.loads((wd / f"{s}_columns.json").read_text())
        o = obs[obs.sample_id == s].copy(); o["bc"] = o.index.str.replace(r"^MJZ\d+_", "", regex=True).str.replace(r"-\d+$", "", regex=True)
        states = {b: name53[st] for b, st in zip(o.bc, o.putative_state)}
        pts = pt_all[pt_all["sample"] == s].set_index("cellBC")["pt_raw"]; src = "combined_pseudotime_percell (best run per sample)"
        if pts.empty:
            fb = pd.read_csv(p["npp53_pt_fallback"] / f"{s}_phylovelo_pseudotime.csv", index_col=0).iloc[:, 0]; pts, src = fb, "phylovelo_out fallback"
        pts = pts[~pts.index.duplicated()]; pts = (pts - pts.min()) / (pts.max() - pts.min() + 1e-9)
        tr_ = obs.loc[obs.sample_id == s, "treatment"].iloc[0]
        stem = out / f"Npp53_{TREAT[tr_]}_{s}_barcode_heatmap"
        present53 |= {states[l] for l in leaves if l in states}
        n_nopt = int(sum(l not in pts.index for l in leaves))
        note53 = f"{n_nopt:,} of {len(leaves):,} cells without PhyloVelo pseudotime placed at median" if n_nopt else None
        nbc, ncol, n_all, hidden = plot_panel(tr, states, cm, keys, pts.to_dict(), lambda x: pal53.get(x, "#d9d9d9"),
                               f"{s} Barcode edits", stem, p, note53)
        joins.append(dict(project="Npp53", sample=s, leaves=len(leaves), matrix_rows=int(len(set(cm.index) & leaves)),
                          states_matched=int(sum(l in states for l in leaves)), pseudotime_matched=int(sum(l in pts.index for l in leaves)),
                          intBC_groups=n_all, intBC_groups_shown=nbc, display_trim_threshold=None, intBC_hidden=hidden, char_cols_shown=ncol, ambiguous_entries=n_amb, pseudotime_source=src, file=stem.name))
        print(joins[-1], flush=True)
    legend(pal53, present53, "Npp53 cell state", out / "Npp53_cell_state_legend")
    pd.DataFrame([j for j in joins if j["project"] == "Npp53"]).to_csv(out / "JOINS.csv", index=False)
    # ---------- Np53Rb1 ----------
    out = p["rb1_out"]; out.mkdir(parents=True, exist_ok=True)
    palr = {k: "#" + "".join(f"{round(v*255):02x}" for v in rgb) for k, rgb in json.loads(p["rb1_palette"].read_text())["colour_mapping"].items()}
    presentr = set()
    for s, cas_dir, tv_dir in p["rb1_runs"]:
        res = cas_dir / "results" / s
        tr = Tree(open(res / f"{s}_tree_greedy.nwk").read(), format=1); leaves = set(tr.get_leaf_names())
        cm_raw = pd.read_csv(res / "character_matrix.csv.gz", index_col=0, dtype=str)
        at = pd.read_csv(cas_dir / "inputs" / s / "allele_table.csv.gz", dtype={"cellBC": str, "intBC": str})
        keys, filt = column_keys(at, p["cut_sites"], p["allele_rep_thresh"])
        verify_keys(cm_raw, keys, filt)
        cm, n_amb = to_numeric_cm(cm_raw)
        meta = pd.read_csv(cas_dir / "inputs" / s / "cell_metadata.csv", index_col=0)
        states = meta.putative_state.to_dict()
        ptdir = tv_dir / "results" / s / ("PhyloVelo_excl" if (tv_dir / "results" / s / "PhyloVelo_excl" / "pseudotime.csv").exists() else "PhyloVelo")
        pts = pd.read_csv(ptdir / "pseudotime.csv", index_col=0)["pseudotime"]
        pts = (pts - pts.min()) / (pts.max() - pts.min() + 1e-9)
        tok = TREAT[meta.treatment.iloc[0]]
        title = "JZ135 (n3000) Barcode edits" if s == "JZ135_n3000" else f"{s} Barcode edits"
        note = None
        if n_amb: note = f"{n_amb:,} ambiguous (multi-allele) entries shown as missing (grey)"
        if ptdir.name == "PhyloVelo_excl": note = "1 isolated cell excluded from pseudotime (placed at the median)"
        stem = out / f"Np53Rb1_{tok}_{s}_barcode_heatmap"
        presentr |= {states[l] for l in leaves if l in states}
        trim = p["rb1_display_trim"].get(s)
        nbc, ncol, n_all, hidden = plot_panel(tr, states, cm, keys, pts.to_dict(), lambda x: palr.get(x, "#d9d9d9"), title, stem, p, note, trim_cov=trim)
        joins.append(dict(project="Np53Rb1", sample=s, leaves=len(leaves), matrix_rows=int(len(set(cm.index) & leaves)),
                          states_matched=int(sum(l in states for l in leaves)), pseudotime_matched=int(sum(l in pts.index for l in leaves)),
                          intBC_groups=n_all, intBC_groups_shown=nbc, display_trim_threshold=trim, intBC_hidden=hidden, char_cols_shown=ncol, ambiguous_entries=n_amb, pseudotime_source=f"{tv_dir.name}/{ptdir.name}", file=stem.name))
        print(joins[-1], flush=True)
    legend(palr, presentr, "Np53Rb1 cell state", out / "Np53Rb1_cell_state_legend")
    pd.DataFrame([j for j in joins if j["project"] == "Np53Rb1"]).to_csv(out / "JOINS.csv", index=False)
    return {"output_path": [str(p["npp53_out"]), str(p["rb1_out"])], "params_used": {k: str(v) for k, v in p.items()}, "summary": joins}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", choices=["npp53_matrices", "plot"], required=True)
    a = ap.parse_args()
    if a.stage == "npp53_matrices": npp53_matrices(DEFAULTS)
    else: run_pipeline()

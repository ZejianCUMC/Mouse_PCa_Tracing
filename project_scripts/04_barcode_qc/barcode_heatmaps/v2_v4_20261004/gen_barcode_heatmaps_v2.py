# =============================================================================
# Per-sample barcode heatmap v2 Pipeline
# =============================================================================
# Description : Revision of 2026-09-25 gen_barcode_heatmaps.py (same trees, cell states, pseudotime,
#               colours, intBC ordering) with Zejian's 2026-10-04 changes:
#               1) cell-state legend drawn inside each figure, on the right;
#               2) display-only rule on EVERY panel: hide intBCs detected in < 10 % of tree cells;
#               3) only intBCs used by the tree (>=1 site kept at 0.99); every kept intBC shows ALL its cut sites (r1/r2/r3) - display matrix rebuilt from the
#                  same allele table with allele_rep_thresh=1.0 (keeps low-diversity sites the tree build
#                  dropped at 0.99). Trees are unchanged (still built from the 0.99 matrix).
# Input       : as v1 (DEFAULTS imported from ../barcode_heatmaps_20260925/gen_barcode_heatmaps.py)
# Output      : <project>/figures/barcode_heatmaps_v2_20261004/*_barcode_heatmap.{png,pdf}, JOINS.csv,
#               _display_matrices/ (full-site matrices + column keys)
# Conda env   : stage "matrices": /home/jrchristin/miniconda3/envs/Cassiopeia_v2/bin/python
#               stage "plot":     /home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python
# Key deps    : cassiopeia 2.0 (matrices), ete3, matplotlib, pandas
# Validated   : 2026-10-04
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import argparse, json, re, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "barcode_heatmaps_20260925"))
import gen_barcode_heatmaps as v1   # read-only reuse: DEFAULTS, column_keys, verify_keys, to_numeric_cm, group_columns

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
DEFAULTS = {
    **v1.DEFAULTS,
    "display_allele_rep_thresh": 1.0,   # keep every cut site for display (tree build used 0.99)
    "min_intbc_coverage": 0.10,         # Zejian 2026-10-04: hide intBCs detected (any site != -1) in < 10 % of tree cells
    "npp53_out": v1.N / "figures/barcode_heatmaps_v2_20261004",
    "rb1_out": v1.R / "figures/barcode_heatmaps_v2_20261004",
    "figsize": (17.5, 11),
    "width_ratios": [0.9, 0.10, 5.4, 1.0],   # tree | state strip | heatmap | legend
    "legend_font": 9, "legend_title_font": 10,
}

TITLE_TREAT = {"Intact": "Intact", "Castration": "Castration", "CastrationENZA": "Castration+ENZA"}   # panel title: "<model> <treatment> (<JZ>)"

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def samples(p):
    """(project, sample, tree path, allele table loader, convert kwargs, output dir)."""
    out = []
    for s in p["npp53_samples"]:
        def load(s=s):   # identical to v8 tree build input (default read_csv, Ionocyte/Cycling removed)
            at = pd.read_csv(p.get("npp53_allele_override", {}).get(s, p["npp53_allele"] / f"{s}_allele_table.csv"))   # optional override (v3 rebuilt inputs)
            return at[~at["putative_state"].isin(set(p["exclude_states_npp53"]))]
        out.append(("Npp53", s, p.get("npp53_tree_override", {}).get(s, p["npp53_trees"] / f"{s}_tree_greedy.nwk"), load, {}, p["npp53_out"]))
    for s, cas_dir, tv_dir in p["rb1_runs"]:
        def load(cas_dir=cas_dir, s=s):   # identical to run_cassiopeia_eight.py input
            return pd.read_csv(cas_dir / "inputs" / s / "allele_table.csv.gz", dtype={"cellBC": str, "intBC": str})
        out.append(("Np53Rb1", s, cas_dir / "results" / s / f"{s}_tree_greedy.nwk", load,
                    dict(collapse_duplicates=True, cut_sites=p["cut_sites"]), p["rb1_out"]))
    return out

def build_matrices(p):
    """Stage 1 (Cassiopeia_v2 env): full-site display matrix per sample, rows restricted to tree leaves."""
    import cassiopeia as cas
    from ete3 import Tree
    for proj, s, tree_f, load, kw, out in samples(p):
        wd = out / "_display_matrices"; wd.mkdir(parents=True, exist_ok=True)
        f = wd / f"{s}_display_matrix.csv"
        assert not f.exists(), f"refusing to overwrite {f}"
        at = load()
        cm, _, _ = cas.pp.convert_alleletable_to_character_matrix(at, allele_rep_thresh=p["display_allele_rep_thresh"], **kw)
        leaves = set(Tree(open(tree_f).read(), format=1).get_leaf_names())
        assert leaves <= set(cm.index), f"{s}: tree leaves missing from display matrix"
        keys, filt = v1.column_keys(at, p["cut_sites"], p["display_allele_rep_thresh"])
        v1.verify_keys(cm, keys, filt)
        tree_keys, _ = v1.column_keys(at, p["cut_sites"], p["allele_rep_thresh"])   # columns the tree was built from (0.99)
        tree_intbcs = sorted({k[:-2] for k in tree_keys})
        keep = [j for j, k in enumerate(keys) if k[:-2] in set(tree_intbcs)]     # only intBCs used by the tree, all their sites
        cm, keys = cm.iloc[:, keep], [keys[j] for j in keep]
        cm = cm.loc[[c for c in cm.index if c in leaves]]
        cm.to_csv(f); (wd / f"{s}_columns.json").write_text(json.dumps(keys))
        print("built", proj, s, cm.shape, flush=True)

def plot_panel(tr, states, cm, keys, pt, pal, title, out_stem, p, note=None):
    """v1 plot_panel geometry + coverage trim on all panels + in-figure cell-state legend (right)."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    from gen_barcode_joint import allele_colormap, cm_to_rgb, rect_layout
    scolor = lambda x: pal.get(x, "#d9d9d9")
    tr.ladderize()
    leaves, y, depth = rect_layout(tr, pt)
    order = [lf.name for lf in leaves]
    cm = cm.reindex(order)
    cm, sites, bounds = v1.group_columns(cm, keys)
    n_all, hidden, thr = len(bounds), 0, p["min_intbc_coverage"]
    keep_cols, keep_sites, keep_bounds = [], [], []
    for a, b in bounds:   # DISPLAY-ONLY: hide intBC groups detected (any site != -1) in < thr of cells
        if float((cm.iloc[:, a:b] >= 0).any(axis=1).mean()) < thr: hidden += 1; continue
        keep_bounds.append((len(keep_cols), len(keep_cols) + (b - a)))
        keep_cols += list(range(a, b)); keep_sites += sites[a:b]
    cm, sites, bounds = cm.iloc[:, keep_cols], keep_sites, keep_bounds
    tnote = f"Display: {hidden} of {n_all} intBCs detected in <{thr:.0%} of cells hidden" if hidden else f"All {n_all} intBCs detected in ≥{thr:.0%} of cells"
    note = tnote if not note else f"{note} · {tnote}"
    amap = allele_colormap(cm); rgb = cm_to_rgb(cm, amap); Nc = len(order)
    fig = plt.figure(figsize=p["figsize"])
    gs = fig.add_gridspec(1, 4, width_ratios=p["width_ratios"], wspace=0.03)
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
    st = [states.get(c, "NA") for c in order]
    strip = np.array([[matplotlib.colors.to_rgb(scolor(s)) for s in st]]).transpose(1, 0, 2)
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
    axH.text(1.0, -0.012, note, transform=axH.transAxes, ha="right", va="top", fontsize=7, color="0.35")
    for s in ["top", "right", "left"]: axH.spines[s].set_visible(False)
    # legend (right): cell states present in this mouse, palette order, with cell counts
    axL = fig.add_subplot(gs[0, 3]); axL.axis("off")
    cnt = pd.Series(st).value_counts()
    ks = [k for k in pal if k in cnt.index] + [k for k in cnt.index if k not in pal]
    h = [Line2D([0], [0], marker="s", color="w", markerfacecolor=scolor(k), markersize=10,
                label=f"{k} ({cnt[k]:,})") for k in ks]
    lg = axL.legend(handles=h, title="Cell state", loc="center left", bbox_to_anchor=(0.02, 0.5), frameon=False,
                    fontsize=p["legend_font"], title_fontsize=p["legend_title_font"], handletextpad=0.4)
    lg._legend_box.align = "left"   # matplotlib 3.x in hotspot env has no alignment= kwarg
    fig.tight_layout()
    # header labels by available pixel width per intBC group (v1 rule)
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

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    sys.path.insert(0, p["ref_dir"])
    from ete3 import Tree
    # palettes / states / pseudotime exactly as v1
    html = p["npp53_viewer"].read_text()
    pal53 = json.loads(re.search(r"const DATA = (\{.*?\});\n", html, re.S).group(1))["palettes"]["state"]
    norm = lambda s: re.sub(r"^c\d+_", "", s).replace("Lum-", "Luminal-")
    name53 = {norm(k): k for k in pal53}
    obs = pd.read_csv(p["npp53_obs"], index_col=0); pt_all = pd.read_csv(p["npp53_pt"])
    palr = {k: "#" + "".join(f"{round(v*255):02x}" for v in rgb) for k, rgb in json.loads(p["rb1_palette"].read_text())["colour_mapping"].items()}
    rb1 = {s: (cd, tv) for s, cd, tv in p["rb1_runs"]}
    joins = {"Npp53": [], "Np53Rb1": []}
    for proj, s, tree_f, _, _, out in samples(p):
        tr = Tree(open(tree_f).read(), format=1); leaves = set(tr.get_leaf_names())
        wd = p.get("display_matrix_dirs", {}).get(proj, out / "_display_matrices")   # optional: reuse matrices from another run (v3)
        cm, n_amb = v1.to_numeric_cm(pd.read_csv(wd / f"{s}_display_matrix.csv", index_col=0))
        keys = json.loads((wd / f"{s}_columns.json").read_text())
        assert set(cm.index) == leaves
        if proj == "Npp53":
            o = obs[obs.sample_id == s].copy(); o["bc"] = o.index.str.replace(r"^MJZ\d+_", "", regex=True).str.replace(r"-\d+$", "", regex=True)
            states = {b: name53[st] for b, st in zip(o.bc, o.putative_state)}
            pts = pt_all[pt_all["sample"] == s].set_index("cellBC")["pt_raw"]; src = "combined_pseudotime_percell (best run per sample)"
            if s in p.get("npp53_pt_override", {}):   # v3: rebuilt PhyloVelo (UMAP) pseudotime
                pts = pd.read_csv(p["npp53_pt_override"][s], index_col=0).iloc[:, 0]; src = str(p["npp53_pt_override"][s])
            if pts.empty:
                pts = pd.read_csv(p["npp53_pt_fallback"] / f"{s}_phylovelo_pseudotime.csv", index_col=0).iloc[:, 0]; src = "phylovelo_out fallback"
            pts = pts[~pts.index.duplicated()]; pts = (pts - pts.min()) / (pts.max() - pts.min() + 1e-9)
            tok = v1.TREAT[obs.loc[obs.sample_id == s, "treatment"].iloc[0]]
            n_nopt = int(sum(l not in pts.index for l in leaves))
            note = f"{n_nopt:,} of {len(leaves):,} cells without PhyloVelo pseudotime placed at median" if n_nopt else None
            title, pal = f"Npp53 {TITLE_TREAT[tok]} ({s})", pal53
        else:
            cas_dir, tv_dir = rb1[s]
            meta = pd.read_csv(cas_dir / "inputs" / s / "cell_metadata.csv", index_col=0)
            states = meta.putative_state.to_dict()
            ptdir = tv_dir / "results" / s / ("PhyloVelo_excl" if (tv_dir / "results" / s / "PhyloVelo_excl" / "pseudotime.csv").exists() else "PhyloVelo")
            pts = pd.read_csv(ptdir / "pseudotime.csv", index_col=0)["pseudotime"]
            pts = (pts - pts.min()) / (pts.max() - pts.min() + 1e-9); src = f"{tv_dir.name}/{ptdir.name}"
            tok = v1.TREAT[meta.treatment.iloc[0]]
            note = f"{n_amb:,} ambiguous (multi-allele) entries shown as missing (grey)" if n_amb else None
            if ptdir.name == "PhyloVelo_excl": note = "1 isolated cell excluded from pseudotime (placed at the median)"
            title, pal = (f"Np53Rb1 {TITLE_TREAT[tok]} (JZ135, n=3000 subsample)" if s == "JZ135_n3000" else f"Np53Rb1 {TITLE_TREAT[tok]} ({s})"), palr
        stem = out / f"{proj}_{tok}_{s}_barcode_heatmap"
        nbc, ncol, n_all, hidden = plot_panel(tr, states, cm, keys, pts.to_dict(), pal, title, stem, p, note)
        joins[proj].append(dict(project=proj, sample=s, leaves=len(leaves), states_matched=int(sum(l in states for l in leaves)),
                                pseudotime_matched=int(sum(l in pts.index for l in leaves)), intBC_groups=n_all,
                                intBC_hidden=hidden, intBC_shown=nbc, min_intbc_coverage=p["min_intbc_coverage"],
                                site_cols_shown=ncol, ambiguous_entries=n_amb, pseudotime_source=src, file=stem.name))
        print(joins[proj][-1], flush=True)
    for proj, out in [("Npp53", p["npp53_out"]), ("Np53Rb1", p["rb1_out"])]:
        pd.DataFrame(joins[proj]).to_csv(out / "JOINS.csv", index=False)
    return {"output_path": [str(p["npp53_out"]), str(p["rb1_out"])], "params_used": {k: str(v) for k, v in p.items()}, "summary": joins}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", choices=["matrices", "plot"], required=True)
    ap.add_argument("--only", nargs="*", help="restrict to these samples (test run)")
    ap.add_argument("--out_root", help="test run: write both projects under this folder instead")
    a = ap.parse_args()
    if a.out_root:
        DEFAULTS["npp53_out"] = Path(a.out_root) / "Npp53"; DEFAULTS["rb1_out"] = Path(a.out_root) / "Np53Rb1"
    if a.only:
        DEFAULTS["npp53_samples"] = [s for s in DEFAULTS["npp53_samples"] if s in a.only]
        DEFAULTS["rb1_runs"] = [r for r in DEFAULTS["rb1_runs"] if r[0] in a.only]
    if a.stage == "matrices": build_matrices(DEFAULTS)
    else: run_pipeline()

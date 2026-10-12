# =============================================================================
# Reviewed-10 relabel + UMAP / dotplot redraw Pipeline
# =============================================================================
# Description : Zejian 2026-10-08 — rename c6 -> Partial-EMT(Lum), c9 -> Partial-EMT(Mes), c7 -> Mes-Fibroblast-like; re-index so
#               the two partial-EMT states are adjacent (old c9 <-> old c7); other names unchanged. Redraw the UMAP in the
#               numbered-cluster style (number on each cluster; "n ● name" legend on the right; L-shaped UMAP axes) and the
#               (numbers drawn in white circles, like ①, on clusters and in the legend) and the
#               marker dotplot in the poster Panel-D style (23-ancestorOT/report/gen_umap_dotplot_joint.py: BuPu, per-gene z of
#               group means clipped ±2, dot size = % expressed, bold coloured "● name" row labels, legends on the right).
#               Cells, counts, UMAP coordinates and the approved 79-gene panel are unchanged; only labels, order and style.
#               Non-destructive: reads reviewed10/, writes only into this folder.
# Input       : ../named_graph_PRIVATE.h5ad, ../final_metadata_PRIVATE.csv, ../dotplot_counts_PRIVATE.mtx.gz, ../dotplot_panels.csv
# Output      : STATE_MAP.csv, cell_labels_PRIVATE.csv, umap_states_numbered.{png,pdf}, marker_dotplot_panelD.{png,pdf},
#               dotplot_values.csv, RELABEL_RECEIPT.json
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/moscot/bin/python
# Validated   : 2026-10-08
# Parameters  : See DEFAULTS dict below
# =============================================================================
import json, hashlib, numpy as np, pandas as pd, anndata as ad, scipy.io as sio
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
N = Path(__file__).resolve().parent; R = N.parent
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial"], "pdf.fonttype": 42, "ps.fonttype": 42})
DEFAULTS = {
    # old final_cluster -> (new cluster, new name, colour). Re-index: old c9 <-> old c7; renames approved by Zejian 2026-10-08.
    "map": {"c1": ("c1", "Lum-general", "#2166ac"), "c2": ("c2", "Lum-Proximal-like", "#4393c3"), "c3": ("c3", "Lum-basal-hybrid", "#35978f"),
            "c4": ("c4", "Basal-general", "#1b7837"), "c5": ("c5", "Basal-Squamous", "#7fbc41"), "c6": ("c6", "Partial-EMT(Lum)", "#E69F00"),
            "c9": ("c7", "Partial-EMT(Mes)", "#d95f02"), "c8": ("c8", "Mes-ECMremodeling", "#9970ab"), "c7": ("c9", "Mes-Fibroblast-like", "#b4232a"),
            "c10": ("c10", "Mes-Proliferation", "#762a83")},
    # gene blocks re-ordered to follow the new state order (approved 79 genes unchanged)
    "block_order": ["Luminal", "ISG", "Proximal-like", "Basal / epithelial", "Squamous", "EMT", "c9 markers / metabolism", "ECM remodeling",
                    "Matrix / CAF-like", "Contractile", "Cycling", "Plasticity / NE"],
    "clip": 2.0, "point_size_umap": 1.2, "dot_min": 6, "dot_max": 210,
}
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run_pipeline(p=DEFAULTS):
    mp = pd.DataFrame([(o, n, nm, col) for o, (n, nm, col) in p["map"].items()], columns=["old_cluster", "new_cluster", "new_name", "colour"])
    old_names = pd.read_csv(R / "state_summary.csv").set_index("cluster").putative_name
    mp["old_name"] = mp.old_cluster.map(old_names); mp["new_state"] = mp.new_cluster + "_" + mp.new_name
    mp["order"] = mp.new_cluster.str[1:].astype(int); mp = mp.sort_values("order"); mp.to_csv(N / "STATE_MAP.csv", index=False)
    a = ad.read_h5ad(R / "named_graph_PRIVATE.h5ad"); meta = pd.read_csv(R / "final_metadata_PRIVATE.csv", index_col=0)
    assert meta.index.equals(a.obs_names) and set(meta.final_cluster) == set(mp.old_cluster)
    m = mp.set_index("old_cluster"); lab = pd.DataFrame({"old_cluster": meta.final_cluster, "old_state": meta.putative_state,
                                                       "new_cluster": meta.final_cluster.map(m.new_cluster), "new_state": meta.final_cluster.map(m.new_state)})
    lab.to_csv(N / "cell_labels_PRIVATE.csv")
    order = mp.new_state.tolist(); col = dict(zip(mp.new_state, mp.colour)); num = {s: s.split("_")[0][1:] for s in order}
    # ── UMAP: numbered clusters + "n ● name" legend, L-shaped axes ─────────────────────────────────────────────
    CIRC = dict(boxstyle="circle,pad=0.1", facecolor="white", edgecolor="black", linewidth=0.8)   # circled numbers, like ①
    xy = a.obsm["X_umap"]; st = lab.new_state.values; rng = np.random.default_rng(42).permutation(len(xy))
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    ax.scatter(xy[rng, 0], xy[rng, 1], c=[col[s] for s in st[rng]], s=p["point_size_umap"], linewidths=0, rasterized=True)
    for s in order:
        cx, cy = np.median(xy[st == s], axis=0); ax.text(cx, cy, num[s], fontsize=13, ha="center", va="center", color="black", bbox=CIRC)
    ax.set_xticks([]); ax.set_yticks([]); [ax.spines[k].set_visible(False) for k in ax.spines]
    x0, x1 = xy[:, 0].min(), xy[:, 0].max(); y0, y1 = xy[:, 1].min(), xy[:, 1].max(); L = 0.16
    ax.set_xlim(x0 - .06 * (x1 - x0), x1 + .02 * (x1 - x0)); ax.set_ylim(y0 - .06 * (y1 - y0), y1 + .02 * (y1 - y0))
    ox, oy = x0 - .05 * (x1 - x0), y0 - .05 * (y1 - y0)
    ax.plot([ox, ox + L * (x1 - x0)], [oy, oy], color="black", lw=1.6, clip_on=False); ax.plot([ox, ox], [oy, oy + L * (y1 - y0)], color="black", lw=1.6, clip_on=False)
    ax.text(ox + L * (x1 - x0) / 2, oy - .035 * (y1 - y0), "UMAP1", ha="center", va="top", fontsize=15)
    ax.text(ox - .02 * (x1 - x0), oy + L * (y1 - y0) / 2, "UMAP2", ha="right", va="center", rotation=90, fontsize=15)
    lx = 1.02
    for i, s in enumerate(order):
        y = .92 - i * .085
        ax.text(lx - .006, y, num[s], transform=ax.transAxes, ha="center", va="center", fontsize=13, bbox=CIRC)
        ax.scatter([lx + .030], [y], transform=ax.transAxes, s=110, c=col[s], clip_on=False, linewidths=0)
        ax.text(lx + .056, y, s.split("_", 1)[1], transform=ax.transAxes, ha="left", va="center", fontsize=15)
    fig.subplots_adjust(left=.06, right=.66, top=.97, bottom=.05)
    for ext in ["png", "pdf"]: fig.savefig(N / f"umap_states_numbered.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    # ── Dotplot, poster Panel-D style ───────────────────────────────────────────────────────────────────────
    C = sio.mmread(str(R / "dotplot_counts_PRIVATE.mtx.gz")).tocsr()                 # genes x cells
    panels = pd.read_csv(R / "dotplot_panels.csv"); genes0 = pd.read_csv(R / "dotplot_genes.csv").gene.tolist()
    assert panels.gene.tolist() == genes0 and C.shape == (len(genes0), len(meta))
    rank = {b: i for i, b in enumerate(p["block_order"])}; assert set(panels.panel) == set(rank)
    panels["k"] = panels.panel.map(rank); panels["pos"] = range(len(panels)); genes = panels.sort_values(["k", "pos"]).gene.tolist()
    gi = [genes0.index(g) for g in genes]
    X = C[gi].T.multiply(1e4 / meta.total_counts.values[:, None]).tocsr()             # CP10k (cells x genes)
    rows = []
    for s in order:
        mk = (st == s); sub = X[mk]
        rows.append(dict(state=s, **{f"mean|{g}": v for g, v in zip(genes, np.asarray(sub.mean(0)).ravel())}, **{f"pct|{g}": v for g, v in zip(genes, 100 * np.asarray((sub > 0).mean(0)).ravel())}))
    T = pd.DataFrame(rows).set_index("state"); mean = T[[f"mean|{g}" for g in genes]]; mean.columns = genes; pct = T[[f"pct|{g}" for g in genes]]; pct.columns = genes
    lm = np.log1p(mean); z = ((lm - lm.mean(0)) / lm.std(0, ddof=1)).clip(-p["clip"], p["clip"])   # Seurat DotPlot scaling (log1p of mean, scale)
    out = pd.concat({"avg_exp_cp10k": mean, "pct_exp": pct, "scaled": z}, axis=1); out.to_csv(N / "dotplot_values.csv")
    ng, nr = len(genes), len(order); fig, ax = plt.subplots(figsize=(0.255 * ng + 4.2, 0.52 * nr + 2.2))
    norm = Normalize(-p["clip"], p["clip"]); cmap = matplotlib.colormaps["BuPu"]; s_of = lambda q: p["dot_min"] + q / 100 * (p["dot_max"] - p["dot_min"])
    for yi, s in enumerate(order[::-1]):
        ax.scatter(range(ng), [yi] * ng, s=[s_of(v) for v in pct.loc[s]], c=[cmap(norm(v)) for v in z.loc[s]], linewidths=.3, edgecolors="0.55")
    ax.set_xticks(range(ng)); ax.set_xticklabels(genes, rotation=90, fontsize=12, style="italic")
    ax.set_yticks(range(nr)); ax.set_yticklabels([f"●  {s.split('_')[0]} {s.split('_', 1)[1]}" for s in order[::-1]], fontsize=14)
    for tk, s in zip(ax.get_yticklabels(), order[::-1]): tk.set_color(col[s]); tk.set_fontweight("bold")
    ax.set_xlim(-.7, ng - .3); ax.set_ylim(-.7, nr - .3); [ax.spines[k].set_visible(False) for k in ("top", "right")]; ax.tick_params(length=0)
    for q in [25, 50, 75, 100]: ax.scatter([], [], s=s_of(q), c="0.5", linewidths=.3, edgecolors="0.55", label=f"{q}")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.03), frameon=False, title="Percent\nexpressed", fontsize=11, title_fontsize=12, labelspacing=1.0)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([]); cax = ax.inset_axes([1.022, .04, .012, .36])
    cb = fig.colorbar(sm, cax=cax, ticks=[-2, -1, 0, 1, 2]); cb.set_label("Scaled avg.\nexpression", fontsize=12); cb.ax.tick_params(labelsize=10)
    fig.tight_layout()
    for ext in ["png", "pdf"]: fig.savefig(N / f"marker_dotplot_panelD.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    rec = dict(status="PASS", cells=int(len(lab)), states=order, cells_per_state=lab.new_state.value_counts().reindex(order).to_dict(),
               umap_source_sha256=sha(R / "named_graph_PRIVATE.h5ad"), genes=len(genes), gene_order=genes,
               note="Labels/order/style only; cells, counts, UMAP coordinates and the 79-gene panel unchanged. reviewed10/ not modified.")
    (N / "RELABEL_RECEIPT.json").write_text(json.dumps(rec, indent=2)); return mp, rec
if __name__ == "__main__":
    mp, rec = run_pipeline(); print(mp[["old_cluster", "old_name", "new_cluster", "new_name"]].to_string(index=False)); print(rec["cells_per_state"])

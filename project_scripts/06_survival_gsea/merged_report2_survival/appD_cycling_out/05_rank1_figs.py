"""Appendix D (cycling cells excluded) figures for the best merged candidate (rank 1, two neighbouring Hotspot modules; Zejian 2026-10-11): (a) pooled Hotspot
local-correlation heatmap in dendrogram order with the module bar and the candidate boxed, (b) score_genes activity per tumor state, (c) UMAP of
the candidate score (own colour scale). Same styles as Fig 4.5 / 4.6. Env: moscot."""
import json, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
O = Path(__file__).resolve().parent; O5 = O.parents[2]; A = O5 / "metaNorm_20261011_all_cycout_bcr"; RV10 = O5.parents[1]; MERGED = O5.parents[5]
assert (MERGED / "shared_counts.h5ad").exists() and (RV10 / "relabel_20261008/STATE_MAP.csv").exists()
plt.rcParams.update({"font.family": "DejaVu Sans"}); RCOL = "#cc1a1a"
T = pd.read_csv(O / "pooled_top.csv", index_col=0); top = T.index[0]; parts = [int(x[1:]) for x in T.loc[top, "hotspot"].split("+")]
genes = json.load(open(O / "pooled_top_mouse.json"))[top]
lid = "META_all_FSU_top1000_n3000"; S = json.load(open(A / "hotspot_out" / f"{lid}_k30_t1000.json"))["sets"]["hs_top1000_min30"]
Z = np.load(A / "hotspot_out" / f"{lid}_k30_t1000_lcz.npy"); zg = (A / "hotspot_out" / f"{lid}_k30_t1000_lcz_genes.txt").read_text().split(); pos = {x: i for i, x in enumerate(zg)}
mem = {x: int(m) for m in S["order"] for x in S["modules"][str(m)]}; ordg = [x for m in S["order"] for x in sorted(S["modules"][str(m)])] + sorted(set(zg) - set(mem)); n = len(ordg)
fig, ax = plt.subplots(figsize=(10.5, 10)); ax.imshow(np.clip(Z[np.ix_([pos[x] for x in ordg], [pos[x] for x in ordg])], -8, 8), cmap="RdBu_r", vmin=-8, vmax=8, interpolation="nearest")
cm = plt.get_cmap("tab20"); hcol = {m: cm(i % 20) for i, m in enumerate(S["order"])}
for i, x in enumerate(ordg):
    m = mem.get(x); ax.add_patch(Rectangle((i - .5, -60), 1, 24, color=hcol[m] if m else (.88, .88, .88), clip_on=False))
    r = m in parts; ax.add_patch(Rectangle((i - .5, -32), 1, 24, color=RCOL if r else "white", clip_on=False)); ax.add_patch(Rectangle((-32, i - .5), 24, 1, color=RCOL if r else "white", clip_on=False))
ax.text(-40, -48, f"Hotspot modules ({len(S['order'])})", ha="right", va="center", fontsize=9); ax.text(-40, -20, "rank-1 candidate", ha="right", va="center", fontsize=9)
idx = [i for i, x in enumerate(ordg) if mem.get(x) in parts]; a0, a1 = min(idx), max(idx); ax.add_patch(Rectangle((a0 - .5, a0 - .5), a1 - a0 + 1, a1 - a0 + 1, fill=False, ec="black", lw=2.2))
ax.annotate(f"rank-1 merged candidate\n({len(genes)} genes)", xy=(a0 - 4, (a0 + a1) / 2), xytext=(max(a0 - 420, 20), min(a1 + 160, n - 60)), fontsize=11, weight="bold", bbox=dict(fc="white", ec="black", lw=.8), arrowprops=dict(arrowstyle="->", lw=1.3))
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(-33, n - .5); ax.set_ylim(n - .5, -61)
ax.set_title(f"Merged Npp53 + Npp53Rb1, cycling cells excluded: Hotspot local correlations of the top-1,000 fitness-up genes (dendrogram order)\n{len(S['order'])} Hotspot modules (min 30 genes) + {n - len(mem)} unassigned (grey); second bar and box = rank-1 candidate", fontsize=10, loc="left")
fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-8, 8), cmap="RdBu_r"), ax=ax, fraction=.035, pad=.01).set_label("local-correlation Z (clipped ±8)")
fig.savefig(O / "figC1_rank1_hotspot.png", dpi=130, bbox_inches="tight"); plt.close(fig)
SM = pd.read_csv(RV10 / "relabel_20261008/STATE_MAP.csv").sort_values("order"); SCOL = dict(zip(SM.new_state, SM.colour)); states = SM.new_state.tolist()
labs = pd.read_csv(RV10 / "relabel_20261008/cell_labels_PRIVATE.csv", usecols=["cell", "new_state"]).set_index("cell").new_state
sh = ad.read_h5ad(MERGED / "shared_counts.h5ad"); sh = sh[sh.obs_names.isin(labs.index)].copy(); sh.var_names = sh.var.gene_symbol.astype(str).values; sh.var_names_make_unique()
sc.pp.normalize_total(sh, target_sum=1e4); sc.pp.log1p(sh); sc.tl.score_genes(sh, [g for g in genes if g in sh.var_names], score_name="s", random_state=42)
v = sh.obs["s"].values; st = labs.reindex(sh.obs_names).values; mm = pd.Series(v).groupby(st).mean().reindex(states); z = (mm - mm.mean()) / mm.std(); z.to_csv(O / "rank1_state_z.csv")
gph = ad.read_h5ad(RV10 / "named_graph_PRIVATE.h5ad", backed="r"); um = pd.DataFrame(gph.obsm["X_umap"], index=gph.obs_names).reindex(sh.obs_names).values
fig = plt.figure(figsize=(15, 5.2)); gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], wspace=.25)
ax = fig.add_subplot(gs[0]); cols = ["#c44e52" if x > 0 else "#4a7fbf" for x in z.values]; ax.bar(range(len(states)), z.values, color=cols)
ax.axhline(0, color="k", lw=.8); ax.set_xticks(range(len(states))); ax.set_xticklabels(states, rotation=35, ha="right", fontsize=9)
for t, s in zip(ax.get_xticklabels(), states): t.set_color(SCOL[s])
for i, x in enumerate(z.values): ax.text(i, x + (.06 if x > 0 else -.06), f"{x:.1f}", ha="center", va="bottom" if x > 0 else "top", fontsize=8.5)
ax.set_ylabel("module score, z across states"); ax.set_title(f"b  Rank-1 merged candidate (no cycling cells) activity per tumor state (all {sh.n_obs:,} merged cells)", loc="left", fontsize=10.5)
for s_ in ["top", "right"]: ax.spines[s_].set_visible(False)
ax = fig.add_subplot(gs[1]); o = np.argsort(v); lo, hi = np.percentile(v, [2, 98])
p = ax.scatter(um[o, 0], um[o, 1], c=np.clip(v[o], lo, hi), s=.4, cmap="viridis_r", vmin=lo, vmax=hi, rasterized=True); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
fig.colorbar(p, ax=ax, fraction=.045, pad=.01).set_label("module score"); ax.set_title(f"c  Rank-1 merged candidate (no cycling cells) score ({len(genes)} genes)", loc="left", fontsize=10.5)
fig.savefig(O / "figC2_rank1_states_umap.png", dpi=150, bbox_inches="tight"); plt.close(fig); print(z.round(2).to_string())

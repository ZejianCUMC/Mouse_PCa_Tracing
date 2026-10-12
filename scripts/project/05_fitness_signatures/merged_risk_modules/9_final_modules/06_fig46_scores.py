"""Fig 4.6 for the final modules (Zejian 2026-10-11: "draw final module score project into UMAP for Npp53 and Npp53Rb1 in Fig 4.6").
Scores every final module on all merged cells exactly as the earlier Fig 4.6 did (normalize_total 1e4 + log1p; scanpy score_genes,
random_state 42), so the two risk-module scores are unchanged. a = mean score per tumor state (z across states); b = UMAPs, one per
final module. Writes final_module_state_z.csv and final_module_scores_PRIVATE.csv.gz. Env: moscot."""
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
O = Path(__file__).resolve().parent; P = O.parents[1]; RV10 = P.parents[2]; MERGED = P.parents[6]
assert (RV10 / "relabel_20261008/STATE_MAP.csv").exists() and (MERGED / "shared_counts.h5ad").exists()
plt.rcParams.update({"font.family": "DejaVu Sans"})
F = pd.read_csv(O / "final_modules.csv"); F["label"] = F.genotype + " " + F.final_module
SM = pd.read_csv(RV10 / "relabel_20261008/STATE_MAP.csv").sort_values("order"); SCOL = dict(zip(SM.new_state, SM.colour)); states = SM.new_state.tolist()
labs = pd.read_csv(RV10 / "relabel_20261008/cell_labels_PRIVATE.csv", usecols=["cell", "new_state"]).set_index("cell").new_state
import os; AUC = os.environ.get("SCORE", "scanpy") == "aucell"   # Zejian 2026-10-11: AUCell scores (11_aucell.R); SCORE=scanpy gives the earlier score_genes version
CACHE = O / ("final_module_aucell_PRIVATE.csv.gz" if AUC else "final_module_scores_PRIVATE.csv.gz")
if AUC:
    SC = pd.read_csv(CACHE, index_col=0); sco = {c: SC[c].values.astype(float) for c in F.label}; cells = SC.index; st = labs.reindex(cells).values
    Z = pd.DataFrame({c: (lambda mm: (mm - mm.mean()) / mm.std())(pd.Series(v).groupby(st).mean().reindex(states)) for c, v in sco.items()}).T
    Z.to_csv(O / "final_module_state_z_aucell.csv")
elif CACHE.exists() and set(pd.read_csv(CACHE, nrows=1, index_col=0).columns) == set(F.label):   # plotting-only rerun: reuse cached scores
    SC = pd.read_csv(CACHE, index_col=0); sco = {c: SC[c].values.astype(float) for c in F.label}; cells = SC.index; st = labs.reindex(cells).values
    Z = pd.DataFrame({c: (lambda mm: (mm - mm.mean()) / mm.std())(pd.Series(v).groupby(st).mean().reindex(states)) for c, v in sco.items()}).T
    Z.to_csv(O / "final_module_state_z.csv")
else:
    sh = ad.read_h5ad(MERGED / "shared_counts.h5ad"); sh = sh[sh.obs_names.isin(labs.index)].copy(); sh.var_names = sh.var.gene_symbol.astype(str).values; sh.var_names_make_unique()
    sc.pp.normalize_total(sh, target_sum=1e4); sc.pp.log1p(sh); st = labs.reindex(sh.obs_names).values
    Zs, sco = {}, {}
    for r in F.itertuples():
        gl = [x for x in r.genes.split() if x in sh.var_names]; sc.tl.score_genes(sh, gl, score_name="tmp", random_state=42)
        v = sh.obs["tmp"].values.copy(); sco[r.label] = v; mm = pd.Series(v).groupby(st).mean().reindex(states); Zs[r.label] = (mm - mm.mean()) / mm.std()
    Z = pd.DataFrame(Zs).T; Z.to_csv(O / "final_module_state_z.csv")   # full precision (avoid double rounding in labels)
    pd.DataFrame(sco, index=sh.obs_names).astype("float32").to_csv(O / "final_module_scores_PRIVATE.csv.gz", compression="gzip")
    cells = sh.obs_names
gph = ad.read_h5ad(RV10 / "named_graph_PRIVATE.h5ad", backed="r"); um = pd.DataFrame(gph.obsm["X_umap"], index=gph.obs_names).reindex(cells).values
ng = F.groupby("genotype").size(); ncol = int(ng.max())
fig = plt.figure(figsize=(4.1 * ncol, 4.6 + 3.9 * len(ng))); gs = fig.add_gridspec(1 + len(ng), ncol, height_ratios=[1.25] + [1] * len(ng), hspace=.62, wspace=.28)
ax = fig.add_subplot(gs[0, :]); im = ax.imshow(Z.values, cmap="RdBu_r", vmin=-2.2, vmax=2.2, aspect="auto")
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]): ax.text(j, i, f"{Z.values[i, j]:.1f}", ha="center", va="center", fontsize=8.5, color="white" if abs(Z.values[i, j]) > 1.4 else "black")
ax.set_yticks(range(len(Z))); ax.set_yticklabels(Z.index, fontsize=9.5); ax.set_xticks(range(len(states))); ax.set_xticklabels(states, rotation=30, ha="right", fontsize=9)
for t, s in zip(ax.get_xticklabels(), states): t.set_color(SCOL[s])
for t in ax.get_yticklabels():
    if "risk" in t.get_text(): t.set_fontweight("bold")
fig.colorbar(im, ax=ax, fraction=.02, pad=.01).set_label(("AUCell" if AUC else "module score") + ", z across states")
ax.set_title(f"a  Final-module activity per tumor state (all {len(cells):,} merged cells)", loc="left", fontsize=11)
for gi, g in enumerate(ng.index):
    rows = list(F[F.genotype == g].itertuples())
    for k, r in enumerate(rows):
        ax = fig.add_subplot(gs[1 + gi, k]); v = sco[r.label]; o = np.argsort(v); lo, hi = np.percentile(v, [2, 98])   # own colour scale per panel (Zejian 2026-10-11, after trying one scale per row)
        p = ax.scatter(um[o, 0], um[o, 1], c=np.clip(v[o], lo, hi), s=.35, cmap="viridis_r", vmin=lo, vmax=hi, rasterized=True); ax.set_xticks([]); ax.set_yticks([])
        if k == 0: ax.set_ylabel("UMAP 2")
        ax.set_xlabel("UMAP 1")
        fig.colorbar(p, ax=ax, fraction=.045, pad=.01).set_label("AUCell" if AUC else "module score", fontsize=8)
        ax.set_title(f"{'b' if gi == 0 else 'c'}{k + 1}  {g} {r.final_module} ({r.n_genes} genes)", loc="left", fontsize=10, fontweight="bold" if "risk" in r.final_module else "normal")
fig.savefig(O / ("fig4_6_final_modules_umap.png"), dpi=150, bbox_inches="tight"); plt.close(fig); print(Z.round(2).to_string())

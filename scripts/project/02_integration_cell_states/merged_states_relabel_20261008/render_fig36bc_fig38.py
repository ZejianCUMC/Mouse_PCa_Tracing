"""Companion figures for the 2026-10-08 relabel (conda env: moscot):
1) umap_treatment_run.png — treatment and run UMAPs in the same style as umap_states_numbered.png (kept from the old Fig 3.6);
2) aucell_heatmap_relabel.png — Fig 3.8 re-rendered from the SAVED seed-42 outputs (selected pathways, z-scores, annotations),
   columns re-indexed/re-ordered to the new cluster IDs. Row z-scores are computed over all 10 states, so re-ordering columns
   changes no value; no statistic is recomputed."""
from pathlib import Path
import json, textwrap, numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
N = Path(__file__).resolve().parent; R = N.parent; P = R / "pathway_statistics"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial"], "pdf.fonttype": 42})
smap = pd.read_csv(N / "STATE_MAP.csv"); new2old = dict(zip(smap.new_cluster, smap.old_cluster)); new = smap.new_cluster.tolist()
# ── 1. treatment / run UMAPs ─────────────────────────────────────────────────────────────────────
a = ad.read_h5ad(R / "named_graph_PRIVATE.h5ad"); xy = a.obsm["X_umap"]; rng = np.random.default_rng(42).permutation(len(xy))
treat = {"Intact": "#4c78a8", "Castration": "#f58518", "Castration+ENZA": "#54a24b"}
runs = {r: matplotlib.colors.to_hex(plt.get_cmap("tab20")(i)) for i, r in enumerate(sorted(a.obs.run_id.astype(str).unique()))}
fig, axes = plt.subplots(1, 2, figsize=(15, 6.2))
for ax, col, pal, title, ncol in [(axes[0], "treatment", treat, "Treatment", 1), (axes[1], "run_id", runs, "Run", 2)]:
    v = a.obs[col].astype(str).values; ax.scatter(xy[rng, 0], xy[rng, 1], c=[pal[x] for x in v[rng]], s=0.8, linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); [ax.spines[k].set_visible(False) for k in ax.spines]
    x0, x1, y0, y1 = xy[:, 0].min(), xy[:, 0].max(), xy[:, 1].min(), xy[:, 1].max(); ox, oy = x0 - .05 * (x1 - x0), y0 - .05 * (y1 - y0)
    ax.plot([ox, ox + .16 * (x1 - x0)], [oy, oy], color="k", lw=1.4, clip_on=False); ax.plot([ox, ox], [oy, oy + .16 * (y1 - y0)], color="k", lw=1.4, clip_on=False)
    ax.text(ox + .08 * (x1 - x0), oy - .035 * (y1 - y0), "UMAP1", ha="center", va="top", fontsize=12)
    ax.text(ox - .02 * (x1 - x0), oy + .08 * (y1 - y0), "UMAP2", ha="right", va="center", rotation=90, fontsize=12)
    h = [plt.Line2D([], [], marker="o", ls="", color=c, markersize=8, label=k) for k, c in pal.items()]
    ax.legend(handles=h, title=title, loc="center left", bbox_to_anchor=(1.0, .5), frameon=False, fontsize=12, title_fontsize=13, ncol=ncol, columnspacing=1.0)
fig.tight_layout(); fig.savefig(N / "umap_treatment_run.png", dpi=200, bbox_inches="tight"); fig.savefig(N / "umap_treatment_run.pdf", dpi=200, bbox_inches="tight"); plt.close(fig)
# ── 2. Fig 3.8 re-indexed ─────────────────────────────────────────────────────────────────────────
z = pd.read_csv(P / "selected_heatmap_zscores.csv", index_col=0); chosen = pd.read_csv(P / "selected_pathways.csv"); ann = pd.read_csv(P / "heatmap_annotations.csv").fillna("")
st = pd.read_csv(P / "statistics_seed42.csv").set_index(["cluster", "signature"])
selected = z.index.tolist(); assert selected == chosen.signature.tolist()
zn = z[[new2old[c] for c in new]]; zn.columns = new
A = ann.pivot(index="signature", columns="cluster", values="annotation").reindex(index=selected)[[new2old[c] for c in new]]; A.columns = new
ns = [int(st.loc[new2old[c]].eligible_samples.iloc[0]) for c in new]
readable = {'HALLMARK_TNFA_SIGNALING_VIA_NFKB':'TNFα signaling via NF-κB','HALLMARK_E2F_TARGETS':'E2F targets','HALLMARK_MYC_TARGETS_V1':'MYC targets V1','Dehm_luminal':'Dehm luminal','Hieronymus_AR_signaling':'Hieronymus AR signaling','Tang_CRPC_SCL':'Tang CRPC-SCL','MetC':'MetC','CCP_31':'CCP-31 (cell-cycle progression)','REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY':'Type I hemidesmosome assembly','REACTOME_ECM_PROTEOGLYCANS':'ECM proteoglycans','REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS':'TGF-β receptor signaling → SMADs','REACTOME_FOLDING_OF_ACTIN_BY_CCT_TRIC':'Actin folding by CCT/TRiC','KEGG_MAPK_SIGNALING_PATHWAY':'MAPK signaling','KEGG_CITRATE_CYCLE_TCA_CYCLE':'Citrate cycle (TCA cycle)'}
short = [readable.get(n, n.replace('HALLMARK_','').replace('REACTOME_','').replace('KEGG_','').replace('_',' ').lower().capitalize()) for n in selected]
fig, ax = plt.subplots(figsize=(12, max(6, len(selected) * .29 + 1.8)))
im = ax.imshow(zn.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_xticks(range(10), new, fontsize=11); ax.set_yticks(range(len(selected)), [textwrap.fill(v, 48) for v in short], fontsize=10)
ax.tick_params(length=0); ax.set_xticks(np.arange(-.5, 10, 1), minor=True); ax.set_yticks(np.arange(-.5, len(selected), 1), minor=True); ax.grid(which="minor", color="white", linewidth=.6); ax.tick_params(which="minor", length=0)
for i in range(len(selected)):
    for j in range(10):
        m = A.iloc[i, j]
        if m: ax.text(j, i, m, ha="center", va="center", fontsize=8 if m == "NT" else 9, color="white" if abs(zn.iloc[i, j]) > 1.15 else "#222222")
off = 0
for coll in ["Hallmark", "ProstateCancerRelated", "Reactome", "KEGG_Legacy"]:
    cnt = int(chosen.collection.eq(coll).sum())
    if cnt:
        if off: ax.axhline(off - .5, color="black", linewidth=1.3)
        ax.text(10.05, off + (cnt - 1) / 2, {"ProstateCancerRelated": "Prostate cancer", "KEGG_Legacy": "KEGG legacy"}.get(coll, coll), va="center", ha="center", rotation=90, fontsize=10, weight="bold"); off += cnt
ax.set_title("Representative AUCell programs · seed 42", fontsize=13, pad=12)
fig.subplots_adjust(left=.40, right=.86, top=.94, bottom=.16)
cax = fig.add_axes([.93, .70, .014, .19]); cb = fig.colorbar(im, cax=cax); cb.set_label("Row z-score\nsample-balanced mean AUCell", fontsize=9)
fig.text(.40, .025, "* q<0.05   ** q<0.01   *** q<0.001 (BH FDR)\n↑ / ↓: higher / lower than other states within sample.\nNT: insufficient samples; color is descriptive.\nPaired samples per c1–c10: " + ", ".join(map(str, ns)), fontsize=9, va="bottom")
fig.savefig(N / "aucell_heatmap_relabel.png", dpi=180, bbox_inches="tight"); plt.close(fig)
(N / "FIG36_38_RECEIPT.json").write_text(json.dumps(dict(status="PASS", fig38_columns_new=new, fig38_columns_from_old=[new2old[c] for c in new], eligible_samples_new_order=ns,
     note="Fig 3.8 values/stars identical to saved seed-42 outputs; columns re-indexed only."), indent=2))
print("PASS", ns)

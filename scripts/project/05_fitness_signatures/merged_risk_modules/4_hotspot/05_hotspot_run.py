"""opt5 step 05 — Hotspot for one list (usage: 05_hotspot_run.py <arm dir> <list_id> <k>). Unchanged maths from opt2-opt4 step 05
(reviewed 2026-10-09, no defect found): danb model on raw counts, umi = all-gene UMI, unweighted EXACT KNN (approx_neighbors=False,
deterministic) on the merged scVI latent; autocorrelation -> top 500 by Z among FDR < 0.05 -> pairwise local-correlation Z ->
create_modules(core_only, fdr 0.05). Module sets saved:
  hs_top{300,500}_min{10,15,20,30} : Hotspot modules (dendrogram order kept)
  ward{3,5}_top500                 : Yang 2022-style clusters (ward linkage on max(Z) - Z of the top-500 local-correlation Z)
  merge5_top500_min20              : adjacent (dendrogram) modules with the highest mean between-module Z merged until 5 remain
Inputs <arm>/fsu_lists.json, <arm>/hotspot_inputs/<cellset>_*; output <arm>/hotspot_out/<list_id>_k<k>.json (+ _lcz.npy, _lcz_genes.txt).
Conda env: hotspot (local and c2b2 HPC envs have identical versions: hotspotsc 1.1.3, numpy 1.21.6, scikit-learn 1.0.2)."""
import sys, json, time, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, hotspot
from pathlib import Path
from scipy.cluster.hierarchy import leaves_list, linkage, fcluster
from scipy.spatial.distance import squareform
A = Path(sys.argv[1]).resolve(); lid, k = sys.argv[2], int(sys.argv[3]); D = A / "hotspot_inputs"; OUT = A / "hotspot_out"; OUT.mkdir(exist_ok=True); t0 = time.time()
import os; TOPN = int(os.environ.get("TOPN", 500)); SUF = "" if TOPN == 500 else f"_t{TOPN}"   # round 2+: top-N autocorrelated genes (500 default, as rounds 0-1)
L = json.load(open(A / "fsu_lists.json"))[lid]; v = L["cells"]
b = np.load(D / f"{v}_bundle_PRIVATE.npz"); X = sp.csr_matrix((b["data"], b["indices"], b["indptr"]), shape=tuple(b["shape"])).tocsc(); allg = b["genes"].astype(str)
assert X.shape[1] == len(allg) and set(L["genes"]) <= set(allg), "input/gene mismatch"
obs = pd.read_csv(D / f"{v}_obs_PRIVATE.csv").set_index("cell"); lat = np.load(D / f"{v}_latent_PRIVATE.npy"); assert len(obs) == X.shape[0] == lat.shape[0]
gi = [i for i, g in enumerate(allg) if g in set(L["genes"])]; Xs = X[:, gi]; gs = allg[gi]; nz = np.asarray((Xs > 0).sum(0)).ravel() > 0
Ad = ad.AnnData(X=Xs[:, nz], obs=obs, var=pd.DataFrame(index=gs[nz])); Ad.layers["counts"] = Ad.X; Ad.obsm["X_scVI"] = lat
hs = hotspot.Hotspot(Ad, layer_key="counts", model="danb", latent_obsm_key="X_scVI", umi_counts_obs_key="total_counts")
hs.create_knn_graph(weighted_graph=False, n_neighbors=k, approx_neighbors=False); res = hs.compute_autocorrelations(jobs=5)
top500 = res.query("FDR<0.05").sort_values("Z", ascending=False).index[:TOPN].tolist(); lcz = hs.compute_local_correlations(genes=top500, jobs=5)   # name kept; holds the top TOPN genes
out = dict(list_id=lid, k=k, top=TOPN, n_genes=int(nz.sum()), n_autocorr_fdr05=int((res.FDR < .05).sum()), n_top=len(top500), sets={})
def ordered(mods):
    order = [hs.local_correlation_z.index[i] for i in leaves_list(hs.linkage)]; mo = [int(m) for m in pd.unique(mods.reindex(order).values) if m > 0]
    return mo, {int(m): sorted(mods[mods == m].index.tolist()) for m in mo}
for top in sorted({300, 500, TOPN}):
    sub = top500[:top]; hs.local_correlation_z = lcz.loc[sub, sub]
    for mn in [10, 15, 20, 30]:
        mods = hs.create_modules(min_gene_threshold=mn, core_only=True, fdr_threshold=0.05); mo, mm = ordered(mods)
        out["sets"][f"hs_top{top}_min{mn}"] = dict(order=mo, modules=mm)
        if top == TOPN and mn == 20 and len(mo) > 0:
            groups = [[m] for m in mo]; Z = lcz.loc[top500, top500]
            mz = lambda a, b: Z.loc[[g for m in a for g in mm[m]], [g for m in b for g in mm[m]]].values.mean()
            while len(groups) > 5:
                s = [mz(groups[i], groups[i + 1]) for i in range(len(groups) - 1)]; i = int(np.argmax(s)); groups[i:i + 2] = [groups[i] + groups[i + 1]]
            out["sets"][f"merge5_top{TOPN}_min20"] = dict(order=list(range(1, len(groups) + 1)), groups=groups, modules={i + 1: sorted(g for m in grp for g in mm[m]) for i, grp in enumerate(groups)})
if len(top500) >= 6:
    dd = lcz.loc[top500, top500].values.copy(); np.fill_diagonal(dd, 0); cond = squareform((dd + dd.T) / 2, checks=False) * -1; cond += -cond.min()
    Zw = linkage(cond, method="ward")
    for nc in [3, 5]:
        cl = fcluster(Zw, t=nc, criterion="maxclust"); order = [top500[i] for i in leaves_list(Zw)]; mo = [int(m) for m in pd.unique(pd.Series(cl, index=top500).reindex(order).values)]
        out["sets"][f"ward{nc}_top{TOPN}"] = dict(order=mo, modules={int(m): sorted(np.array(top500)[cl == m].tolist()) for m in mo})
np.save(OUT / f"{lid}_k{k}{SUF}_lcz.npy", lcz.loc[top500, top500].values.astype(np.float32)); (OUT / f"{lid}_k{k}{SUF}_lcz_genes.txt").write_text("\n".join(top500) + "\n")
out["seconds"] = round(time.time() - t0, 1); (OUT / f"{lid}_k{k}{SUF}.json").write_text(json.dumps(out)); print(A.name, lid, k, TOPN, out["seconds"], "s", flush=True)

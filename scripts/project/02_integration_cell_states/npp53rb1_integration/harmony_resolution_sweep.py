"""Compare broad-lineage resolutions without refitting Harmony or UMAP.

Selection is a subsequent marker-review decision, not an automatic purity maximum.
Existing fastMNN results and the source AnnData are never modified.
"""
import json
import os
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

run = Path(os.environ['RB1_RUN'])
out = run / 'harmony_major_types'
out.mkdir(exist_ok=True)
if (out / 'sweep_complete.json').exists():
    raise RuntimeError('Completed sweep exists; refusing overwrite')
source = sc.read_h5ad(run / 'major_types_full.h5ad', backed='r')
assert source.n_obs == 56686 and source.obs.sample_id.nunique() == 22
assert not source.obs.sample_id.isin(['JZ257', 'JZ275', 'JZ278', 'JZ281']).any()
markers = ['Krt5','Krt14','Trp63','Krt17','Krt8','Krt18','Krt19','Epcam',
           'Cdh1','Vim','Zeb1','Fn1','Col1a1','Col1a2','Col3a1','Dcn','Lum',
           'Pdgfra','Lyz2','Tyrobp','Csf1r','Adgre1','Ptprc','Cd74',
           'Cd3d','Cd3e','Trbc2','Nkg7','Klrd1','Cd79a','Ms4a1',
           'Pecam1','Cdh5','Emcn','Vwf','Rgs5','Pdgfrb','Acta2','Tagln',
           'Syp','Chga','Chgb','Ascl1','Insm1','Mki67','Top2a','GFP','mCherry']
markers = [g for g in markers if g in source.var_names]
a = source[:, markers].to_memory()
a.obsm['X_umap'] = source.obsm['X_umap_Harmony_theta1'].copy()
source.file.close()
assert np.isfinite(a.obsm['X_Harmony_theta1']).all()
sc.pp.neighbors(a, use_rep='X_Harmony_theta1', n_neighbors=15, random_state=19)
expr = a.to_df()
resolutions = [.1, .2, .3, .5, .8]
seeds = [19, 29, 39]
rows = []
fig, axes = plt.subplots(2, 3, figsize=(24, 15))
for res, ax in zip(resolutions, axes.flat):
    key = 'leiden_r' + str(res).replace('.', 'p')
    labels = []
    for seed in seeds:
        sc.tl.leiden(a, resolution=res, key_added=key, flavor='igraph',
                     n_iterations=2, directed=False, random_state=seed)
        labels.append(a.obs[key].astype(str).copy())
    a.obs[key] = pd.Categorical(labels[0])
    group = a.obs[key]
    sizes = group.value_counts()
    cross = pd.crosstab(group, a.obs.expression_anchor)
    known = cross.drop(columns='Unresolved', errors='ignore')
    cross.to_csv(out / f'{key}_anchors.csv')
    expr.groupby(group, observed=True).mean().to_csv(out / f'{key}_marker_means.csv')
    (expr > 0).groupby(group, observed=True).mean().to_csv(out / f'{key}_marker_fractions.csv')
    pd.crosstab(group, a.obs.sample_id).to_csv(out / f'{key}_samples.csv')
    rows.append(dict(resolution=res, clusters=len(sizes), smallest_cluster=int(sizes.min()),
                     clusters_below_100=int((sizes < 100).sum()),
                     anchor_purity_descriptive=float(known.max(axis=1).sum()/known.values.sum()),
                     ARI_19_29=adjusted_rand_score(labels[0],labels[1]),
                     ARI_19_39=adjusted_rand_score(labels[0],labels[2]),
                     ARI_29_39=adjusted_rand_score(labels[1],labels[2])))
    xy = a.obsm['X_umap']
    for i, cluster in enumerate(sorted(sizes.index, key=int)):
        mask = np.asarray(group == cluster)
        points = xy[mask]
        ax.scatter(*points.T, s=3, color=plt.get_cmap('tab20')(i % 20), linewidths=0)
        median = np.median(points, axis=0)
        pos = points[np.argmin(((points-median)**2).sum(axis=1))]
        ax.text(*pos, cluster, fontsize=17, weight='bold', ha='center',
                path_effects=[pe.withStroke(linewidth=3, foreground='white')])
    ax.set_title(f'Leiden {res}: {len(sizes)} clusters', fontsize=20)
    ax.set(xticks=[], yticks=[], xlabel='UMAP_1', ylabel='UMAP_2', aspect='equal')
    ax.spines[['top','right','bottom','left']].set_visible(False)
    pd.DataFrame(rows).to_csv(out / 'resolution_diagnostics.csv', index=False)
    print(rows[-1], flush=True)
axes.flat[-1].axis('off')
fig.suptitle('Fixed Harmony θ1 UMAP • resolution review • 56,686 cells', fontsize=23)
fig.tight_layout()
fig.savefig(out / 'resolution_review.png', dpi=160)
plt.close(fig)
a.obs.to_csv(out / 'resolution_cell_annotations.csv')
a.write_h5ad(out / 'harmony_resolution_review.h5ad')
(out / 'sweep_complete.json').write_text(json.dumps({
    'source':str(run/'major_types_full.h5ad'), 'cells':a.n_obs,
    'representation':'X_Harmony_theta1','batch_key':'library',
    'UMAP':'unchanged X_umap_Harmony_theta1','resolutions':resolutions,'seeds':seeds,
    'neighbors':15,'leiden_iterations':2,'selection':'pending marker review',
    'caution':'Anchor purity is circular and rises with granularity; not an optimization objective. Reporter expression is not used for clustering or resolution selection.'
}, indent=2))
print('Sweep complete; manual marker review required.', flush=True)

import scanpy as sc, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
a=sc.read_h5ad(f"{OUT}/npp53_scvi.h5ad")
md=pd.read_csv(f"{OUT}/npp53_putative_celltype.csv").set_index("cell")
a.obs["putative"]=md["putative"].reindex(a.obs_names).astype(str).values
a.obs["cluster"]=md["cluster"].reindex(a.obs_names).astype(str).values
a=a[a.obs["cluster"]!="5"].copy()                       # drop low-quality cluster 5
sc.pp.neighbors(a, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(a)
ordr=["Epithelial_Luminal","Epithelial_Basal","Mesenchymal_EMT","Fibroblast","Endothelial","Macrophage_Myeloid"]
present=[c for c in ordr if c in set(a.obs["putative"])]
a.obs["putative"]=pd.Categorical(a.obs["putative"], categories=present)
pal={"Epithelial_Luminal":"#1f78b4","Epithelial_Basal":"#a6cee3","Mesenchymal_EMT":"#e31a1c",
     "Fibroblast":"#ff7f00","Endothelial":"#33a02c","Macrophage_Myeloid":"#6a3d9a"}
fig=sc.pl.umap(a, color="putative", palette=[pal[c] for c in present], show=False, return_fig=True,
               title="scVI — putative cell type (low-quality removed)")
fig.savefig(f"{OUT}/UMAP_scVI_putative_clean.png", dpi=150, bbox_inches="tight"); plt.close(fig)
a.write(f"{OUT}/npp53_scvi_clean.h5ad")
print("DONE scVI rebuild; cells:", a.n_obs)

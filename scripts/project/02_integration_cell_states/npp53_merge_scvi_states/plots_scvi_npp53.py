import scanpy as sc, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
a=sc.read_h5ad(f"{OUT}/npp53_scvi.h5ad")
put=pd.read_csv(f"{OUT}/npp53_putative_celltype.csv").set_index("cell")["putative"]
a.obs["putative"]=put.reindex(a.obs_names).astype(str).values
ordr=["Epithelial_Luminal","Epithelial_Basal","Mesenchymal_EMT","Fibroblast","Endothelial","Macrophage_Myeloid","Bcell"]
present=[c for c in ordr if c in set(a.obs["putative"])]
a.obs["putative"]=pd.Categorical(a.obs["putative"], categories=present)
pal={"Epithelial_Luminal":"#1f78b4","Epithelial_Basal":"#a6cee3","Mesenchymal_EMT":"#e31a1c",
     "Fibroblast":"#ff7f00","Endothelial":"#33a02c","Macrophage_Myeloid":"#6a3d9a","Bcell":"#b15928"}
fig=sc.pl.umap(a, color="putative", palette=[pal[c] for c in present], show=False, return_fig=True,
               title="scVI — putative cell type")
fig.savefig(f"{OUT}/UMAP_scVI_putative.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("DONE scVI putative plot")

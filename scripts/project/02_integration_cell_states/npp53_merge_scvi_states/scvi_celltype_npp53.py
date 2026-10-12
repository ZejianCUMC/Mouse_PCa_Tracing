import scanpy as sc, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
a=sc.read_h5ad(f"{OUT}/npp53_scvi_clean.h5ad")
md=pd.read_csv(f"{OUT}/npp53_putative_celltype_clean.csv").set_index("cell")
# rebuild celltype from cluster mapping (cl6->Partial-EMT, cl13+8->Mesenchymal) via the labeled RDS csv if present
lab=pd.read_csv(f"{OUT}/npp53_reporter_pos_cells.csv") if False else None
# simpler: read celltype straight from labeled object export
ctcsv=f"{OUT}/npp53_celltype_cells.csv"
import os
cc=pd.read_csv(ctcsv).set_index("cell") if os.path.exists(ctcsv) else None
if cc is not None:
    a.obs["celltype"]=cc["celltype"].reindex(a.obs_names).astype(str).values
else:
    a.obs["celltype"]=md["putative"].reindex(a.obs_names).astype(str).values
ordr=["Epithelial_Luminal","Epithelial_Basal","Partial-EMT","Mesenchymal","Endothelial","Macrophage_Myeloid"]
present=[c for c in ordr if c in set(a.obs["celltype"])]
a.obs["celltype"]=pd.Categorical(a.obs["celltype"], categories=present)
pal={"Epithelial_Luminal":"#1f78b4","Epithelial_Basal":"#a6cee3","Partial-EMT":"#fb9a99",
     "Mesenchymal":"#e31a1c","Endothelial":"#33a02c","Macrophage_Myeloid":"#6a3d9a"}
fig=sc.pl.umap(a, color="celltype", palette=[pal[c] for c in present], show=False, return_fig=True, title="scVI")
fig.savefig(f"{OUT}/UMAP_scVI_celltype.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("DONE scVI celltype")

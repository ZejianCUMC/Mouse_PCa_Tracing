import scanpy as sc, pandas as pd
OUT="/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
a=sc.read_h5ad(f"{OUT}/npp53_scvi_clean.h5ad")
key="X_umap" if "X_umap" in a.obsm else [k for k in a.obsm if "umap" in k.lower()][0]
df=pd.DataFrame(a.obsm[key][:,:2], index=a.obs_names, columns=["scVI_1","scVI_2"])
df.index.name="cell"; df.to_csv(f"{OUT}/scvi_umap_clean.csv")
print("exported scVI umap:", df.shape, "key=", key)

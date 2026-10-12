import scanpy as sc, pandas as pd
OUT="/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
a=sc.read_h5ad(f"{OUT}/npp53_scvi_clean.h5ad")
key="X_scVI" if "X_scVI" in a.obsm else [k for k in a.obsm if "scvi" in k.lower()][0]
df=pd.DataFrame(a.obsm[key], index=a.obs_names, columns=[f"scvi_{i+1}" for i in range(a.obsm[key].shape[1])])
df.index.name="cell"; df.to_csv(f"{OUT}/scvi_latent_clean.csv")
print("scVI latent:", df.shape, "key=", key)

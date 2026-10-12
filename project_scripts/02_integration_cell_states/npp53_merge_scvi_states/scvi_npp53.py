#!/usr/bin/env python
# scVI batch correction on the merged Npp53 cohort (from merge_npp53.R exports). GPU.
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp
import scvi, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
# build anndata from mtx (genes x cells) -> cells x genes
M=sio.mmread(f"{OUT}/npp53_counts.mtx").tocsr()            # genes x cells
genes=[l.strip() for l in open(f"{OUT}/npp53_genes.txt")]
bcs=[l.strip() for l in open(f"{OUT}/npp53_barcodes.txt")]
meta=pd.read_csv(f"{OUT}/npp53_meta.csv").set_index("cell").loc[bcs]
adata=ad.AnnData(X=sp.csr_matrix(M.T), obs=meta.copy(), var=pd.DataFrame(index=genes))
adata.layers["counts"]=adata.X.copy()
print("adata:", adata.shape, "batches:", adata.obs['batch'].value_counts().to_dict())
# HVG on raw counts (seurat_v3), batch-aware
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                            batch_key="batch", layer="counts", subset=True)
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model=scvi.model.SCVI(adata, n_latent=30, n_layers=2)
model.train(max_epochs=200, early_stopping=True)
adata.obsm["X_scVI"]=model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(adata)
# save
model.save(f"{OUT}/scvi_model", overwrite=True)
adata.write(f"{OUT}/npp53_scvi.h5ad")
# QC umap by batch + treatment
for col in ["batch","treatment"]:
    fig=sc.pl.umap(adata, color=col, show=False, return_fig=True, title=f"scVI · {col}")
    fig.savefig(f"{OUT}/UMAP_scVI_by{col}.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("DONE scVI ->", OUT)

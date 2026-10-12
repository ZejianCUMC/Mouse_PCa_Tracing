suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_annotated.RDS")); DefaultAssay(so) <- "RNA"
so  <- so[, so$seurat_clusters != "5"]                 # drop low-quality cluster 5
so$putative <- droplevels(so$putative)
cat("kept", ncol(so), "cells; types:", paste(levels(so$putative), collapse=", "), "\n")
# rebuild UMAP from each method's latent (subset embeddings)
so <- RunUMAP(so, reduction="pca",     dims=1:30, reduction.name="umap.uncorr",  verbose=FALSE)
so <- RunUMAP(so, reduction="harmony", dims=1:30, reduction.name="umap.harmony", verbose=FALSE)
so <- RunUMAP(so, reduction="mnn",     dims=1:30, reduction.name="umap.mnn",     verbose=FALSE)
saveRDS(so, file.path(OUT, "npp53_clean.RDS"))
ord <- intersect(c("Epithelial_Luminal","Epithelial_Basal","Mesenchymal_EMT","Fibroblast",
                   "Endothelial","Macrophage_Myeloid"), levels(so$putative))
so$putative <- factor(so$putative, levels=ord)
pal <- c(Epithelial_Luminal="#1f78b4", Epithelial_Basal="#a6cee3", Mesenchymal_EMT="#e31a1c",
         Fibroblast="#ff7f00", Endothelial="#33a02c", Macrophage_Myeloid="#6a3d9a")
for (r in c("umap.uncorr","umap.harmony","umap.mnn")) {
  p <- DimPlot(so, reduction=r, group.by="putative", cols=pal, raster=FALSE) +
       ggtitle(paste0(r, " — putative cell type (low-quality removed)"))
  ggsave(file.path(OUT, paste0("UMAP_", r, "_putative_clean.png")), p, width=9, height=6.5, dpi=150)
}
feats <- c("Epcam","Krt8","Krt18","Pbsn","Ar","Krt5","Krt14","Trp63","Vim","Col3a1","Zeb1","Fn1",
           "Dcn","Lum","Pdgfra","Pecam1","Cldn5","Cdh5","Lyz2","C1qa","Cd68","Csf1r")
feats <- feats[feats %in% rownames(so)]
p <- DotPlot(so, features=feats, group.by="putative") + RotatedAxis() +
     theme(axis.text.x=element_text(size=8)) + ggtitle("Marker genes by putative cell type (clean)")
ggsave(file.path(OUT, "DotPlot_putative_markers_clean.png"), p, width=11, height=5, dpi=150)
cat("=== composition (clean) ===\n"); print(table(so$putative))
write.csv(data.frame(cell=colnames(so), putative=so$putative, batch=so$batch,
                     sample_id=so$sample_id, treatment=so$treatment),
          file.path(OUT,"npp53_putative_celltype_clean.csv"), row.names=FALSE)
cat("DONE rebuild\n")

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_annotated.RDS")); DefaultAssay(so) <- "RNA"
ord <- c("Epithelial_Luminal","Epithelial_Basal","Mesenchymal_EMT","Fibroblast",
         "Endothelial","Macrophage_Myeloid","Bcell")
ord <- intersect(ord, levels(so$putative)); so$putative <- factor(so$putative, levels=ord)
pal <- c(Epithelial_Luminal="#1f78b4", Epithelial_Basal="#a6cee3", Mesenchymal_EMT="#e31a1c",
         Fibroblast="#ff7f00", Endothelial="#33a02c", Macrophage_Myeloid="#6a3d9a", Bcell="#b15928")
for (r in c("umap.uncorr","umap.harmony","umap.mnn")) {
  p <- DimPlot(so, reduction=r, group.by="putative", cols=pal, raster=FALSE) +
       ggtitle(paste0(r, " — putative cell type"))
  ggsave(file.path(OUT, paste0("UMAP_", r, "_putative.png")), p, width=9, height=6.5, dpi=150)
}
feats <- c("Epcam","Krt8","Krt18","Pbsn","Ar", "Krt5","Krt14","Trp63",
           "Vim","Col3a1","Zeb1","Fn1", "Dcn","Lum","Pdgfra",
           "Pecam1","Cldn5","Cdh5", "Lyz2","C1qa","Cd68","Csf1r",
           "Cd79a","Ms4a1","Igkc","Jchain")
feats <- feats[feats %in% rownames(so)]
p <- DotPlot(so, features=feats, group.by="putative") + RotatedAxis() +
     theme(axis.text.x=element_text(size=8)) + ggtitle("Marker genes by putative cell type")
ggsave(file.path(OUT, "DotPlot_putative_markers.png"), p, width=12, height=5.5, dpi=150)
cat("DONE plots ->", OUT, "\n")

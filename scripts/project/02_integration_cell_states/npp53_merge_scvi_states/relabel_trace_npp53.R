suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_clean.RDS")); DefaultAssay(so) <- "RNA"

# ---- relabel: cl6 -> Partial-EMT ; cl13 + Fibroblast(8) -> Mesenchymal ; rest unchanged ----
m <- c("0"="Epithelial_Basal","2"="Epithelial_Basal","4"="Epithelial_Basal","7"="Epithelial_Basal",
       "1"="Epithelial_Luminal","3"="Epithelial_Luminal","11"="Epithelial_Luminal","12"="Epithelial_Luminal","15"="Epithelial_Luminal",
       "6"="Partial-EMT","13"="Mesenchymal","8"="Mesenchymal",
       "9"="Macrophage_Myeloid","10"="Macrophage_Myeloid","14"="Endothelial")
lev <- c("Epithelial_Luminal","Epithelial_Basal","Partial-EMT","Mesenchymal","Endothelial","Macrophage_Myeloid")
so$celltype <- factor(unname(m[as.character(so$seurat_clusters)]), levels=lev)
pal <- c(Epithelial_Luminal="#1f78b4", Epithelial_Basal="#a6cee3", `Partial-EMT`="#fb9a99",
         Mesenchymal="#e31a1c", Endothelial="#33a02c", Macrophage_Myeloid="#6a3d9a")
saveRDS(so, file.path(OUT, "npp53_clean_labeled.RDS"))
write.csv(data.frame(cell=colnames(so), celltype=so$celltype, cluster=so$seurat_clusters),
          file.path(OUT, "npp53_celltype_cells.csv"), row.names=FALSE)

# ---- UMAPs by celltype, method-only subtitle ----
methnm <- c(umap.uncorr="Uncorrected", umap.harmony="Harmony", umap.mnn="fastMNN")
for (r in names(methnm)) {
  p <- DimPlot(so, reduction=r, group.by="celltype", cols=pal, raster=FALSE) + ggtitle(methnm[[r]])
  ggsave(file.path(OUT, paste0("UMAP_", r, "_celltype.png")), p, width=8.6, height=6.2, dpi=150)
}
# marker dotplot by celltype
feats <- c("Epcam","Krt8","Krt18","Pbsn","Ar","Krt5","Krt14","Trp63","Vim","Col3a1","Zeb1","Fn1",
           "Dcn","Lum","Pdgfra","Pecam1","Cldn5","Cdh5","Lyz2","C1qa","Cd68","Csf1r")
feats <- feats[feats %in% rownames(so)]
p <- DotPlot(so, features=feats, group.by="celltype") + RotatedAxis() +
     theme(axis.text.x=element_text(size=8)) + ggtitle("Marker genes by cell type")
ggsave(file.path(OUT, "DotPlot_celltype_markers.png"), p, width=11, height=4.5, dpi=150)

# ---- tracing features ----
trc <- intersect(c("GFP","mCherry","Cas9","Cas9-full","Cre","Cre-full"), rownames(so))
cat("tracing features:", paste(trc, collapse=", "), "\n")
# dotplot by cluster and by celltype
p <- DotPlot(so, features=trc, group.by="seurat_clusters") + RotatedAxis() +
     ggtitle("Tracing/reporter features by cluster")
ggsave(file.path(OUT, "DotPlot_tracing_by_cluster.png"), p, width=7.5, height=6, dpi=150)
p <- DotPlot(so, features=trc, group.by="celltype") + RotatedAxis() +
     ggtitle("Tracing/reporter features by cell type")
ggsave(file.path(OUT, "DotPlot_tracing_by_celltype.png"), p, width=7, height=4, dpi=150)
# featureplots GFP + mCherry on Harmony
for (f in intersect(c("GFP","mCherry"), trc)) {
  p <- FeaturePlot(so, features=f, reduction="umap.harmony", raster=FALSE, order=TRUE) + ggtitle(paste("Harmony —", f))
  ggsave(file.path(OUT, paste0("UMAP_harmony_", f, ".png")), p, width=7.5, height=6, dpi=150)
}
# per-cluster reporter positivity (pct of cells with >0 counts)
cnt <- SeuratObject::LayerData(so, assay="RNA", layer="counts")
pctpos <- function(gene) { v <- cnt[gene,]; tapply(v>0, so$seurat_clusters, mean) }
tab <- data.frame(cluster=levels(so$seurat_clusters))
for (g in trc) tab[[paste0("pct_",g)]] <- round(pctpos(g)[tab$cluster],3)
tab$celltype <- m[tab$cluster]
tab$reporter_pos <- (tab$pct_GFP >= 0.25) | (if("pct_mCherry" %in% names(tab)) tab$pct_mCherry >= 0.25 else FALSE)
cat("\n=== per-cluster reporter positivity (fraction of cells >0) ===\n"); print(tab)
write.csv(tab, file.path(OUT, "tracing_reporter_by_cluster.csv"), row.names=FALSE)

# ---- extract reporter+ clusters ----
pos_clusters <- tab$cluster[tab$reporter_pos]
cat("\nreporter+ clusters:", paste(pos_clusters, collapse=","), "\n")
so_rep <- so[, so$seurat_clusters %in% pos_clusters]
cat("extracted reporter+ cells:", ncol(so_rep), "/", ncol(so), "\n")
cat("by celltype:\n"); print(table(so_rep$celltype))
saveRDS(so_rep, file.path(OUT, "npp53_reporter_pos.RDS"))
write.csv(data.frame(cell=colnames(so_rep), celltype=so_rep$celltype, cluster=so_rep$seurat_clusters,
                     sample_id=so_rep$sample_id, treatment=so_rep$treatment),
          file.path(OUT, "npp53_reporter_pos_cells.csv"), row.names=FALSE)
cat("DONE relabel+trace\n")

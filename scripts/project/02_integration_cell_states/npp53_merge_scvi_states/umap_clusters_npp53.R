suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_clean_labeled.RDS"))
p <- DimPlot(so, reduction="umap.harmony", group.by="seurat_clusters", label=TRUE,
             label.size=4.5, repel=TRUE, raster=FALSE) + ggtitle("Harmony — cluster_id") +
     theme(legend.position="right")
ggsave(file.path(OUT, "UMAP_harmony_clusterid.png"), p, width=8.6, height=6.2, dpi=150)
cat("DONE cluster umap\n")

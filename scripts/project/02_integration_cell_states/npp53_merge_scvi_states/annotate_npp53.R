# Cluster the Harmony-integrated Npp53 object + putative cell-type annotation by canonical marker scoring.
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_merged_integrated.RDS"))
DefaultAssay(so) <- "RNA"

# cluster on Harmony
so <- FindNeighbors(so, reduction = "harmony", dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 0.5, verbose = FALSE)

# canonical mouse-prostate marker sets (broad compartments)
mk <- list(
  Epithelial_Luminal = c("Epcam","Krt8","Krt18","Krt19","Cdh1","Pbsn","Nkx3-1","Ar","Krt7","Cldn3"),
  Epithelial_Basal   = c("Krt5","Krt14","Trp63","Krt15"),
  Mesenchymal_EMT    = c("Vim","Col3a1","Col1a1","Zeb1","Zeb2","Fn1","Igfbp2","Twist1","Snai2"),
  Fibroblast         = c("Dcn","Lum","Pdgfra","Gsn","Col1a2","Sfrp2"),
  SmoothMuscle       = c("Acta2","Myh11","Tagln","Myl9","Des"),
  Endothelial        = c("Pecam1","Cldn5","Cdh5","Flt1","Emcn","Egfl7"),
  Macrophage_Myeloid = c("Lyz2","C1qa","C1qb","C1qc","Adgre1","Cd68","Itgam","Csf1r"),
  Tcell_NK           = c("Cd3e","Cd3d","Cd3g","Cd8a","Nkg7","Gzmb","Klrb1c"),
  Bcell              = c("Cd79a","Cd79b","Ms4a1","Igkc"),
  Neuroendocrine     = c("Chga","Chgb","Syp","Ncam1","Ascl1")
)
mk <- lapply(mk, function(g) intersect(g, rownames(so)))
mk <- mk[sapply(mk, length) >= 2]
cat("compartments scored:", paste(names(mk), collapse=", "), "\n")

# deterministic module score = mean of per-gene z-scored normalized expression (no AddModuleScore binning)
d <- SeuratObject::LayerData(so, assay = "RNA", layer = "data")
zrow <- function(m){ m <- as.matrix(m); mu <- rowMeans(m); s <- apply(m, 1, sd); s[s==0] <- 1; (m - mu)/s }
scores <- sapply(names(mk), function(nm) colMeans(zrow(d[mk[[nm]], , drop=FALSE])))  # cells x compartments
rownames(scores) <- colnames(so)
for (nm in colnames(scores)) so@meta.data[[paste0("sc_", nm)]] <- scores[, nm]

# per-cluster mean score -> argmax label
cl  <- so$seurat_clusters
mat <- apply(scores, 2, function(x) tapply(x, cl, mean))   # clusters x compartments
cluster_label <- setNames(colnames(mat)[max.col(mat, ties.method="first")], rownames(mat))
so$putative <- factor(unname(cluster_label[as.character(cl)]))
cat("\n=== cluster -> putative ===\n"); print(data.frame(cluster=rownames(mat), putative=cluster_label))
cat("\n=== putative cell-type counts ===\n"); print(table(so$putative))
cat("\n=== putative x batch ===\n"); print(table(so$putative, so$batch))

saveRDS(so, file.path(OUT, "npp53_annotated.RDS"))
write.csv(data.frame(cell=colnames(so), cluster=cl, putative=so$putative, batch=so$batch,
                     sample_id=so$sample_id, treatment=so$treatment),
          file.path(OUT, "npp53_putative_celltype.csv"), row.names=FALSE)

pal <- c(Epithelial_Luminal="#2166ac", Epithelial_Basal="#4393c3", Mesenchymal_EMT="#b4232a",
         Fibroblast="#E69F00", SmoothMuscle="#762a83", Endothelial="#1b7837",
         Macrophage_Myeloid="#9970ab", Tcell_NK="#f781bf", Bcell="#a6761d", Neuroendocrine="#000000")
p1 <- DimPlot(so, reduction="umap.harmony", group.by="putative", cols=pal, raster=FALSE) +
      ggtitle("Npp53 (Harmony) — putative cell type")
ggsave(file.path(OUT,"UMAP_harmony_putativeCellType.png"), p1, width=9, height=6.5, dpi=150)
p2 <- DimPlot(so, reduction="umap.harmony", group.by="seurat_clusters", label=TRUE, raster=FALSE) +
      ggtitle("Npp53 (Harmony) — clusters")
ggsave(file.path(OUT,"UMAP_harmony_clusters.png"), p2, width=8.5, height=6.5, dpi=150)
cat("DONE annotate ->", OUT, "\n")

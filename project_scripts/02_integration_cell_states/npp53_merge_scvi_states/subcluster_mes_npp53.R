suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_clean_labeled.RDS")); DefaultAssay(so) <- "RNA"
mes <- so[, so$seurat_clusters %in% c("8","13")]   # Mesenchymal compartment
cat("Mesenchymal cells:", ncol(mes), "\n"); print(table(mes$seurat_clusters))

# ---- Harmony subclustering ----
mes <- FindNeighbors(mes, reduction="harmony", dims=1:30, verbose=FALSE)
mes <- FindClusters(mes, resolution=0.5, verbose=FALSE); mes$sub_harmony <- mes$seurat_clusters
mes <- RunUMAP(mes, reduction="harmony", dims=1:30, reduction.name="umap.sub", verbose=FALSE)

# ---- scVI subclustering (import latent) ----
lat <- read.csv(file.path(OUT, "scvi_latent_clean.csv"), row.names=1)
lat <- as.matrix(lat[colnames(mes), , drop=FALSE])
mes[["scvi"]] <- CreateDimReducObject(embeddings=lat, key="scvi_", assay=DefaultAssay(mes))
mes <- FindNeighbors(mes, reduction="scvi", dims=1:ncol(lat), verbose=FALSE)
mes <- FindClusters(mes, resolution=0.5, verbose=FALSE); mes$sub_scvi <- mes$seurat_clusters
mes <- RunUMAP(mes, reduction="scvi", dims=1:ncol(lat), reduction.name="umap.subscvi", verbose=FALSE)

# ---- reporter + host/fibroblast markers per subcluster ----
cnt <- SeuratObject::LayerData(mes, assay="RNA", layer="counts")
dat <- SeuratObject::LayerData(mes, assay="RNA", layer="data")
mk  <- intersect(c("GFP","mCherry","Dcn","Lum","Pdgfra","Pdgfrb","Col1a1","Col3a1","Acta2","Gsn","Sfrp2",
                   "Vim","Zeb1","Fn1","Epcam","Krt8","Krt18","Ptprc","Pecam1"), rownames(mes))
report <- function(grp, gname) {
  cat("\n===== subclusters by", gname, "=====\n")
  tb <- data.frame(subcluster=levels(grp), n=as.integer(table(grp)))
  tb$pct_GFP     <- round(tapply(cnt["GFP",]>0, grp, mean)[tb$subcluster],3)
  tb$pct_mCherry <- round(tapply(cnt["mCherry",]>0, grp, mean)[tb$subcluster],3)
  for (g in intersect(c("Dcn","Lum","Pdgfra","Col1a1","Acta2","Vim","Epcam","Ptprc"), mk))
    tb[[paste0("m_",g)]] <- round(tapply(dat[g,], grp, mean)[tb$subcluster],2)
  print(tb); write.csv(tb, file.path(OUT, paste0("mes_subcluster_",gname,".csv")), row.names=FALSE); tb
}
th <- report(mes$sub_harmony, "harmony"); ts <- report(mes$sub_scvi, "scvi")

# ---- plots ----
pal <- scales::hue_pal()(max(nlevels(mes$sub_harmony), nlevels(mes$sub_scvi)))
save_umap <- function(red, grp, ttl, fn) {
  p1 <- DimPlot(mes, reduction=red, group.by=grp, label=TRUE, raster=FALSE) + ggtitle(ttl)
  p2 <- FeaturePlot(mes, features="GFP", reduction=red, raster=FALSE, order=TRUE) + ggtitle("GFP")
  p3 <- FeaturePlot(mes, features="mCherry", reduction=red, raster=FALSE, order=TRUE) + ggtitle("mCherry")
  p4 <- FeaturePlot(mes, features="Dcn", reduction=red, raster=FALSE, order=TRUE) + ggtitle("Dcn (fibroblast)")
  ggsave(file.path(OUT, fn), (p1|p2)/(p3|p4), width=12, height=9, dpi=140)
}
save_umap("umap.sub",     "sub_harmony", "Mesenchymal subclusters — Harmony (res 0.5)", "mes_subcluster_harmony.png")
save_umap("umap.subscvi", "sub_scvi",    "Mesenchymal subclusters — scVI (res 0.5)",   "mes_subcluster_scvi.png")
# marker dotplot (harmony subclusters)
p <- DotPlot(mes, features=mk, group.by="sub_harmony") + RotatedAxis() + ggtitle("Mesenchymal subclusters (Harmony) — markers")
ggsave(file.path(OUT, "mes_subcluster_dotplot.png"), p, width=11, height=4.5, dpi=140)

saveRDS(mes, file.path(OUT, "npp53_mesenchymal_sub.RDS"))
write.csv(data.frame(cell=colnames(mes), orig_cluster=mes$seurat_clusters,
                     sub_harmony=mes$sub_harmony, sub_scvi=mes$sub_scvi),
          file.path(OUT, "mes_subcluster_assign.csv"), row.names=FALSE)
cat("\nDONE subcluster mes\n")

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr); library(patchwork) })
M <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
P <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_plots"
so <- readRDS(file.path(M,"npp53_edited.RDS")); DefaultAssay(so) <- "RNA"
nd <- ncol(Embeddings(so,"scvi"))
so <- FindNeighbors(so, reduction="scvi", dims=1:nd, verbose=FALSE)

# canonical prostate + tracing marker panel (grouped by identity)
feats <- intersect(c(
  # luminal core
  "Krt8","Krt18","Epcam","Cldn3","Ppp1r1b","Pigr",
  # luminal AR axis
  "Ar","Nkx3-1","Pbsn","Fkbp5","Tmprss2","Msmb",
  # luminal club / hillock / squamous / distal / agr2
  "Rbp2","Cst6","Tacstd2","Psca","Aldh1a3","Krt79","Gabrp","Krt4","Krt6b","Krt10","Sprr1b","Defb42","Eppin","Spink10","Cuzd1","Agr2","Reg3g","Cxcl15","Cxcl17",
  # ionocyte
  "Foxi1","Atp6v1g3","Aqp6","Atp6v0b",
  # basal
  "Krt5","Krt14","Trp63","Krt15","Serpina3n","Fst","Ndufa4l2","Col17a1","Ovol1","Dsg3",
  # cycling
  "Mki67","Top2a","Ccnb1","Aspm","Pbk","Aurkb","Ccna2",
  # EMT / mesenchymal
  "Vim","Cdh6","Vcam1","Cd44","Zeb1","Zeb2","Fn1","Twist1","Sox2","Snai2",
  # mature mesenchymal
  "Col3a1","Col1a1","Col5a2","Vcan","Dpt","Bgn","Myoc","Igfbp2","Bmp6","Lrp2","Mmp13","Hoxb13","Lrrc15","Aspn","Col11a1","Mfap4","Acan",
  # tracing
  "GFP","mCherry"), rownames(so))

# for each candidate res: cluster, save DotPlot BY CLUSTER (raw)
for (r in c(0.4, 0.5, 0.8, 1.0)) {
  so <- FindClusters(so, resolution=r, verbose=FALSE)
  cl <- so$seurat_clusters
  cat(sprintf("res=%.1f k=%d\n", r, length(unique(cl))))
  Idents(so) <- "seurat_clusters"
  p <- DotPlot(so, features=feats, group.by="seurat_clusters", cluster.idents=TRUE) +
       RotatedAxis() + theme(axis.text.x=element_text(size=6.5)) +
       ggtitle(paste0("scVI res",r," (per cluster; hierarchical order)"))
  ggsave(file.path(P, paste0("edited_scvi_res",r,"_dotplot_cluster.png")), p, width=20, height=6+0.16*length(unique(cl)), dpi=140)
}

# ---- FindAllMarkers at the biologically-recommended res 0.5 ----
so <- FindClusters(so, resolution=0.5, verbose=FALSE)
so$scvi_res05 <- so$seurat_clusters
Idents(so) <- "scvi_res05"
mk <- FindAllMarkers(so, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.4, verbose=FALSE)
write.csv(mk, file.path(M,"scvi_res05_all_markers.csv"), row.names=FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n=15) %>% summarise(genes=paste(gene, collapse=", "))
cat("=== top markers per scVI res0.5 cluster ===\n"); print(as.data.frame(top), right=FALSE)
write.csv(as.data.frame(top), file.path(M,"scvi_res05_top15_markers.csv"), row.names=FALSE)

# also do res 0.8 (current primary) for completeness
so <- FindClusters(so, resolution=0.8, verbose=FALSE)
so$scvi_res08 <- so$seurat_clusters
Idents(so) <- "scvi_res08"
mk8 <- FindAllMarkers(so, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.4, verbose=FALSE)
write.csv(mk8, file.path(M,"scvi_res08_all_markers.csv"), row.names=FALSE)
top8 <- mk8 %>% group_by(cluster) %>% slice_max(avg_log2FC, n=15) %>% summarise(genes=paste(gene, collapse=", "))
cat("\n=== top markers per scVI res0.8 cluster ===\n"); print(as.data.frame(top8), right=FALSE)
write.csv(as.data.frame(top8), file.path(M,"scvi_res08_top15_markers.csv"), row.names=FALSE)

saveRDS(so, file.path(M,"npp53_edited.RDS"))
cat("DONE dotplots + markers\n")

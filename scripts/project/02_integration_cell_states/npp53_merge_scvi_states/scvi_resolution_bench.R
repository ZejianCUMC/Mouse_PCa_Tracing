suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(clustree)
  library(cluster); library(mclust); library(pheatmap); library(dplyr) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT,"npp53_edited.RDS")); DefaultAssay(so) <- "RNA"
nd  <- ncol(Embeddings(so,"scvi"))
so  <- FindNeighbors(so, reduction="scvi", dims=1:nd, verbose=FALSE)
grid <- c(0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2)
for (r in grid){ so <- FindClusters(so, resolution=r, verbose=FALSE); so[[paste0("scvi_res.",r)]] <- so$seurat_clusters }

# (A) clustree stability across resolutions
pc <- clustree(so, prefix="scvi_res.") + ggtitle("scVI clustering tree")
ggsave(file.path(OUT,"edited_scvi_clustree.png"), pc, width=9, height=11, dpi=140)

# (B) silhouette (scVI latent) + ARI vs Harmony 18-state — per resolution
set.seed(1); ix <- sort(sample(seq_len(ncol(so)), 6000)); idx <- colnames(so)[ix]
D <- dist(Embeddings(so,"scvi")[ix,])
bench <- data.frame(res=grid, k=NA_integer_, mean_sil=NA_real_, ARI_vs_Harmony=NA_real_)
for (i in seq_along(grid)) {
  lab <- so@meta.data[[paste0("scvi_res.",grid[i])]]
  cl  <- as.integer(as.factor(lab[ix]))
  sl  <- silhouette(cl, D)
  bench$k[i]        <- length(unique(lab))
  bench$mean_sil[i] <- if (is.matrix(sl)) mean(sl[,3]) else NA
  bench$ARI_vs_Harmony[i] <- adjustedRandIndex(lab, so$sub08)
}
cat("=== resolution benchmark (scVI latent) ===\n"); print(bench, row.names=FALSE)
write.csv(bench, file.path(OUT,"scvi_resolution_bench.csv"), row.names=FALSE)
sc <- max(bench$mean_sil)/max(bench$k)
p1 <- ggplot(bench, aes(res)) +
  geom_line(aes(y=mean_sil,color="mean silhouette")) + geom_point(aes(y=mean_sil,color="mean silhouette")) +
  geom_line(aes(y=ARI_vs_Harmony,color="ARI vs Harmony")) + geom_point(aes(y=ARI_vs_Harmony,color="ARI vs Harmony")) +
  geom_text(aes(y=mean_sil,label=paste0("k=",k)), vjust=-0.9, size=3) +
  labs(y="metric", color="", title="Resolution selection — scVI latent") + theme_bw()
ggsave(file.path(OUT,"edited_scvi_resolution_metrics.png"), p1, width=8.5, height=5, dpi=140)

# (C) annotate scVI clusters (choose res with k nearest 18 for like-for-like) via signature argmax
rbest <- grid[which.min(abs(bench$k-18))]
cat("\nscVI res chosen for subtype comparison (k nearest 18):", rbest, "k=", bench$k[which.min(abs(bench$k-18))], "\n")
so$scvi_best <- so@meta.data[[paste0("scvi_res.",rbest)]]
sig <- list(
  Luminal=c("Pigr","Ppp1r1b","Krt8","Krt18"), Club=c("Rbp2","Cst6","Aldh1a3","Krt79"),
  Squamous=c("Krt6b","Krt10","Sprr1b"), Ionocyte=c("Foxi1","Atp6v1g3","Aqp6"),
  Basal=c("Krt5","Krt14","Trp63"), Cycling=c("Mki67","Top2a","Ccnb1"),
  PartialEMT=c("Cdh6","Vcam1"), Mesenchymal=c("Col5a2","Vcan","Twist1","Lrrc15","Igfbp2"))
sig <- lapply(sig,function(g) intersect(g,rownames(so))); d<-SeuratObject::LayerData(so,assay="RNA",layer="data")
zrow<-function(m){m<-as.matrix(m);mu<-rowMeans(m);s<-apply(m,1,sd);s[s==0]<-1;(m-mu)/s}
sco <- sapply(names(sig),function(nm) colMeans(zrow(d[sig[[nm]], ,drop=FALSE])))
hint <- tapply(1:ncol(so), so$scvi_best, function(ix) names(sig)[which.max(colMeans(sco[ix,,drop=FALSE]))])
cat("\n=== scVI cluster -> signature hint ===\n"); print(data.frame(scvi_cluster=names(hint), hint=unlist(hint)))

# (D) concordance heatmap: Harmony 18-state (rows) x scVI clusters (cols), row-normalised
ct <- table(so$substate, so$scvi_best); ctn <- sweep(ct,1,pmax(rowSums(ct),1),"/")
png(file.path(OUT,"edited_scvi_vs_harmony_heatmap.png"), width=1500, height=1200, res=150)
pheatmap(ctn, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=round(ct,0),
         number_format="%.0f", fontsize_number=7, main=paste0("Harmony 18-state (rows) vs scVI res",rbest," (cols) — row fraction"),
         color=colorRampPalette(c("white","#4575b4"))(50))
dev.off()
cat("\nOverall ARI (scVI res",rbest,"vs Harmony sub08):", round(adjustedRandIndex(so$scvi_best,so$sub08),3),"\n")

# scVI-clustered UMAP at chosen res
p2 <- DimPlot(so, reduction="umap.scvi", group.by="scvi_best", raster=FALSE, label=TRUE) + ggtitle(paste0("scVI res",rbest))
ggsave(file.path(OUT,"edited_umap_scvi_clusters.png"), p2, width=9.5, height=6.4, dpi=150)
saveRDS(so, file.path(OUT,"npp53_edited.RDS"))
cat("DONE scvi resolution bench\n")

# Bootstrap cluster stability across scVI resolutions on a step-0.1 grid.
suppressPackageStartupMessages({ library(Seurat); library(mclust); library(ggplot2) })
M <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
P <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_plots"
so <- readRDS(file.path(M,"npp53_edited.RDS")); DefaultAssay(so) <- "RNA"
nd <- ncol(Embeddings(so,"scvi"))
so <- FindNeighbors(so, reduction="scvi", dims=1:nd, verbose=FALSE)

grid <- seq(0.1, 1.2, 0.1)
full <- list()
for (r in grid){
  so <- FindClusters(so, resolution=r, verbose=FALSE)
  v <- as.character(so$seurat_clusters); names(v) <- colnames(so)
  full[[as.character(r)]] <- v
}
B <- 20; A <- matrix(NA_real_, length(grid), B)
set.seed(1); samp <- lapply(1:B, function(x) sample(colnames(so), round(0.8*ncol(so))))
for (bi in 1:B) {
  cells <- samp[[bi]]
  sub <- so[, cells]; sub <- FindNeighbors(sub, reduction="scvi", dims=1:nd, verbose=FALSE)
  cn <- colnames(sub)
  for (j in seq_along(grid)) {
    sub <- FindClusters(sub, resolution=grid[j], verbose=FALSE)
    A[j,bi] <- adjustedRandIndex(as.character(sub$seurat_clusters), full[[as.character(grid[j])]][cn])
  }
  cat("bootstrap", bi, "done\n")
}
st <- data.frame(res=grid, k=sapply(grid,function(r) length(unique(full[[as.character(r)]]))),
                 mean_ari=rowMeans(A), se_ari=apply(A,1,sd)/sqrt(B))
cat("=== bootstrap stability (scVI, step 0.1, B=20) ===\n"); print(st, row.names=FALSE)
write.csv(st, file.path(M,"scvi_bootstrap_stability.csv"), row.names=FALSE)

p <- ggplot(st, aes(res, mean_ari)) +
  geom_line(color="#2166ac") + geom_point(color="#2166ac", size=3) +
  geom_errorbar(aes(ymin=mean_ari-se_ari, ymax=mean_ari+se_ari), width=.02, color="#2166ac") +
  geom_text(aes(label=paste0("k=",k)), vjust=-1.2, size=3.2) +
  scale_x_continuous(breaks=seq(0.1, 1.2, 0.1), limits=c(0.05, 1.25)) +
  labs(y="bootstrap ARI (reproducibility)", x="scVI clustering resolution",
       title="Bootstrap cluster stability across scVI resolutions",
       subtitle="B=20, 80% subsample; error bars = SE (step 0.1)") + theme_classic()
ggsave(file.path(P,"edited_scvi_bootstrap_stability.png"), p, width=9, height=5, dpi=150)
cat("DONE\n")

# chooseR-style quantitative resolution selection: subsample 80%, re-cluster, measure
# reproducibility of each resolution's clustering by ARI vs the full-data clustering; 1-SE rule.
suppressPackageStartupMessages({ library(Seurat); library(mclust); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT,"npp53_edited.RDS")); DefaultAssay(so) <- "RNA"
nd  <- ncol(Embeddings(so,"scvi")); grid <- c(0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2)
so  <- FindNeighbors(so, reduction="scvi", dims=1:nd, verbose=FALSE)
full <- list()
for (r in grid){ so <- FindClusters(so, resolution=r, verbose=FALSE); v <- as.character(so$seurat_clusters); names(v) <- colnames(so); full[[as.character(r)]] <- v }
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
imax <- which.max(st$mean_ari); thr <- st$mean_ari[imax]-st$se_ari[imax]
st$within_1SE <- st$mean_ari >= thr
pick <- st$res[max(which(st$within_1SE))]          # most granular resolution within 1 SE of the max
cat("=== bootstrap stability (scVI latent, B=20, 80% subsample) ===\n"); print(st, row.names=FALSE)
cat("\nmax-stability res:", st$res[imax], "| 1-SE threshold:", round(thr,3),
    "| chosen (most granular within 1SE):", pick, "\n")
write.csv(st, file.path(OUT,"scvi_bootstrap_stability.csv"), row.names=FALSE)
p <- ggplot(st, aes(res, mean_ari)) +
  geom_hline(yintercept=thr, linetype="dashed", color="#d05252") +
  geom_line(color="#2166ac") + geom_point(aes(color=within_1SE), size=3) +
  geom_errorbar(aes(ymin=mean_ari-se_ari, ymax=mean_ari+se_ari), width=.02, color="#2166ac") +
  geom_text(aes(label=paste0("k=",k)), vjust=-1.1, size=3) +
  scale_color_manual(values=c(`TRUE`="#4caf50",`FALSE`="#999999")) +
  labs(y="bootstrap ARI (reproducibility)", x="scVI resolution",
       title="Bootstrap cluster stability (chooseR-style, 1-SE rule)",
       subtitle=paste0("dashed = max - 1SE; chosen res = ",pick)) + theme_bw()
ggsave(file.path(OUT,"edited_scvi_bootstrap_stability.png"), p, width=8.5, height=5, dpi=140)
cat("DONE bootstrap stability\n")

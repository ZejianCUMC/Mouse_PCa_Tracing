# Per-cluster bootstrap Jaccard at scVI res 0.4 AND res 0.8 -> two-panel Fig 4.2.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
M <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
P <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_plots"
so <- readRDS(file.path(M,"npp53_edited.RDS")); DefaultAssay(so) <- "RNA"
nd <- ncol(Embeddings(so,"scvi"))
so <- FindNeighbors(so, reduction="scvi", dims=1:nd, verbose=FALSE)

# use the same "state_res0p4" and "state_res0p8" columns saved by prior job
compute_jac <- function(r) {
  so <- FindClusters(so, resolution=r, verbose=FALSE)
  full <- setNames(as.character(so$seurat_clusters), colnames(so))
  cls_ids <- sort(unique(full))
  B <- 20; set.seed(1)
  jac <- matrix(NA_real_, length(cls_ids), B, dimnames=list(cls_ids, NULL))
  for (b in 1:B) {
    cells <- sample(colnames(so), round(0.8*ncol(so)))
    sub <- so[, cells]
    sub <- FindNeighbors(sub, reduction="scvi", dims=1:nd, verbose=FALSE)
    sub <- FindClusters(sub, resolution=r, verbose=FALSE)
    bc <- setNames(as.character(sub$seurat_clusters), colnames(sub))
    for (cid in cls_ids) {
      A <- names(full)[full==cid & names(full) %in% names(bc)]
      if (!length(A)) next
      best <- 0
      for (bc_id in unique(bc[A])) {
        Bset <- names(bc)[bc==bc_id]
        jj <- length(intersect(A,Bset)) / length(union(A,Bset))
        if (jj>best) best <- jj
      }
      jac[cid,b] <- best
    }
    cat(sprintf("res=%.1f boot %d done\n", r, b))
  }
  # label per cluster from putative-state assignment
  lab_col <- if (r==0.4) "state_res0p4" else if (r==0.8) "state_res0p8" else NULL
  if (is.null(lab_col) || !(lab_col %in% colnames(so@meta.data))) lab_col <- NULL
  lab_by_cl <- sapply(cls_ids, function(cid) {
    ix <- which(full==cid)
    if (!is.null(lab_col)) { x <- so@meta.data[[lab_col]][ix]; names(sort(table(x), decreasing=TRUE))[1] } else "?"
  })
  data.frame(res=paste0("scVI res ",r), cluster=cls_ids, label=lab_by_cl,
             n=as.integer(table(full)[cls_ids]),
             mean_jaccard=rowMeans(jac, na.rm=TRUE),
             sd_jaccard=apply(jac,1,sd,na.rm=TRUE))
}
df04 <- compute_jac(0.4); df08 <- compute_jac(0.8)
df <- rbind(df04, df08)
write.csv(df, file.path(M,"scvi_per_cluster_jaccard_04_08.csv"), row.names=FALSE)

make_panel <- function(sub) {
  sub <- sub[order(sub$mean_jaccard),]
  sub$labshow <- paste0(sub$label," (cl",sub$cluster,", n=",sub$n,")")
  sub$labshow <- factor(sub$labshow, levels=sub$labshow)
  ggplot(sub, aes(labshow, mean_jaccard, fill=mean_jaccard>=0.75)) +
    geom_col() + geom_errorbar(aes(ymin=mean_jaccard-sd_jaccard, ymax=mean_jaccard+sd_jaccard), width=.35, color="#333") +
    geom_hline(yintercept=0.75, linetype="dashed", color="#d05252") +
    geom_hline(yintercept=0.60, linetype="dotted", color="#e8a23d") +
    scale_fill_manual(values=c(`TRUE`="#4caf50", `FALSE`="#a08a4c"), guide="none") +
    coord_flip(ylim=c(0,1.02)) +
    labs(x="", y="mean bootstrap Jaccard", title=unique(sub$res)) + theme_classic()
}
p4 <- make_panel(df04); p8 <- make_panel(df08)
composite <- (p4 | p8) + plot_annotation(title="Per-cluster bootstrap Jaccard stability",
                                         subtitle=">=0.75 stable (green), >=0.60 sensible (dotted), <0.60 unstable/dissolving; B=20, 80% subsample")
ggsave(file.path(P,"edited_percluster_jaccard.png"), composite, width=15, height=6.6, dpi=150)
cat("DONE per-cluster Jaccard (res 0.4 + 0.8)\n")

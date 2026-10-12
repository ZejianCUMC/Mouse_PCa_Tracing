suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr); library(cluster); library(scales) })
M <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
P <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_plots"
so <- readRDS(file.path(M,"npp53_edited.RDS")); DefaultAssay(so) <- "RNA"
nd <- ncol(Embeddings(so,"scvi"))

# ---- 1. Fig 4.1 re-render: theme_classic, x 0-1.4 step 0.2, no 1-SE line ----
st <- read.csv(file.path(M,"scvi_bootstrap_stability.csv"))
p41 <- ggplot(st, aes(res, mean_ari)) +
  geom_line(color="#2166ac") + geom_point(color="#2166ac", size=3) +
  geom_errorbar(aes(ymin=mean_ari-se_ari, ymax=mean_ari+se_ari), width=.02, color="#2166ac") +
  geom_text(aes(label=paste0("k=",k)), vjust=-1.2, size=3.2) +
  scale_x_continuous(breaks=seq(0, 1.4, 0.2), limits=c(0.1, 1.3)) +
  labs(y="bootstrap ARI (reproducibility)", x="scVI clustering resolution",
       title="Bootstrap cluster stability across scVI resolutions",
       subtitle="B=20, 80% subsample; error bars = SE") + theme_classic()
ggsave(file.path(P,"edited_scvi_bootstrap_stability.png"), p41, width=8.5, height=5, dpi=150)

# ---- 2. §4.4 silhouette per-cluster diagram at candidate resolutions ----
# do it on scVI latent, cluster at res 0.4 and 0.8, plot silhouette-per-cell histogram
so <- FindNeighbors(so, reduction="scvi", dims=1:nd, verbose=FALSE)
set.seed(1); ix <- sort(sample(seq_len(ncol(so)), 6000))
Dm <- dist(Embeddings(so,"scvi")[ix,])
sil_df <- list()
for (r in c(0.4, 0.8)) {
  so <- FindClusters(so, resolution=r, verbose=FALSE)
  cl <- as.integer(as.factor(so$seurat_clusters[ix]))
  s <- silhouette(cl, Dm)
  if (!is.matrix(s)) next
  df <- data.frame(cluster=factor(s[,"cluster"]), sil_width=s[,"sil_width"], res=paste0("res ",r))
  sil_df[[as.character(r)]] <- df
  cat(sprintf("res=%.1f mean_sil=%.3f\n", r, mean(s[,"sil_width"])))
}
sil <- do.call(rbind, sil_df)
p44 <- ggplot(sil, aes(sil_width)) +
  geom_histogram(binwidth=.02, fill="#4caf50", color="#333", alpha=.85) +
  facet_wrap(~res, ncol=2) +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x="silhouette width (per cell)", y="cell count",
       title="Silhouette width distribution on scVI latent",
       subtitle="Values ~0 = continuous manifold (no discrete cluster boundaries)") + theme_classic()
ggsave(file.path(P,"edited_scvi_silhouette.png"), p44, width=10, height=4.5, dpi=150)

# ---- 3. res 0.4 signature-argmax annotation + dotplot by putative state ----
# same signature list used previously
sig <- list(
  `Luminal-Secretory`=c("Pigr","Ppp1r1b","Krt8","Krt18","Cldn3"),
  `Luminal-Club`     =c("Rbp2","Cst6","Tacstd2","Psca"),
  `Luminal-Distal`   =c("Defb42","Eppin","Spink10","Defb28","Cuzd1"),
  `Luminal-Agr2`     =c("Agr2","Reg3g","Cxcl15","Cxcl17"),
  `Luminal-Club-Hillock`=c("Aldh1a3","Krt79","Gabrp","Krt4","Serping1"),
  `Luminal-Squamous` =c("Krt6b","Krt10","Sprr1b","Krtdap","Calm4"),
  `Ionocyte`         =c("Foxi1","Atp6v1g3","Aqp6","Atp6v0b"),
  `Basal`            =c("Krt5","Krt14","Trp63"),
  `Basal-Inflammatory`=c("Serpina3n","Fst","Igfbp5","Apoe","Id3"),
  `Basal-Hypoxic`    =c("Ndufa4l2","Col17a1","Cxcl14","Smoc2"),
  `Basal-Differentiating`=c("Ovol1","Dsg3","Il24","Cntfr","Trem2"),
  `Cycling`          =c("Mki67","Top2a","Pbk","Aurkb","Ccnb1","Aspm"),
  `Partial-EMT`      =c("Cdh6","Vcam1","Gdnf","Vim","Cd44"),
  `Mes-Igfbp2`       =c("Igfbp2","Bmp6","Lrp2","Pcolce2"),
  `Mes-Collagen`     =c("Col5a2","Vcan","Dpt","Bgn","Myoc"),
  `Mes-EMT-TF`       =c("Twist1","Sox2","Mmp13","Hoxb13","Col24a1"),
  `Mes-CAF-like`     =c("Lrrc15","Aspn","Col11a1","Mfap4","Acan"))
sig <- lapply(sig, function(g) intersect(g, rownames(so))); sig <- sig[sapply(sig,length)>=2]
d <- SeuratObject::LayerData(so, assay="RNA", layer="data")
zrow <- function(m){ m<-as.matrix(m); mu<-rowMeans(m); s<-apply(m,1,sd); s[s==0]<-1; (m-mu)/s }
sco <- sapply(names(sig), function(nm) colMeans(zrow(d[sig[[nm]], , drop=FALSE])))

# canonical marker panel
feats <- intersect(c(
  "Krt8","Krt18","Epcam","Ppp1r1b","Pigr","Ar","Nkx3-1","Pbsn","Fkbp5",
  "Rbp2","Cst6","Aldh1a3","Krt79","Gabrp","Krt6b","Krt10","Sprr1b","Defb42","Eppin","Cuzd1","Agr2","Reg3g","Cxcl15","Cxcl17","Foxp2","Pitx1","Hoxa9","Krt13",
  "Foxi1","Atp6v1g3","Aqp6","Krt5","Krt14","Trp63","Serpina3n","Fst","Ndufa4l2","Col17a1","Ovol1","Dsg3","Ccl2","Spp1",
  "Mki67","Top2a","Ccnb1","Aspm",
  "Vim","Cdh6","Vcam1","Cd44",
  "Col3a1","Col1a1","Col5a2","Vcan","Dpt","Bgn","Igfbp2","Bmp6","Mmp13","Hoxb13","Lrrc15","Aspn","Col11a1","Prrx1","Twist1","Sox2",
  "GFP","mCherry"), rownames(so))

annotate_and_dotplot <- function(r, tag) {
  so <- FindClusters(so, resolution=r, verbose=FALSE)
  cl <- so$seurat_clusters
  mat <- apply(sco, 2, function(x) tapply(x, cl, mean))
  hint <- setNames(colnames(mat)[max.col(mat, ties.method="first")], rownames(mat))
  # collapse repeated labels but note pooling
  ann <- factor(unname(hint[as.character(cl)]),
                levels=intersect(names(sig), unique(unname(hint))))
  colname <- paste0("state_res",tag)
  so@meta.data[[colname]] <<- ann
  cat("res=",r," clusters=",length(unique(cl))," states=",length(unique(ann)),"\n")
  # cluster -> putative state map (for memo)
  df <- data.frame(cluster=rownames(mat), putative=hint, n=as.integer(table(cl)))
  write.csv(df, file.path(M, paste0("scvi_res",r,"_clusters_to_states.csv")), row.names=FALSE)
  # dotplot by putative STATE (not cluster)
  Idents(so) <- colname
  p <- DotPlot(so, features=feats, group.by=colname, cluster.idents=TRUE) +
       RotatedAxis() + theme(axis.text.x=element_text(size=7)) +
       ggtitle(paste0("scVI res ",r," — marker expression by putative cell state (",length(unique(ann))," states)"))
  ggsave(file.path(P, paste0("edited_scvi_res",r,"_dotplot_by_state.png")), p, width=17, height=6, dpi=140)
  so
}
so <- annotate_and_dotplot(0.4, "0p4")
so <- annotate_and_dotplot(0.8, "0p8")

saveRDS(so, file.path(M,"npp53_edited.RDS"))
cat("DONE v3 updates\n")

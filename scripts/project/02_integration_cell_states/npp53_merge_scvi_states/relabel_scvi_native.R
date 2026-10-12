# Relabel Harmony states + build fully-independent scVI-native annotation for direct comparison.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(pheatmap); library(mclust); library(dplyr) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_edited.RDS")); DefaultAssay(so) <- "RNA"

# ---- 1. Relabel Harmony sub-states (add Luminal- prefix; Ionocyte-like -> Ionocyte) ----
rn <- c("Club-Hillock"="Luminal-Club-Hillock",
        "Squamous"     ="Luminal-Squamous",
        "Ionocyte-like"="Ionocyte")
oldlev <- levels(so$substate)
newlev <- ifelse(oldlev %in% names(rn), rn[oldlev], oldlev)
newvec <- as.character(so$substate)
newvec[newvec %in% names(rn)] <- rn[newvec[newvec %in% names(rn)]]
so$substate <- factor(newvec, levels=newlev)
pal_h <- setNames(c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#fdae6b","#6a51a3",
                    "#238b45","#00441b","#41ab5d","#a1d99b","#969696","#525252",
                    "#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"), newlev)
cat("=== Harmony sub-state counts (after rename) ===\n"); print(table(so$substate))

# ---- 2. scVI-native clustering + annotation (independent) ----
nd <- ncol(Embeddings(so, "scvi"))
so <- FindNeighbors(so, reduction="scvi", dims=1:nd, verbose=FALSE)
so <- FindClusters(so,  resolution=0.8, verbose=FALSE)
so$scvi_cl <- so$seurat_clusters

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
  `Cycling-S`        =c("Mki67","Top2a","Pbk","Aurkb"),
  `Cycling-G2M`      =c("Ccnb1","Aspm","Ccna2","Cenpe"),
  `Partial-EMT`      =c("Cdh6","Vcam1","Gdnf","Vim","Cd44"),
  `Mes-Igfbp2`       =c("Igfbp2","Bmp6","Lrp2","Pcolce2"),
  `Mes-Collagen`     =c("Col5a2","Vcan","Dpt","Bgn","Myoc"),
  `Mes-EMT-TF`       =c("Twist1","Sox2","Mmp13","Hoxb13","Col24a1"),
  `Mes-CAF-like`     =c("Lrrc15","Aspn","Col11a1","Mfap4","Acan"))
sig <- lapply(sig, function(g) intersect(g, rownames(so)))
sig <- sig[sapply(sig, length) >= 2]
d <- SeuratObject::LayerData(so, assay="RNA", layer="data")
zrow <- function(m){ m<-as.matrix(m); mu<-rowMeans(m); s<-apply(m,1,sd); s[s==0]<-1; (m-mu)/s }
sco <- sapply(names(sig), function(nm) colMeans(zrow(d[sig[[nm]], , drop=FALSE])))
cl  <- so$scvi_cl
mat <- apply(sco, 2, function(x) tapply(x, cl, mean))
hint <- setNames(colnames(mat)[max.col(mat, ties.method="first")], rownames(mat))
# guarantee a stable factor order matching Harmony's (unused labels dropped)
so$substate_scvi <- factor(unname(hint[as.character(cl)]),
                           levels=intersect(newlev, unique(unname(hint))))
cat("\n=== scVI-native cluster -> label ===\n")
print(data.frame(cluster=rownames(mat), label=hint, n=as.integer(table(cl))))
cat("\n=== scVI-native sub-state counts ===\n"); print(table(so$substate_scvi))
pal_s <- pal_h[levels(so$substate_scvi)]

# ---- 3. UMAPs ----
save1 <- function(red, grp, ttl, fn, pal=NULL) {
  p <- DimPlot(so, reduction=red, group.by=grp, cols=pal, raster=FALSE, label=TRUE, label.size=3, repel=TRUE) + ggtitle(ttl)
  ggsave(file.path(OUT, fn), p, width=10.5, height=6.6, dpi=150) }
save1("umap.ed",   "substate",      "Harmony res0.8",         "edited_umap_substate.png",              pal_h)
save1("umap.scvi", "substate",      "scVI UMAP",    "edited_umap_scvi_HarmLabels.png",       pal_h)
save1("umap.scvi", "substate_scvi", "scVI res0.8",        "edited_umap_scvi_scviLabels.png",       pal_s)

# ---- 4. Dotplots for both label systems ----
feats <- intersect(c("Pigr","Ppp1r1b","Rbp2","Cst6","Defb42","Eppin","Agr2","Reg3g",
  "Aldh1a3","Krt79","Gabrp","Krt6b","Krt10","Sprr1b","Foxi1","Atp6v1g3","Aqp6",
  "Krt5","Krt14","Trp63","Serpina3n","Fst","Ndufa4l2","Col17a1","Ovol1","Dsg3",
  "Mki67","Top2a","Ccnb1","Aspm","Cdh6","Vcam1","Igfbp2","Bmp6","Col5a2","Vcan","Dpt",
  "Twist1","Sox2","Mmp13","Lrrc15","Aspn","Col11a1","Ar","Pbsn","GFP","mCherry"), rownames(so))
p3 <- DotPlot(so, features=feats, group.by="substate") + RotatedAxis() + theme(axis.text.x=element_text(size=7)) + ggtitle("Harmony res0.8")
ggsave(file.path(OUT,"edited_dotplot_substate.png"), p3, width=16, height=6.5, dpi=150)
p4 <- DotPlot(so, features=feats, group.by="substate_scvi") + RotatedAxis() + theme(axis.text.x=element_text(size=7)) + ggtitle("scVI res0.8")
ggsave(file.path(OUT,"edited_dotplot_substate_scvi.png"), p4, width=16, height=6.5, dpi=150)

# ---- 5. Concordance: Harmony (rows) vs scVI-native (cols) ----
ct  <- table(so$substate, so$substate_scvi)
ctn <- sweep(ct, 1, pmax(rowSums(ct),1), "/")
png(file.path(OUT,"edited_harmony_vs_scvi_labels_heatmap.png"), width=1700, height=1300, res=150)
pheatmap(ctn, cluster_rows=FALSE, cluster_cols=FALSE,
         display_numbers=matrix(as.integer(ct), nrow=nrow(ct)),
         number_format="%d", fontsize_number=7,
         main="Harmony sub-state (rows) vs scVI-native sub-state (cols) - row fraction",
         color=colorRampPalette(c("white","#4575b4"))(50))
dev.off()
ARI <- adjustedRandIndex(as.character(so$substate), as.character(so$substate_scvi))
cat("\nARI(Harmony sub-state, scVI-native sub-state) =", round(ARI,3), "\n")

# per-Harmony-state agreement (fraction of cells whose scVI label matches)
agree <- sapply(levels(so$substate), function(x) {
  ix <- so$substate == x
  if (!any(ix)) return(NA_real_)
  mean(as.character(so$substate_scvi[ix]) == x)
})
cat("\n=== per-Harmony-state agreement (fraction with matching scVI label) ===\n")
print(round(agree, 3))
write.csv(data.frame(harmony_state=names(agree), n=as.integer(table(so$substate)[names(agree)]), scvi_match_frac=round(agree,3)),
          file.path(OUT,"harmony_vs_scvi_per_state_agreement.csv"), row.names=FALSE)

# save
saveRDS(so, file.path(OUT,"npp53_edited.RDS"))
write.csv(data.frame(cell=colnames(so), harmony_substate=so$substate, scvi_substate=so$substate_scvi,
                     scvi_cluster=so$scvi_cl, run=so$run, sample_id=so$sample_id, treatment=so$treatment),
          file.path(OUT,"npp53_edited_labels_dual.csv"), row.names=FALSE)
cat("DONE relabel + scVI-native annotation\n")

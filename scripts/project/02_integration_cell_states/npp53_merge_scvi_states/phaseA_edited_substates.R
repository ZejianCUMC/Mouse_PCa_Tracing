# Phase A: build the edited (GFP+ & mCherry+, host-removed) tumour object across ALL samples,
# re-embed, cluster at finer resolution, and score candidate tumour sub-states for annotation.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr); library(harmony) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_clean_labeled.RDS")); DefaultAssay(so) <- "RNA"
ed  <- read.csv(file.path(OUT, "npp53_edited_cells_preview.csv"))   # 20,978 edited cells (all runs)
so  <- so[, colnames(so) %in% ed$cell]
cat("edited cells:", ncol(so), "\n"); print(table(so$run))

# ---- re-embed tumour-only cells ----
DefaultAssay(so) <- "RNA"
so <- NormalizeData(so) |> FindVariableFeatures(nfeatures=2000) |> ScaleData() |> RunPCA(npcs=50, verbose=FALSE)
so <- RunHarmony(so, group.by.vars="batch", reduction.use="pca")
so <- RunUMAP(so, reduction="harmony", dims=1:30, reduction.name="umap.ed", verbose=FALSE)
so <- FindNeighbors(so, reduction="harmony", dims=1:30, verbose=FALSE)
so <- FindClusters(so, resolution=0.8, verbose=FALSE); so$sub08 <- so$seurat_clusters
so <- FindClusters(so, resolution=1.0, verbose=FALSE); so$sub10 <- so$seurat_clusters
Idents(so) <- "sub08"

# ---- candidate sub-state signatures (mouse prostate tumour) ----
sig <- list(
  Luminal_ARhigh = c("Ar","Nkx3-1","Pbsn","Fkbp5","Tmprss2","Krt8","Krt18"),
  Luminal_ARlow  = c("Krt8","Krt18","Krt4","Tacstd2","Psca"),
  Basal          = c("Krt5","Krt14","Trp63","Krt15"),
  Proliferating  = c("Mki67","Top2a","Ccnb1","Cdk1","Birc5"),
  PartialEMT     = c("Vim","Cd44","Krt8","Vcan","Sparc"),
  Mesenchymal    = c("Col3a1","Col1a1","Zeb1","Zeb2","Fn1","Igfbp2"),
  Hillock_Club   = c("Krt4","Krt13","Ltf","Scgb1a1","Aqp3"),
  Neuroendocrine = c("Chga","Chgb","Syp","Ascl1","Ncam1"),
  Interferon     = c("Isg15","Ifit1","Ifit3","Irf7","Stat1") )
sig <- lapply(sig, function(g) intersect(g, rownames(so))); sig <- sig[sapply(sig,length)>=2]
d <- SeuratObject::LayerData(so, assay="RNA", layer="data")
zrow <- function(m){ m<-as.matrix(m); mu<-rowMeans(m); s<-apply(m,1,sd); s[s==0]<-1; (m-mu)/s }
scores <- sapply(names(sig), function(nm) colMeans(zrow(d[sig[[nm]], ,drop=FALSE])))
for (nm in colnames(scores)) so@meta.data[[paste0("sig_",nm)]] <- scores[,nm]
# per-cluster mean signature -> argmax hint
cl <- so$sub08; mat <- apply(scores, 2, function(x) tapply(x, cl, mean))
hint <- setNames(colnames(mat)[max.col(mat, ties.method="first")], rownames(mat))
cat("\n=== sub08 cluster -> signature argmax hint ===\n")
print(data.frame(cluster=rownames(mat), hint=hint, n=as.integer(table(cl))))
cat("\n=== per-cluster mean signature scores ===\n"); print(round(mat,2))

# ---- top DE markers per cluster ----
mk <- FindAllMarkers(so, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.4, verbose=FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n=12) %>% summarise(genes=paste(gene,collapse=", "))
cat("\n=== top markers per sub08 cluster ===\n"); print(as.data.frame(top), right=FALSE)
write.csv(mk, file.path(OUT,"edited_sub08_markers.csv"), row.names=FALSE)

# ---- plots ----
save1 <- function(g,ttl,fn,w=8.6,h=6.2){ p<-DimPlot(so,reduction="umap.ed",group.by=g,label=(g=="sub08"),raster=FALSE)+ggtitle(ttl); ggsave(file.path(OUT,fn),p,width=w,height=h,dpi=150) }
save1("sub08","Edited tumour — clusters (res 0.8)","edited_umap_sub08.png")
save1("celltype","Edited tumour — coarse cell type","edited_umap_celltype.png")
save1("treatment","Edited tumour — treatment","edited_umap_treatment.png")
save1("run","Edited tumour — run/batch","edited_umap_run.png")
mkfeat <- intersect(c("Ar","Nkx3-1","Pbsn","Krt8","Krt5","Krt14","Trp63","Mki67","Top2a","Vim","Cd44",
                      "Col3a1","Zeb1","Fn1","Krt4","Chga","Syp","Isg15","GFP","mCherry"), rownames(so))
p <- DotPlot(so, features=mkfeat, group.by="sub08") + RotatedAxis() + ggtitle("Edited tumour sub08 — markers")
ggsave(file.path(OUT,"edited_dotplot_sub08.png"), p, width=13, height=5.5, dpi=150)

saveRDS(so, file.path(OUT,"npp53_edited.RDS"))
# per-cell raw cellBC (strip run prefix + -1 suffix) for Phase B barcode join
raw <- sub("-\\d+$","", sub("^MJZ\\d+_","", colnames(so)))
write.csv(data.frame(cell=colnames(so), rawBC=raw, run=so$run, sample_id=so$sample_id,
                     treatment=so$treatment, celltype=so$celltype, sub08=so$sub08, sub10=so$sub10),
          file.path(OUT,"npp53_edited_meta.csv"), row.names=FALSE)
cat("\nDONE phaseA edited substates\n")

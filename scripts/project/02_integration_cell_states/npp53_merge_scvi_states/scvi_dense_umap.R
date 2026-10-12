suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_edited.RDS")); DefaultAssay(so) <- "RNA"
lev <- levels(so$substate)
pal <- setNames(c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#fdae6b","#6a51a3",
  "#238b45","#00441b","#41ab5d","#a1d99b","#969696","#525252",
  "#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"), lev)
# denser scVI UMAP: min.dist 0.3 -> 0.1, n.neighbors 30 -> 50
so <- RunUMAP(so, reduction="scvi", dims=1:ncol(Embeddings(so,"scvi")),
              min.dist=0.1, n.neighbors=50, reduction.name="umap.scvi", verbose=FALSE)
dim2 <- function(red,ttl,fn){ p<-DimPlot(so,reduction=red,group.by="substate",cols=pal,raster=FALSE,label=TRUE,label.size=3,repel=TRUE)+ggtitle(ttl); ggsave(file.path(OUT,fn),p,width=10.5,height=6.6,dpi=150)
  p2<-DimPlot(so,reduction=red,group.by="substate",cols=pal,raster=FALSE)+ggtitle(ttl); ggsave(file.path(OUT,sub("\\.png$","_nolabel.png",fn)),p2,width=10.5,height=6.6,dpi=150) }
dim2("umap.scvi","scVI UMAP","edited_umap_substate_scvi.png")
saveRDS(so, file.path(OUT,"npp53_edited.RDS"))
cat("DONE dense scVI umap (min.dist=0.1, n.neighbors=50)\n")

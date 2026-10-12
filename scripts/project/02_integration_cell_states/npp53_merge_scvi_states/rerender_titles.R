suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so <- readRDS(file.path(OUT,"npp53_edited.RDS"))
lev <- levels(so$substate)
pal <- setNames(c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#fdae6b","#6a51a3",
                  "#238b45","#00441b","#41ab5d","#a1d99b","#969696","#525252",
                  "#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"), lev)
save1 <- function(red, grp, ttl, fn, pal=NULL) { p <- DimPlot(so, reduction=red, group.by=grp, cols=pal, raster=FALSE, label=TRUE, label.size=3, repel=TRUE) + ggtitle(ttl); ggsave(file.path(OUT, fn), p, width=10.5, height=6.6, dpi=150) }
save1("umap.ed",   "substate",      "Harmony res0.8", "edited_umap_substate.png", pal)
save1("umap.scvi", "substate",      "scVI UMAP",      "edited_umap_scvi_HarmLabels.png", pal)
lev_s <- levels(so$substate_scvi); pal_s <- pal[lev_s]
save1("umap.scvi", "substate_scvi", "scVI res0.8",    "edited_umap_scvi_scviLabels.png", pal_s)
cat("DONE re-render titles\n")

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_edited.RDS"))
lev <- levels(so$substate)
pal <- setNames(c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#fdae6b","#6a51a3",
  "#238b45","#00441b","#41ab5d","#a1d99b","#969696","#525252",
  "#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"), lev)
TT <- "Harmony res0.8"
p <- DimPlot(so, reduction="umap.ed", group.by="substate", cols=pal, raster=FALSE, label=TRUE, label.size=3, repel=TRUE) + ggtitle(TT)
ggsave(file.path(OUT,"edited_umap_substate.png"), p, width=10.5, height=6.6, dpi=150)
p2 <- DimPlot(so, reduction="umap.ed", group.by="substate", cols=pal, raster=FALSE) + ggtitle(TT)
ggsave(file.path(OUT,"edited_umap_substate_nolabel.png"), p2, width=10.5, height=6.6, dpi=150)
feats <- intersect(c("Pigr","Ppp1r1b","Rbp2","Cst6","Defb42","Eppin","Agr2","Reg3g",
  "Aldh1a3","Krt79","Gabrp","Krt6b","Krt10","Sprr1b","Foxi1","Atp6v1g3","Aqp6",
  "Krt5","Krt14","Trp63","Serpina3n","Fst","Ndufa4l2","Col17a1","Ovol1","Dsg3",
  "Mki67","Top2a","Ccnb1","Aspm","Cdh6","Vcam1","Igfbp2","Bmp6","Col5a2","Vcan","Dpt",
  "Twist1","Sox2","Mmp13","Lrrc15","Aspn","Col11a1","Ar","Pbsn","GFP","mCherry"), rownames(so))
p3 <- DotPlot(so, features=feats, group.by="substate") + RotatedAxis() +
      theme(axis.text.x=element_text(size=7)) + ggtitle(TT)
ggsave(file.path(OUT,"edited_dotplot_substate.png"), p3, width=16, height=6.5, dpi=150)
cat("DONE rerender\n")

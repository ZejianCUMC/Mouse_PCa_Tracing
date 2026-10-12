suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_edited.RDS")); DefaultAssay(so) <- "RNA"

# sub08 cluster -> detailed sub-state (markers + signature based)
map <- c(
 "1"="Luminal-Secretory","2"="Luminal-Club","7"="Luminal-Distal","17"="Luminal-Agr2",
 "5"="Club-Hillock","10"="Squamous","15"="Ionocyte-like",
 "0"="Basal-Inflammatory","4"="Basal-Hypoxic","6"="Basal","9"="Basal-Differentiating",
 "3"="Cycling-S","14"="Cycling-G2M",
 "8"="Partial-EMT","13"="Mes-Igfbp2","12"="Mes-Collagen","11"="Mes-EMT-TF","16"="Mes-CAF-like")
lev <- c("Luminal-Secretory","Luminal-Club","Luminal-Distal","Luminal-Agr2","Club-Hillock","Squamous",
         "Ionocyte-like","Basal","Basal-Inflammatory","Basal-Hypoxic","Basal-Differentiating",
         "Cycling-S","Cycling-G2M","Partial-EMT","Mes-Igfbp2","Mes-Collagen","Mes-EMT-TF","Mes-CAF-like")
so$substate <- factor(unname(map[as.character(so$sub08)]), levels=lev)
pal <- setNames(c(
 "#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#fdae6b",   # luminal/club/squamous
 "#6a51a3",                                                     # ionocyte
 "#238b45","#00441b","#41ab5d","#a1d99b",                        # basal
 "#969696","#525252",                                           # cycling
 "#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"), lev)        # pEMT->mes
cat("=== substate counts ===\n"); print(table(so$substate))

p <- DimPlot(so, reduction="umap.ed", group.by="substate", cols=pal, raster=FALSE, label=TRUE, label.size=3, repel=TRUE) +
     ggtitle("Edited tumour — detailed sub-states")
ggsave(file.path(OUT,"edited_umap_substate.png"), p, width=10.5, height=6.6, dpi=150)
p2 <- DimPlot(so, reduction="umap.ed", group.by="substate", cols=pal, raster=FALSE) + ggtitle("Edited tumour — detailed sub-states")
ggsave(file.path(OUT,"edited_umap_substate_nolabel.png"), p2, width=10.5, height=6.6, dpi=150)

feats <- intersect(c("Pigr","Ppp1r1b","Rbp2","Cst6","Defb42","Eppin","Agr2","Reg3g",
  "Aldh1a3","Krt79","Gabrp","Krt6b","Krt10","Sprr1b","Foxi1","Atp6v1g3","Aqp6",
  "Krt5","Krt14","Trp63","Serpina3n","Fst","Ndufa4l2","Col17a1","Ovol1","Dsg3",
  "Mki67","Top2a","Ccnb1","Aspm","Cdh6","Vcam1","Igfbp2","Bmp6","Col5a2","Vcan","Dpt",
  "Twist1","Sox2","Mmp13","Lrrc15","Aspn","Col11a1","Ar","Pbsn","GFP","mCherry"), rownames(so))
p3 <- DotPlot(so, features=feats, group.by="substate") + RotatedAxis() +
      theme(axis.text.x=element_text(size=7)) + ggtitle("Detailed sub-state markers")
ggsave(file.path(OUT,"edited_dotplot_substate.png"), p3, width=16, height=6.5, dpi=150)

saveRDS(so, file.path(OUT,"npp53_edited.RDS"))
write.csv(data.frame(cell=colnames(so), rawBC=sub("-\\d+$","",sub("^MJZ\\d+_","",colnames(so))),
                     run=so$run, sample_id=so$sample_id, treatment=so$treatment,
                     celltype=so$celltype, sub08=so$sub08, substate=so$substate),
          file.path(OUT,"npp53_edited_meta.csv"), row.names=FALSE)
cat("DONE detailed annotation\n")

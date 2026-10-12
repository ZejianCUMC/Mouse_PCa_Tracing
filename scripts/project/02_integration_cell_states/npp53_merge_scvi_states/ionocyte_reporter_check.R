suppressPackageStartupMessages({ library(Seurat) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so <- readRDS(file.path(OUT,"npp53_edited.RDS")); DefaultAssay(so) <- "RNA"
cnt <- SeuratObject::LayerData(so,assay="RNA",layer="counts")
dat <- SeuratObject::LayerData(so,assay="RNA",layer="data")
g <- so$substate
tab <- data.frame(substate=levels(g), n=as.integer(table(g)))
tab$pctGFP     <- round(tapply(cnt["GFP",]>0, g, mean)[tab$substate],3)
tab$pctmCherry <- round(tapply(cnt["mCherry",]>0, g, mean)[tab$substate],3)
tab$meanGFP     <- round(tapply(dat["GFP",], g, mean)[tab$substate],2)
tab$meanmCherry <- round(tapply(dat["mCherry",], g, mean)[tab$substate],2)
tab$medGFPcount <- tapply(cnt["GFP",], g, median)[tab$substate]
tab$medmChcount <- tapply(cnt["mCherry",], g, median)[tab$substate]
print(tab, row.names=FALSE)
cat("\nOverall edited-set mean GFP:", round(mean(dat["GFP",]),2),
    " mean mCherry:", round(mean(dat["mCherry",]),2), "\n")

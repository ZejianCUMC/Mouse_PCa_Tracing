# Promote substate_scvi to the primary cell-state annotation and re-render Harmony UMAP with the scVI-native labels.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT,"npp53_edited.RDS"))
# Add a stable "primary" column pointing to scVI-native
so$substate_primary <- so$substate_scvi
lev_s <- levels(so$substate_scvi)
pal <- c("Luminal-Secretory"="#08519c","Luminal-Club"="#3182bd","Luminal-Distal"="#6baed6",
         "Luminal-Agr2"="#9ecae1","Luminal-Club-Hillock"="#c6dbef","Luminal-Squamous"="#fdae6b",
         "Ionocyte"="#6a51a3","Basal"="#238b45","Basal-Inflammatory"="#00441b",
         "Basal-Hypoxic"="#41ab5d","Basal-Differentiating"="#a1d99b",
         "Cycling-S"="#969696","Cycling-G2M"="#525252",
         "Partial-EMT"="#fb6a4a","Mes-Igfbp2"="#ef3b2c","Mes-Collagen"="#cb181d",
         "Mes-EMT-TF"="#a50f15","Mes-CAF-like"="#67000d")
pal_s <- pal[lev_s]
# also render Harmony UMAP colored by scVI-native labels for a full 4-panel view
p <- DimPlot(so, reduction="umap.ed",   group.by="substate_scvi", cols=pal_s, raster=FALSE, label=TRUE, label.size=3, repel=TRUE) + ggtitle("Harmony UMAP")
ggsave(file.path(OUT,"edited_umap_harmony_scviLabels.png"), p, width=10.5, height=6.6, dpi=150)
saveRDS(so, file.path(OUT,"npp53_edited.RDS"))
# refresh dual-label csv with primary column
write.csv(data.frame(cell=colnames(so), primary=so$substate_primary,
                     substate_scvi=so$substate_scvi, substate_harmony=so$substate,
                     scvi_cluster=so$scvi_cl, run=so$run,
                     sample_id=so$sample_id, treatment=so$treatment),
          file.path(OUT,"npp53_edited_labels_dual.csv"), row.names=FALSE)
cat("=== primary now = substate_scvi (16 states) ===\n")
print(table(so$substate_primary))
cat("DONE lock scVI-native as primary\n")

suppressPackageStartupMessages({ library(Seurat) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so <- readRDS(file.path(OUT, "npp53_clean_labeled.RDS")); DefaultAssay(so) <- "RNA"
asn <- read.csv(file.path(OUT, "mes_subcluster_assign.csv"))
host_h <- asn$cell[asn$sub_harmony %in% c(1,3,7,8)]   # Harmony res0.5 host-fibroblast subclusters
cat("Harmony host-fib subclusters 1,3,7,8:", length(host_h), "cells\n")
rm_main <- colnames(so)[so$seurat_clusters %in% c("9","10","14")]  # Macro+Macro+Endo
cat("main clusters 9,10,14:", length(rm_main), "cells\n")
keep0 <- setdiff(colnames(so), union(host_h, rm_main))
cat("after removing host+immune+endo:", length(keep0), "cells\n")
cnt <- SeuratObject::LayerData(so, assay="RNA", layer="counts")
gfp <- cnt["GFP", ] > 0; mch <- cnt["mCherry", ] > 0
dpos <- colnames(so)[gfp & mch]
final <- intersect(keep0, dpos)
cat("\nGFP+ & mCherry+ (double-positive) among kept:", length(final), "\n")
sub <- so[, final]
cat("\n=== final edited cells by celltype ===\n"); print(table(sub$celltype))
cat("\n=== by sample_id ===\n"); print(table(sub$sample_id))
cat("\n=== by run ===\n"); print(table(sub$run))
cat("\n=== by treatment ===\n"); print(table(sub$treatment))
write.csv(data.frame(cell=final, sample_id=sub$sample_id, run=sub$run,
                     treatment=sub$treatment, celltype=sub$celltype),
          file.path(OUT, "npp53_edited_cells_preview.csv"), row.names=FALSE)
cat("DONE count\n")

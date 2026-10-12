suppressPackageStartupMessages({ library(Seurat) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so <- readRDS(file.path(OUT, "npp53_reporter_pos.RDS"))
asn <- read.csv(file.path(OUT, "mes_subcluster_assign.csv"))
host <- asn$cell[asn$sub_scvi %in% c(1,3)]               # scVI reporter- host fibroblasts
cat("reporter_pos in:", ncol(so), "| host fibroblasts (scVI sub 1+3):", length(host), "\n")
tumor <- so[, !(colnames(so) %in% host)]
cat("refined tumor-tracing cells:", ncol(tumor), "\n"); print(table(tumor$celltype))
saveRDS(tumor, file.path(OUT, "npp53_tumor_tracing.RDS"))
write.csv(data.frame(cell=colnames(tumor), celltype=tumor$celltype, cluster=tumor$seurat_clusters,
                     sample_id=tumor$sample_id, treatment=tumor$treatment, run=tumor$run),
          file.path(OUT, "npp53_tumor_tracing_cells.csv"), row.names=FALSE)
write.csv(data.frame(cell=host), file.path(OUT, "npp53_host_fibroblast_cells.csv"), row.names=FALSE)
cat("DONE refine ->", file.path(OUT,"npp53_tumor_tracing.RDS"), "\n")

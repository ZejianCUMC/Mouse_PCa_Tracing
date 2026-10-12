suppressPackageStartupMessages({ library(Seurat); library(loupeR) })
OUT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_merge"
so  <- readRDS(file.path(OUT, "npp53_clean_labeled.RDS")); DefaultAssay(so) <- "RNA"
counts <- SeuratObject::LayerData(so, assay="RNA", layer="counts")
bc <- colnames(counts)
sv <- read.csv(file.path(OUT, "scvi_umap_clean.csv"), row.names=1)
stopifnot(all(bc %in% rownames(sv)))
getE <- function(r) { e <- Embeddings(so, r)[bc, 1:2, drop=FALSE]; colnames(e) <- c(paste0(r,"_1"), paste0(r,"_2")); e }
projections <- list(
  UMAP_Uncorrected = getE("umap.uncorr"),
  UMAP_Harmony     = getE("umap.harmony"),
  UMAP_fastMNN     = getE("umap.mnn"),
  UMAP_scVI        = as.matrix(sv[bc, 1:2]) )
clusters <- list(
  putative_cell_type = factor(so$celltype),
  sample_id          = factor(as.character(so$sample_id)),
  run_id             = factor(as.character(so$run)),
  treatment          = factor(as.character(so$treatment)) )
cat("cells:", length(bc), "| projections:", paste(names(projections),collapse=","), "\n")
print(sapply(clusters, nlevels))
create_loupe(counts, clusters=clusters, projections=projections,
             output_dir=OUT, output_name="npp53_hq_49619cells", force=TRUE)
cat("DONE cloupe ->", file.path(OUT,"npp53_hq_49619cells.cloupe"), "\n")

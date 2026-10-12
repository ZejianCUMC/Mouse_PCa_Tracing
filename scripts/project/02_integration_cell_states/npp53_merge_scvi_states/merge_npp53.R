# Merge Npp53 cohort (24 QC-pass new samples + JZ136) on cellbender counts, then batch-correct with
# Harmony + fastMNN. Exports raw counts (mtx) + meta for the separate scVI job.
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(Matrix); library(harmony); library(SeuratWrappers)
})
ROOT   <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing"
OBJDIR <- file.path(ROOT, "seurat_objs")
OUT    <- file.path(ROOT, "1-mouse_based/Npp53_tracing/data/npp53_merge"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

DROP        <- c("JZ215", "JZ233")                    # <1000 HQ cells
NEW_COHORTS <- c("MJZ019","MJZ021","MJZ022","MJZ023","MJZ024")

join_if <- function(o) if (length(SeuratObject::Layers(o, assay="RNA")) > 1) JoinLayers(o) else o
get_counts <- function(o) SeuratObject::LayerData(o, assay="RNA", layer="counts")   # single joined layer
cts_list <- list(); meta_list <- list()

for (co in NEW_COHORTS) {
  so <- readRDS(file.path(OBJDIR, paste0(co,"_per_sample_DF"), paste0(co,"_seurat_DF_clean.RDS")))
  so <- join_if(so)
  so <- so[, !(so$sample_id %in% DROP)]
  m  <- get_counts(so); colnames(m) <- paste0(co, "_", colnames(m))
  cts_list[[co]] <- m
  meta_list[[co]] <- data.frame(cell=colnames(m), sample_id=so$sample_id, run=co, batch=co,
                                treatment=so$treatment, model="Npp53", stringsAsFactors=FALSE)
  cat(co, ": ", ncol(m), " cells\n", sep="")
}
# JZ136 from MJZ008 cellbender-adjusted (already Singlet + doublet-removed)
m8 <- readRDS(file.path(ROOT, "MJZ008_cellbender_adj.RDS"))
DefaultAssay(m8) <- "RNA"
m8 <- DietSeurat(m8, assays = "RNA")     # drop SCT so JoinLayers/LayerData act on RNA
m8 <- join_if(m8)
j  <- m8[, m8$SampleID == "JZ136"]
mj <- get_counts(j); colnames(mj) <- paste0("MJZ008_", colnames(mj))
cts_list[["JZ136"]] <- mj
meta_list[["JZ136"]] <- data.frame(cell=colnames(mj), sample_id="JZ136", run="MJZ008", batch="MJZ008",
                                   treatment="Intact", model="Npp53", stringsAsFactors=FALSE)
cat("JZ136 : ", ncol(mj), " cells\n", sep="")

# common genes + combine
genes <- Reduce(intersect, lapply(cts_list, rownames))
cat("common genes:", length(genes), "\n")
combined <- do.call(cbind, lapply(cts_list, function(m) m[genes, , drop=FALSE]))
meta <- do.call(rbind, meta_list); rownames(meta) <- meta$cell
stopifnot(all(colnames(combined) == rownames(meta)))
cat("MERGED:", ncol(combined), "cells x", nrow(combined), "genes\n")
print(table(meta$batch))

so <- CreateSeuratObject(counts = combined, meta.data = meta)
so <- NormalizeData(so) |> FindVariableFeatures(nfeatures=2000) |> ScaleData() |> RunPCA(npcs=50, verbose=FALSE)
so <- RunUMAP(so, reduction="pca", dims=1:30, reduction.name="umap.uncorr", verbose=FALSE)
# Harmony
so <- RunHarmony(so, group.by.vars="batch", reduction.use="pca")
so <- RunUMAP(so, reduction="harmony", dims=1:30, reduction.name="umap.harmony", verbose=FALSE)
# fastMNN — RunFastMNN returns a NEW object; graft its 'mnn' embedding back into `so`
mnn_obj <- RunFastMNN(object.list = SplitObject(so, split.by="batch"))
emb <- Embeddings(mnn_obj, "mnn")[colnames(so), , drop=FALSE]
so[["mnn"]] <- CreateDimReducObject(embeddings = emb, key = "mnn_", assay = DefaultAssay(so))
so <- RunUMAP(so, reduction="mnn", dims=1:min(30, ncol(emb)), reduction.name="umap.mnn", verbose=FALSE)

# save + export FIRST (so scVI can proceed even if plotting hiccups)
saveRDS(so, file.path(OUT, "npp53_merged_integrated.RDS"))
writeMM(combined, file.path(OUT, "npp53_counts.mtx"))
writeLines(genes,               file.path(OUT, "npp53_genes.txt"))
writeLines(colnames(combined),  file.path(OUT, "npp53_barcodes.txt"))
write.csv(meta, file.path(OUT, "npp53_meta.csv"), row.names=FALSE)
cat("saved merged RDS + counts export ->", OUT, "\n")

# QC UMAPs (batch) — guarded
suppressPackageStartupMessages(library(ggplot2))
for (r in c("umap.uncorr","umap.harmony","umap.mnn")) {
  tryCatch({
    p <- DimPlot(so, reduction=r, group.by="batch") + ggtitle(r)
    ggsave(file.path(OUT, paste0("UMAP_", r, "_byBatch.png")), p, width=7, height=6, dpi=150)
  }, error=function(e) cat("plot", r, "failed:", conditionMessage(e), "\n"))
}
cat("DONE merge + Harmony + fastMNN; exported counts for scVI ->", OUT, "\n")

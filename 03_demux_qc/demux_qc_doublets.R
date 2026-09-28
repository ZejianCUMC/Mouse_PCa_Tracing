# =============================================================================
# 03 · Hashtag demultiplexing + cell QC + DoubletFinder  (Seurat, R)
# =============================================================================
# Description : One 10x run pools several mice, each labelled with a TotalSeq-B
#               hashtag antibody. This script
#                 1. reads Cell Ranger filtered GEX + HTO counts,
#                 2. cell QC: nFeature_RNA > 200, percent.mt < 10,
#                 3. HTODemux (CLR, positive.quantile 0.85) -> keep Singlets,
#                 4. maps hashtag -> mouse (sample_id) -> treatment,
#                 5. swaps in CellBender ambient-corrected counts,
#                 6. DoubletFinder per mouse (catches same-hashtag doublets),
#                 7. keeps ${KEEP_SAMPLES} and exports counts for Python.
#               Same recipe as the Npp53 cohort (htodemux template + doubletfinder_per_sample.R).
#               Newer runs with high ambient HTO use CellBender-HTO + demuxmix instead
#               (see HANDSOFF.md, "Going further").
# Input       : ${OUT_CR}/${RUN}/outs/filtered_feature_bc_matrix/,
#               ${OUT_CB}/${RUN}_cellbender_filtered.h5, ${WORK}/metadata/hash_map.csv
# Output      : ${OUT_QC}/  {RUN}_qc_clean.RDS, counts.mtx, genes.txt, barcodes.txt,
#               cell_meta.csv, demux/QC tables + plots
# Conda env   : HPC ${ENVS}/ST  (Seurat 5, DoubletFinder, DropletUtils)
# Validated   : 2026-09-27 (logic from the Npp53 cohort scripts, 2026-08)
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat); library(DoubletFinder); library(DropletUtils); library(dplyr)
  library(Matrix); library(ggplot2)
})
options(expressions = 500000)

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
DEFAULTS <- list(
  min_features   = 200,     # cell QC: minimum detected genes
  max_pct_mt     = 10,      # cell QC: maximum % mitochondrial UMIs
  mt_pattern     = "^mt-",  # mouse mitochondrial genes
  hto_quantile   = 0.85,    # HTODemux positive.quantile
  df_pcs         = 40,      # DoubletFinder: PCs
  df_resolution  = 0.5,     # DoubletFinder: clustering resolution for homotypic model
  df_rate_per_cell = 8e-6,  # 10x multiplet rate ~0.8% per 1,000 cells
  df_pK_fallback = 0.09,    # used when paramSweep fails
  seed           = as.integer(Sys.getenv("SEED", "42"))
)

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
env <- function(k) { v <- Sys.getenv(k); if (v == "") stop("missing env var ", k, " (source config/config.sh)"); v }

read_cellbender <- function(h5) {
  # CellBender h5 is not the 10x schema -> read with DropletUtils, keep GEX rows, use symbols
  sce <- DropletUtils::read10xCounts(h5)
  m   <- as(SummarizedExperiment::assay(sce, "counts"), "CsparseMatrix")
  rd  <- SummarizedExperiment::rowData(sce)
  keep <- !is.na(rd$Symbol) & !duplicated(rd$Symbol) & rd$Symbol != ""
  if ("Type" %in% colnames(rd)) keep <- keep & (rd$Type == "Gene Expression" | is.na(rd$Type))
  m <- m[keep, ]; rownames(m) <- rd$Symbol[keep]
  colnames(m) <- SummarizedExperiment::colData(sce)$Barcode
  m
}

run_doubletfinder <- function(so, p) {
  so <- SCTransform(so, verbose = FALSE)
  so <- RunPCA(so, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:p$df_pcs, verbose = FALSE)
  so <- FindClusters(so, resolution = p$df_resolution, verbose = FALSE)
  homotypic <- modelHomotypic(so$seurat_clusters)
  nExp      <- round(ncol(so) * p$df_rate_per_cell * ncol(so))
  nExp_adj  <- round(nExp * (1 - homotypic))
  pK <- p$df_pK_fallback
  ps <- tryCatch({
    sw <- paramSweep(so, PCs = 1:p$df_pcs, sct = TRUE)
    find.pK(summarizeSweep(sw, GT = FALSE))
  }, error = function(e) { message("  paramSweep failed -> pK = ", pK); NULL })
  if (!is.null(ps) && any(!is.na(ps$BCmetric))) pK <- as.numeric(as.character(ps$pK[which.max(ps$BCmetric)]))
  message("  pK = ", pK, "  nExp = ", nExp, "  nExp.adj = ", nExp_adj)
  so <- doubletFinder(so, PCs = 1:p$df_pcs, pN = 0.25, pK = pK, nExp = nExp,     sct = TRUE)
  so <- doubletFinder(so, PCs = 1:p$df_pcs, pN = 0.25, pK = pK, nExp = nExp_adj, sct = TRUE)
  cls <- grep("^DF.classifications", colnames(so@meta.data), value = TRUE)
  ifelse(so@meta.data[[cls[1]]] == "Singlet" & so@meta.data[[cls[2]]] == "Singlet", "Singlet",
         ifelse(so@meta.data[[cls[1]]] == "Doublet" & so@meta.data[[cls[2]]] == "Doublet",
                "Doublet_high", "Doublet_low"))
}

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
run_pipeline <- function(input_path, output_dir, ...) {
  p <- modifyList(DEFAULTS, list(...))
  set.seed(p$seed)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  RUN  <- env("RUN")
  keep_samples <- strsplit(env("KEEP_SAMPLES"), ",")[[1]]

  # 1. Cell Ranger filtered matrix (list: Gene Expression + Antibody Capture)
  m  <- Read10X(file.path(input_path, "outs/filtered_feature_bc_matrix"))
  so <- CreateSeuratObject(counts = m$`Gene Expression`)
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = p$mt_pattern)
  n0 <- ncol(so)
  # 2. QC
  p_qc <- VlnPlot(so, c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)
  ggsave(file.path(output_dir, "QC_violin_before_filter.png"), p_qc, width = 10, height = 4, dpi = 150)
  so <- subset(so, nFeature_RNA > p$min_features & percent.mt < p$max_pct_mt)
  n1 <- ncol(so)

  # 3. HTODemux
  hto   <- m$`Antibody Capture`[, colnames(so)]
  so[["HTO"]] <- CreateAssayObject(counts = hto)
  so <- NormalizeData(so, assay = "HTO", normalization.method = "CLR")
  so <- HTODemux(so, assay = "HTO", positive.quantile = p$hto_quantile)
  write.csv(as.data.frame(table(global = so$HTO_classification.global, hash = so$hash.ID)),
            file.path(output_dir, "HTODemux_calls.csv"), row.names = FALSE)
  Idents(so) <- "hash.ID"
  p_r <- RidgePlot(so, assay = "HTO", features = rownames(so[["HTO"]]), ncol = 3)
  ggsave(file.path(output_dir, "HTO_ridgeplot.png"), p_r, width = 12, height = 8, dpi = 150)
  so <- subset(so, HTO_classification.global == "Singlet")
  n2 <- ncol(so)

  # 4. hashtag -> mouse -> treatment
  hm <- read.csv(file.path(env("WORK"), "metadata/hash_map.csv"), stringsAsFactors = FALSE)
  idx <- match(as.character(so$hash.ID), hm$hash_id)
  so$sample_id <- hm$sample_id[idx]; so$treatment <- hm$treatment[idx]; so$model <- hm$model[idx]
  so$run <- RUN
  so <- subset(so, sample_id %in% keep_samples)
  n3 <- ncol(so)

  # 5. CellBender-corrected counts for the kept cells
  cb    <- read_cellbender(file.path(env("OUT_CB"), paste0(RUN, "_cellbender_filtered.h5")))
  joint <- intersect(colnames(so), colnames(cb))
  meta  <- so@meta.data[joint, c("sample_id", "treatment", "model", "run", "hash.ID", "percent.mt")]
  so_cb <- CreateSeuratObject(counts = cb[, joint], meta.data = meta)
  n4 <- ncol(so_cb)

  # 6. DoubletFinder per mouse
  so_cb$DF <- NA_character_
  for (s in unique(so_cb$sample_id)) {
    message("DoubletFinder: ", s)
    cells <- colnames(so_cb)[so_cb$sample_id == s]
    so_cb$DF[cells] <- run_doubletfinder(subset(so_cb, cells = cells), p)
  }
  write.csv(as.data.frame(table(sample = so_cb$sample_id, DF = so_cb$DF)),
            file.path(output_dir, "DoubletFinder_calls.csv"), row.names = FALSE)
  so_cb <- subset(so_cb, DF == "Singlet")

  # 7. save + export for python (genes x cells mtx, keeps GFP / mCherry rows)
  saveRDS(so_cb, file.path(output_dir, paste0(RUN, "_qc_clean.RDS")))
  cts <- SeuratObject::LayerData(so_cb, assay = "RNA", layer = "counts")
  writeMM(cts, file.path(output_dir, "counts.mtx"))
  writeLines(rownames(cts), file.path(output_dir, "genes.txt"))
  writeLines(colnames(cts), file.path(output_dir, "barcodes.txt"))
  md <- so_cb@meta.data; md$cell <- colnames(so_cb)
  write.csv(md, file.path(output_dir, "cell_meta.csv"), row.names = FALSE)

  funnel <- data.frame(step = c("CellRanger cells", "QC pass", "HTO singlet", "kept mice",
                                "in CellBender", "DoubletFinder singlet"),
                       cells = c(n0, n1, n2, n3, n4, ncol(so_cb)))
  write.csv(funnel, file.path(output_dir, "cell_funnel.csv"), row.names = FALSE)
  print(funnel); print(table(so_cb$sample_id, so_cb$treatment))
  list(output_path = output_dir, params_used = p, summary = funnel)
}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if (sys.nframe() == 0) {
  invisible(run_pipeline(file.path(env("OUT_CR"), env("RUN")), env("OUT_QC")))
}

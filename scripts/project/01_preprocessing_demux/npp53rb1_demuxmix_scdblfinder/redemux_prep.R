suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
lib <- commandArgs(TRUE)[1]
BASE <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing"
co   <- file.path(BASE, "counts", lib, "outs")
out  <- file.path(BASE, "seurat_objs", "redemux"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
cat("=== ", lib, " ===\n")

# RNA QC from cellranger filtered (identical thresholds to original template)
m  <- Read10X(file.path(co, "filtered_feature_bc_matrix"))
so <- CreateSeuratObject(counts = m$`Gene Expression`)
so <- PercentageFeatureSet(so, pattern="^mt-", col.name="percent.mt")
so <- subset(so, percent.mt < 10 & nFeature_RNA > 200)
cat("post-QC cells:", ncol(so), "\n")

# cellbender ambient-corrected HTO (pre-extracted to CSV by extract_cbHTO.py)
hto_df <- read.csv(file.path(out, paste0(lib,"_cbHTO.csv")), check.names=FALSE, row.names=1)
hto <- t(as.matrix(hto_df))        # features(B30x) x cells(barcodes)
joint <- intersect(colnames(so), colnames(hto))
cat("RNA n cellbender-HTO cells:", length(joint), "\n")
so  <- so[, joint]; hto <- hto[, joint, drop=FALSE]
hto <- hto[rowSums(hto) > 0, , drop=FALSE]
so[["HTO"]] <- CreateAssayObject(counts = hto)
so <- NormalizeData(so, assay="HTO", normalization.method="CLR")

# Method 1: HTODemux (cellbender HTO)
so <- HTODemux(so, assay="HTO", positive.quantile=0.85)
h_glob <- as.character(so$HTO_classification.global); h_id <- as.character(so$hash.ID)
cat("HTODemux(cellbender):\n"); print(table(h_glob))

# Method 3: MULTIseqDemux (autoThresh)
so <- MULTIseqDemux(so, assay="HTO", autoThresh=TRUE)
m_id <- as.character(so$MULTI_ID)
m_glob <- ifelse(m_id %in% c("Doublet","Negative"), m_id, "Singlet")
cat("MULTIseqDemux:\n"); print(table(m_glob))

df <- data.frame(cell=colnames(so), nFeature_RNA=so$nFeature_RNA,
                 htodemux_global=h_glob, htodemux_id=h_id,
                 multiseq_global=m_glob, multiseq_id=m_id, stringsAsFactors=FALSE)
write.csv(df, file.path(out, paste0(lib,"_htodemux_multiseq.csv")), row.names=FALSE)
saveRDS(hto, file.path(out, paste0(lib,"_cbHTO_counts.rds")))
saveRDS(setNames(so$nFeature_RNA, colnames(so)), file.path(out, paste0(lib,"_nFeature.rds")))
cat("DONE prep", lib, "\n")

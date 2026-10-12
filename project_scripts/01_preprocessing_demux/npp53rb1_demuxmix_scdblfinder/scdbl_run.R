suppressPackageStartupMessages({ library(Matrix); library(SingleCellExperiment); library(scDblFinder); library(DropletUtils) })
lib <- commandArgs(TRUE)[1]
BASE <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing"
co  <- file.path(BASE,"counts",lib,"outs","filtered_feature_bc_matrix")
out <- file.path(BASE,"seurat_objs","redemux")
cat("=== scDblFinder", lib, "===\n")
sce <- read10xCounts(co, col.names=TRUE)
ft  <- rowData(sce)$Type
sce <- sce[ft == "Gene Expression", ]
cnt <- counts(sce); rn <- rowData(sce)$Symbol
nFeat <- Matrix::colSums(cnt > 0)
mt    <- grepl("^mt-", rn)
pmt   <- Matrix::colSums(cnt[mt, , drop=FALSE]) / pmax(Matrix::colSums(cnt), 1) * 100
keep  <- pmt < 10 & nFeat > 200
cat("cells pre-QC:", ncol(sce), " post-QC:", sum(keep), "\n")
sce <- sce[, keep]
set.seed(1)
sce <- scDblFinder(sce)
df <- data.frame(cell=colnames(sce), scdbl_class=as.character(sce$scDblFinder.class),
                 scdbl_score=as.numeric(sce$scDblFinder.score), stringsAsFactors=FALSE)
write.csv(df, file.path(out, paste0(lib,"_scdblfinder.csv")), row.names=FALSE)
cat("scDblFinder class:\n"); print(table(df$scdbl_class))
cat("doublet rate:", round(100*mean(df$scdbl_class=="doublet"),1), "%\n")

suppressPackageStartupMessages({ library(demuxmix); library(Matrix) })
lib <- commandArgs(TRUE)[1]
out <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/seurat_objs/redemux"
hto <- as.matrix(readRDS(file.path(out, paste0(lib,"_cbHTO_counts.rds"))))
nf  <- readRDS(file.path(out, paste0(lib,"_nFeature.rds")))
nf  <- as.numeric(nf[colnames(hto)])
cat("=== demuxmix", lib, " HTOs:", nrow(hto), " cells:", ncol(hto), "===\n")
set.seed(1)
dmm <- demuxmix(hto, rna=nf, maxIter=2000)            # regmix model with RNA (nFeature) covariate
cls <- dmmClassify(dmm)                  # cols: HTO, Prob, Type
cls$cell <- rownames(cls)
write.csv(cls, file.path(out, paste0(lib,"_demuxmix.csv")), row.names=FALSE)
cat("Type:\n"); print(table(cls$Type))
cat("singlet HTO:\n"); print(table(cls$HTO[cls$Type=="singlet"]))
cat("DONE demuxmix", lib, "\n")

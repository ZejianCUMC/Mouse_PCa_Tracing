suppressPackageStartupMessages({library(demuxmix); library(Matrix)})
O <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/seurat_objs/redemux"
lib <- "MJZ028"
cat("demuxmix args:\n"); print(args(demuxmix))
hto <- as.matrix(readRDS(file.path(O,paste0(lib,"_cbHTO_counts.rds"))))
nf  <- as.numeric(readRDS(file.path(O,paste0(lib,"_nFeature.rds")))[colnames(hto)])
warns <- NULL
set.seed(1)
dmm <- withCallingHandlers(
  demuxmix(hto, rna=nf, maxIter=2000),
  warning=function(w){ warns[[length(warns)+1]] <<- conditionMessage(w); invokeRestart("muffleWarning") })
cat("\n=== warnings captured (maxIter=2000) ===\n"); if(length(warns)) print(unlist(warns)) else cat("NONE — all models converged\n")
cat("\n=== summary(dmm) ===\n"); print(summary(dmm))
cat("\n=== dmmOverlap ===\n"); print(round(dmmOverlap(dmm),4))
cls <- dmmClassify(dmm); cls$cell <- rownames(cls)
write.csv(cls, file.path(O,paste0(lib,"_demuxmix_refit.csv")), row.names=FALSE)
cat("\n=== new Type table ===\n"); print(table(cls$Type))
cat("=== new singlet-HTO ===\n"); print(table(cls$HTO[cls$Type=="singlet"]))
# compare to original
old <- read.csv(file.path(O,paste0(lib,"_demuxmix.csv")))
m <- merge(data.frame(cell=cls$cell,Type=cls$Type,HTO=cls$HTO),
           data.frame(cell=old$cell,Type0=old$Type,HTO0=old$HTO), by="cell")
cat("\n=== refit vs original agreement ===\n")
cat("Type match:", round(mean(m$Type==m$Type0)*100,2),"%\n")
sg <- m[m$Type=="singlet" & m$Type0=="singlet",]
cat("HTO match among cells singlet in both:", round(mean(sg$HTO==sg$HTO0)*100,2),"% (n=",nrow(sg),")\n")
cat("singlets: original",sum(old$Type=="singlet")," refit",sum(cls$Type=="singlet"),"\n")

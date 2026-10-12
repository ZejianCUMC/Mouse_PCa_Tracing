suppressPackageStartupMessages({library(demuxmix); library(Matrix)})
O <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/seurat_objs/redemux"
for(lib in c("MJZ025","MJZ026","MJZ027","MJZ028")){
  cat("\n=================== ", lib, " ===================\n")
  hto <- as.matrix(readRDS(file.path(O,paste0(lib,"_cbHTO_counts.rds"))))
  nf  <- as.numeric(readRDS(file.path(O,paste0(lib,"_nFeature.rds")))[colnames(hto)])
  cat("orientation:", nrow(hto),"HTOs x",ncol(hto),"cells  | rownames:",paste(rownames(hto),collapse=","),"\n")
  cat("counts integer?", all(hto==round(hto)), " | count range:",paste(range(hto),collapse="-"),"\n")
  cat("rna covariate: any NA?", any(is.na(nf))," range:",paste(range(nf),collapse="-"),"\n")
  set.seed(1); dmm <- demuxmix(hto, rna=nf)
  cat("--- summary(dmm): per-HTO fitted model ---\n"); print(summary(dmm))
  cat("--- dmmOverlap (pos/neg distribution overlap; lower=cleaner) ---\n"); print(round(dmmOverlap(dmm),4))
  # reproducibility vs saved csv
  cls <- dmmClassify(dmm)
  saved <- read.csv(file.path(O,paste0(lib,"_demuxmix.csv")))
  m <- merge(data.frame(cell=rownames(cls),Type=cls$Type,HTO=cls$HTO),
             data.frame(cell=saved$cell,Type2=saved$Type,HTO2=saved$HTO), by="cell")
  cat("reproduces saved csv? Type match:",mean(m$Type==m$Type2)*100,"%  HTO match:",mean(m$HTO==m$HTO2,na.rm=TRUE)*100,"%\n")
}

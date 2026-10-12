suppressPackageStartupMessages({library(Seurat); library(edgeR); library(fgsea)})
Z <- "/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/24PNAS_Zaidi_MSK_PCa_atlas"
FSdir <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/FitnessSignature"
OUT <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/progression_genesets/zaidi_gsea"
cat("[1] load tumor.cells.rds\n"); t0<-Sys.time()
obj <- readRDS(file.path(Z,"GSE264573_msk.integrated.remove.cellcycle.tumor.cells.rds"))
cat("   loaded", ncol(obj), "cells x", nrow(obj), "genes in", round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
md <- obj@meta.data
scol <- if("subtype" %in% colnames(md)) "subtype" else grep("subtype|disease|state",colnames(md),value=T,ignore.case=T)[1]
pcol <- if("patient" %in% colnames(md)) "patient" else grep("patient",colnames(md),value=T,ignore.case=T)[1]
cat("   subtype col=",scol," patient col=",pcol,"\n"); cat("   subtype table:\n"); print(table(md[[scol]]))
sub <- as.character(md[[scol]]); pat <- as.character(md[[pcol]])
keep <- sub %in% c("CSPC","CRPC")                 # EXCLUDE NEPC
cat("[2] kept CSPC+CRPC cells:", sum(keep), " (excluded NEPC:", sum(sub=="NEPC"),")\n")
suppressPackageStartupMessages(library(Matrix))
cnt <- tryCatch(GetAssayData(obj, assay="RNA", slot="counts"), error=function(e) GetAssayData(obj, layer="counts"))
cnt <- cnt[, keep]; patk <- factor(pat[keep]); subk <- sub[keep]
mm <- sparse.model.matrix(~0+patk)                 # cells x patients indicator
colnames(mm) <- levels(patk)
pb <- as.matrix(cnt %*% mm)                        # genes x patients pseudobulk
ncell <- table(patk); psub <- tapply(subk, patk, function(x) names(which.max(table(x))))
use <- names(ncell)[ncell>=50]
pb <- pb[, use, drop=F]; psub <- psub[use]; ncell <- ncell[use]
cat("[3] pseudobulk patients (>=50 cells):\n"); print(data.frame(patient=use, subtype=psub, n_cells=as.integer(ncell)))
cat("   CSPC=",sum(psub=="CSPC")," CRPC=",sum(psub=="CRPC"),"\n")
y <- DGEList(pb); y <- calcNormFactors(y)
lc <- cpm(y, log=TRUE, prior.count=1)
rnk <- rowMeans(lc[,psub=="CRPC",drop=F]) - rowMeans(lc[,psub=="CSPC",drop=F])   # + = up in CRPC
rnk <- rnk[is.finite(rnk)]; rnk <- sort(rnk, decreasing=TRUE)
write.csv(data.frame(gene=names(rnk), rank_CRPCvsCSPC=rnk), file.path(OUT,"rank_CRPC_vs_CSPC.csv"), row.names=F)
FSU <- readLines(file.path(FSdir,"FitnessPosSignature_human.txt")); FSD <- readLines(file.path(FSdir,"FitnessNegSignature_human.txt"))
FSU <- intersect(FSU, names(rnk)); FSD <- intersect(FSD, names(rnk))
cat("[4] FSU in-data:",length(FSU)," FSD in-data:",length(FSD),"\n")
set.seed(1)
res <- fgsea(list(FSU_highfitness=FSU, FSD_lowfitness=FSD), rnk, minSize=10, maxSize=3000, eps=0)
res$leadingEdge <- sapply(res$leadingEdge, function(x) paste(head(x,15),collapse=","))
write.csv(res, file.path(OUT,"fgsea_FSU_FSD_CRPCvsCSPC.csv"), row.names=F)
cat("\n===== GSEA RESULT (ranking: + = up in CRPC[late], - = up in CSPC[early]) =====\n")
print(res[,c("pathway","ES","NES","pval","padj","size")])
cat("\nInterpretation: FSU NES>0 => FSU up in CRPC(late) ; FSD NES<0 => FSD up in CSPC(early)\n")
cat("ALL DONE\n")

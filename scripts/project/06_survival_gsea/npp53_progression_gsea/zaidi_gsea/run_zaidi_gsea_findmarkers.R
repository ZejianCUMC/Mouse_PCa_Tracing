suppressPackageStartupMessages({library(Seurat); library(clusterProfiler)})
Z <- "/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/24PNAS_Zaidi_MSK_PCa_atlas"
FSdir <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/FitnessSignature"
OUT <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/progression_genesets/zaidi_gsea"
cat("[1] load tumor.cells.rds\n"); t0<-Sys.time()
obj <- readRDS(file.path(Z,"GSE264573_msk.integrated.remove.cellcycle.tumor.cells.rds"))
cat("   loaded",ncol(obj),"cells in",round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
sub <- as.character(obj@meta.data$subtype)
obj$grp <- sub
obj <- subset(obj, subset = grp %in% c("CSPC","CRPC"))     # non-NE: exclude NEPC
cat("[2] cells CSPC=",sum(obj$grp=="CSPC")," CRPC=",sum(obj$grp=="CRPC")," (NEPC excluded)\n")
DefaultAssay(obj) <- "RNA"; Idents(obj) <- "grp"
cat("[3] FindMarkers Wilcoxon CRPC vs CSPC (all genes)\n")
dem <- FindMarkers(obj, ident.1="CRPC", ident.2="CSPC", test.use="wilcox",
                   logfc.threshold=0, min.pct=0.01)
dem$gene <- rownames(dem)
write.csv(dem, file.path(OUT,"FindMarkers_CRPC_vs_CSPC.csv"), row.names=FALSE)
cat("   DE genes:",nrow(dem),"\n")
gl <- dem$avg_log2FC; names(gl) <- dem$gene           # + = up in CRPC(late)
gl <- sort(gl, decreasing=TRUE)
FSU <- readLines(file.path(FSdir,"FitnessPosSignature_human.txt"))
FSD <- readLines(file.path(FSdir,"FitnessNegSignature_human.txt"))
t2g <- rbind(data.frame(term="FSU_highfitness", gene=FSU),
             data.frame(term="FSD_lowfitness",  gene=FSD))
t2g <- t2g[t2g$gene %in% names(gl),]
cat("[4] FSU in-data:",sum(t2g$term=="FSU_highfitness")," FSD in-data:",sum(t2g$term=="FSD_lowfitness"),"\n")
set.seed(1)
eg <- GSEA(gl, TERM2GENE=t2g, minGSSize=10, maxGSSize=3000, pvalueCutoff=1, eps=0, seed=TRUE)
res <- as.data.frame(eg)
write.csv(res[,c("ID","setSize","enrichmentScore","NES","pvalue","p.adjust")],
          file.path(OUT,"clusterprofiler_GSEA_FSU_FSD.csv"), row.names=FALSE)
saveRDS(eg, file.path(OUT,"gsea_egobj.rds"))
cat("\n===== clusterProfiler GSEA (rank=avg_log2FC; + up in CRPC[late], - up in CSPC[early]) =====\n")
print(res[,c("ID","setSize","NES","pvalue","p.adjust")])
for(id in res$ID) cat(id,"leadingEdge:",substr(res$core_enrichment[res$ID==id],1,180),"\n")
cat("\nInterpretation: FSU NES>0 => FSU up in CRPC(late); FSD NES<0 => FSD up in CSPC(early)\nALL DONE\n")

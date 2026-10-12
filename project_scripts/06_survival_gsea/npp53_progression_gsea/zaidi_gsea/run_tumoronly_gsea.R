suppressPackageStartupMessages({library(Seurat); library(clusterProfiler); library(fgsea)})
Z<-"/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/24PNAS_Zaidi_MSK_PCa_atlas"
FSdir<-"/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/FitnessSignature"
OUT<-getwd()
cat("[1] load\n"); obj<-readRDS(file.path(Z,"GSE264573_msk.integrated.remove.cellcycle.tumor.cells.rds"))
obj$grp<-as.character(obj@meta.data$subtype)
obj<-subset(obj, subset=grp %in% c("CSPC","CRPC")); DefaultAssay(obj)<-"RNA"
cat("   CSPC=",sum(obj$grp=="CSPC")," CRPC=",sum(obj$grp=="CRPC"),"\n")
## ---- lineage panels ----
pan<-list(Epi=c("EPCAM","KRT8","KRT18","KRT19","KLK3","AR"),
          Imm=c("PTPRC","C1QA","C1QB","LYZ","CD3D","CD3E","NKG7","MS4A1","CD68","CD14","AIF1"),
          Fib=c("COL1A1","DCN","LUM","PDGFRB"), Endo=c("PECAM1","VWF","CLDN5"), RBC=c("HBB","HBA1","HBA2"))
pan<-lapply(pan,function(g) g[g %in% rownames(obj)])
sc<-sapply(names(pan),function(k){ e<-FetchData(obj,vars=pan[[k]],layer="data"); rowMeans(e) })
cl<-as.character(obj$pheno_cluster)
prof<-aggregate(sc, list(cluster=cl), mean)
assign_lin<-apply(prof[,names(pan)],1,function(r) names(pan)[which.max(r)])
prof$lineage<-assign_lin
cat("[2] per-cluster lineage (Epi kept):\n"); print(prof[order(prof$lineage),c("cluster","lineage",names(pan))],digits=2)
epi_cl<-prof$cluster[prof$lineage=="Epi"]
cat("   epithelial clusters:",length(epi_cl),"of",nrow(prof),"\n")
keep_cl <- cl %in% epi_cl
## per-cell immune/RBC backstop (unambiguous non-tumor)
bad<-c("PTPRC","C1QA","C1QB","LYZ","CD3D","CD3E","NKG7","MS4A1","CD68","CD14","HBB","HBA1","HBA2")
bad<-bad[bad %in% rownames(obj)]
be<-FetchData(obj,vars=bad,layer="data"); contam<-rowSums(be>0)>=1
keep <- keep_cl & !contam
cat("[3] cells: total",length(keep)," drop non-epi clusters",sum(!keep_cl)," drop immune/RBC+ (in epi cl)",sum(keep_cl&contam)," -> KEEP",sum(keep),"\n")
obj2<-subset(obj, cells=colnames(obj)[keep])
cat("   clean CSPC=",sum(obj2$grp=="CSPC")," CRPC=",sum(obj2$grp=="CRPC"),"\n")
## ---- DEG ----
Idents(obj2)<-"grp"
cat("[4] FindMarkers Wilcoxon CRPC vs CSPC (tumor-only)\n")
dem<-FindMarkers(obj2, ident.1="CRPC", ident.2="CSPC", test.use="wilcox", logfc.threshold=0, min.pct=0.01)
dem$gene<-rownames(dem); write.csv(dem, file.path(OUT,"FindMarkers_CRPC_vs_CSPC_tumoronly.csv"), row.names=FALSE)
gl<-dem$avg_log2FC; names(gl)<-dem$gene; gl<-sort(gl,decreasing=TRUE)
cat("   top UP-CRPC:",paste(head(names(gl),12),collapse=","),"\n")
cat("   top UP-CSPC:",paste(tail(names(gl),12),collapse=","),"\n")
FSU<-intersect(readLines(file.path(FSdir,"FitnessPosSignature_human.txt")),names(gl))
FSD<-intersect(readLines(file.path(FSdir,"FitnessNegSignature_human.txt")),names(gl))
## clusterProfiler for FSU (std), one-sided fgsea for each
t2g<-rbind(data.frame(term="FSU_highfitness",gene=FSU),data.frame(term="FSD_lowfitness",gene=FSD))
set.seed(1); eg<-GSEA(gl,TERM2GENE=t2g,minGSSize=10,maxGSSize=3000,pvalueCutoff=1,eps=0,seed=TRUE)
cp<-as.data.frame(eg)
set.seed(1)
fu<-fgsea(list(FSU=FSU),gl,minSize=10,maxSize=3000,eps=0,scoreType="std")
fd<-fgsea(list(FSD=FSD),gl,minSize=10,maxSize=3000,eps=0,scoreType="neg")
fu2<-fgsea(list(FSU=FSU),gl,minSize=10,maxSize=3000,eps=0,scoreType="pos")
res<-data.frame(pathway=c("FSU_highfitness","FSD_lowfitness"),
   NES=c(fu$NES,fd$NES), pval=c(fu$pval,fd$pval), size=c(fu$size,fd$size),
   NES_oneside=c(fu2$NES,fd$NES), pval_oneside=c(fu2$pval,fd$pval))
write.csv(res, file.path(OUT,"GSEA_tumoronly_FSU_FSD.csv"), row.names=FALSE)
cat("\n===== TUMOR-ONLY GSEA (+ up CRPC[late], - up CSPC[early]) =====\n"); print(res)
cat("\nclusterProfiler(FSU std): "); print(cp[,c("ID","NES","pvalue")])
cat("FSD leadingEdge:",paste(head(fd$leadingEdge[[1]],15),collapse=","),"\n")
cat("FSU leadingEdge:",paste(head(fu2$leadingEdge[[1]],15),collapse=","),"\n")
cat("ALL DONE\n")

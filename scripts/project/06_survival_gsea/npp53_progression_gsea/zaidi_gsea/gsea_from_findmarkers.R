suppressPackageStartupMessages(library(clusterProfiler))
OUT<-getwd(); FSdir<-"/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/FitnessSignature"
d<-read.csv("FindMarkers_CRPC_vs_CSPC.csv")
gl<-d$avg_log2FC; names(gl)<-d$gene; gl<-sort(gl,decreasing=TRUE)
FSU<-readLines(file.path(FSdir,"FitnessPosSignature_human.txt")); FSD<-readLines(file.path(FSdir,"FitnessNegSignature_human.txt"))
t2g<-rbind(data.frame(term="FSU_highfitness",gene=FSU),data.frame(term="FSD_lowfitness",gene=FSD))
t2g<-t2g[t2g$gene %in% names(gl),]
set.seed(1)
eg<-GSEA(gl,TERM2GENE=t2g,minGSSize=10,maxGSSize=3000,pvalueCutoff=1,eps=0,nPermSimple=10000,seed=TRUE)
res<-as.data.frame(eg)
write.csv(res[,c("ID","setSize","enrichmentScore","NES","pvalue","p.adjust","core_enrichment")],file.path(OUT,"clusterprofiler_GSEA_FSU_FSD_final.csv"),row.names=FALSE)
saveRDS(eg,file.path(OUT,"gsea_egobj_final.rds"))
cat("\n===== clusterProfiler GSEA FINAL (rank=avg_log2FC; + up CRPC[late], - up CSPC[early]) =====\n")
print(res[,c("ID","setSize","NES","pvalue","p.adjust")])
for(id in res$ID){le<-strsplit(res$core_enrichment[res$ID==id],"/")[[1]]; cat("\n",id," leadingEdge(",length(le),"):",paste(head(le,20),collapse=", "),"\n")}
cat("\nALL DONE\n")

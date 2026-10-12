suppressPackageStartupMessages({library(Seurat); library(clusterProfiler); library(enrichplot); library(ggplot2); library(cowplot); library(ggplotify)})
Z<-"/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/24PNAS_Zaidi_MSK_PCa_atlas"
FSdir<-"/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/FitnessSignature"
OUT<-getwd()
cat("[1] load tumor.cells.rds\n"); t0<-Sys.time()
o<-readRDS(file.path(Z,"GSE264573_msk.integrated.remove.cellcycle.tumor.cells.rds"))
o$grp<-as.character(o@meta.data$subtype); o<-subset(o, subset=grp %in% c("CSPC","CRPC"))   # non-NE
DefaultAssay(o)<-"RNA"; Idents(o)<-"grp"
cat("   CSPC=",sum(o$grp=="CSPC")," CRPC=",sum(o$grp=="CRPC")," (",round(difftime(Sys.time(),t0,units="mins"),1),"min)\n")
cat("[2] FindMarkers CRPC vs CSPC, min.pct=0.2\n")
dem<-FindMarkers(o, ident.1="CRPC", ident.2="CSPC", test.use="wilcox", logfc.threshold=0, min.pct=0.2)
dem$gene<-rownames(dem); write.csv(dem,file.path(OUT,"FindMarkers_CRPC_vs_CSPC_pct02.csv"),row.names=FALSE)
gl<-dem$avg_log2FC; names(gl)<-dem$gene; gl<-sort(gl,decreasing=TRUE); n<-length(gl)
cat("   genes:",n," top UP-CRPC:",paste(head(names(gl),8),collapse=","),"\n")
FSU<-intersect(readLines(file.path(FSdir,"FitnessPosSignature_human.txt")),names(gl))
FSD<-intersect(readLines(file.path(FSdir,"FitnessNegSignature_human.txt")),names(gl))
cat("   FSU in-data:",length(FSU)," FSD in-data:",length(FSD),"\n")
t2g<-rbind(data.frame(term="FSU",gene=FSU),data.frame(term="FSD",gene=FSD))
set.seed(1); eg<-GSEA(gl,TERM2GENE=t2g,minGSSize=10,maxGSSize=3000,pvalueCutoff=1,eps=0,seed=TRUE); r<-eg@result
write.csv(r[,c("ID","setSize","NES","pvalue","p.adjust")],file.path(OUT,"GSEA_zaidi_pct02_FSU_FSD.csv"),row.names=FALSE)
cat("\n===== Zaidi GSEA (min.pct=0.2; +NES=up in CRPC/late) =====\n"); print(r[,c("ID","setSize","NES","pvalue","p.adjust")])
# gseaplot2 split figure (match Cheng style)
mkplot<-function(id,col){nes<-r$NES[r$ID==id]; fdr<-r$p.adjust[r$ID==id]
  p<-gseaplot2(eg, geneSetID=id, title=sprintf("%s   NES=%.2f  FDR=%.3f",id,nes,fdr), color=col, base_size=11, pvalue_table=FALSE, rel_heights=c(1.5,0.4,0.9))
  p[[1]]<-p[[1]] + geom_hline(yintercept=0, linetype="dashed", colour="grey45", linewidth=0.4)
  p[[3]]<-p[[3]] + scale_x_continuous(expand=expansion(mult=0.01)) + xlab(NULL) +
    annotate("text",x=1,y=-Inf,label="CRPC",hjust=0,vjust=2.2,fontface="bold",size=3.6,colour="#b3261e") +
    annotate("text",x=n,y=-Inf,label="CSPC (Primary)",hjust=1,vjust=2.2,fontface="bold",size=3.6,colour="#1f5fa8") +
    coord_cartesian(clip="off") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=margin(5,8,20,5))
  as.ggplot(p)}
pf<-mkplot("FSU","#b45309"); pd<-mkplot("FSD","#2563eb")
comb<-plot_grid(pf,pd,ncol=2)
mt<-ggdraw()+draw_label("FSU/FSD in Zaidi MSK Cohort (PNAS, 2024)",fontface="bold",size=16)
comb<-plot_grid(mt,comb,ncol=1,rel_heights=c(0.07,1))
ggsave(file.path(OUT,"gseaplot2_zaidi_FSU_FSD_split.png"),comb,width=13,height=6.5,dpi=160,bg="white")
ggsave(file.path(OUT,"gseaplot2_zaidi_FSU_FSD_split.pdf"),comb,width=13,height=6.5)
cat("saved figure\nALL DONE\n")

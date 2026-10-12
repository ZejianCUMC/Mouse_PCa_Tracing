suppressPackageStartupMessages({library(Seurat); library(clusterProfiler); library(enrichplot); library(ggplot2); library(cowplot); library(ggplotify)})
Z<-"/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/24PNAS_Zaidi_MSK_PCa_atlas"
FSdir<-"/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/FitnessSignature"
OUT<-getwd()
o<-readRDS(file.path(Z,"GSE264573_msk.integrated.remove.cellcycle.tumor.cells.rds"))
sub<-as.character(o@meta.data$subtype)
o$grp<-ifelse(sub=="CSPC","Primary",ifelse(sub %in% c("CRPC","NEPC"),"Advanced",NA))
o<-subset(o, subset=grp %in% c("Primary","Advanced")); DefaultAssay(o)<-"RNA"; Idents(o)<-"grp"
cat("Primary(CSPC)=",sum(o$grp=="Primary")," Advanced(CRPC+NEPC)=",sum(o$grp=="Advanced")," [CRPC",sum(sub=="CRPC"),"+NEPC",sum(sub=="NEPC"),"]\n")
dem<-FindMarkers(o, ident.1="Advanced", ident.2="Primary", test.use="wilcox", logfc.threshold=0, min.pct=0.2)
dem$gene<-rownames(dem); write.csv(dem,file.path(OUT,"FindMarkers_Advanced_vs_Primary_withNE_pct02.csv"),row.names=FALSE)
gl<-dem$avg_log2FC; names(gl)<-dem$gene; gl<-sort(gl,decreasing=TRUE); n<-length(gl)
FSU<-intersect(readLines(file.path(FSdir,"FitnessPosSignature_human.txt")),names(gl))
FSD<-intersect(readLines(file.path(FSdir,"FitnessNegSignature_human.txt")),names(gl))
t2g<-rbind(data.frame(term="FSU",gene=FSU),data.frame(term="FSD",gene=FSD))
set.seed(1); eg<-GSEA(gl,TERM2GENE=t2g,minGSSize=10,maxGSSize=3000,pvalueCutoff=1,eps=0,seed=TRUE); r<-eg@result
write.csv(r[,c("ID","setSize","NES","pvalue","p.adjust")],file.path(OUT,"GSEA_zaidi_withNE_FSU_FSD.csv"),row.names=FALSE)
cat("\n===== Zaidi +NEPC GSEA (min.pct=0.2; +NES=up in Advanced/late incl NE) =====\n"); print(r[,c("ID","setSize","NES","pvalue","p.adjust")])
mkplot<-function(id,col){nes<-r$NES[r$ID==id]; fdr<-r$p.adjust[r$ID==id]
  p<-gseaplot2(eg, geneSetID=id, title=sprintf("%s   NES=%.2f  FDR=%.3g",id,nes,fdr), color=col, base_size=11, pvalue_table=FALSE, rel_heights=c(1.5,0.4,0.9))
  p[[1]]<-p[[1]] + geom_hline(yintercept=0, linetype="dashed", colour="grey45", linewidth=0.4)
  p[[3]]<-p[[3]] + scale_x_continuous(expand=expansion(mult=0.01)) + xlab(NULL) +
    annotate("text",x=1,y=-Inf,label="CRPC+NEPC (Advanced)",hjust=0,vjust=2.2,fontface="bold",size=3.4,colour="#b3261e") +
    annotate("text",x=n,y=-Inf,label="CSPC (Primary)",hjust=1,vjust=2.2,fontface="bold",size=3.4,colour="#1f5fa8") +
    coord_cartesian(clip="off") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=margin(5,8,20,5))
  as.ggplot(p)}
pf<-mkplot("FSU","#b45309"); pd<-mkplot("FSD","#2563eb"); comb<-plot_grid(pf,pd,ncol=2)
mt<-ggdraw()+draw_label("FSU/FSD in Zaidi MSK Cohort incl. NEPC (PNAS, 2024)",fontface="bold",size=16)
comb<-plot_grid(mt,comb,ncol=1,rel_heights=c(0.07,1))
ggsave(file.path(OUT,"gseaplot2_zaidi_withNE_FSU_FSD_split.png"),comb,width=13,height=6.5,dpi=160,bg="white")
ggsave(file.path(OUT,"gseaplot2_zaidi_withNE_FSU_FSD_split.pdf"),comb,width=13,height=6.5)
cat("saved\nALL DONE\n")

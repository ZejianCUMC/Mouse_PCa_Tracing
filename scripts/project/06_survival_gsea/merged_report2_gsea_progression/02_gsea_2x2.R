# Progression GSEA of the current FSU / FSD (01_genesets.py) for one tumor model (usage: 02_gsea_2x2.R <Npp53|Npp53Rb1>). Copy of
# Npp53_tracing/.../progression_genesets/mk_joint_report.R with only the gene sets and outputs changed: same ranked CRPC-vs-primary tumor-cell
# FindMarkers lists (Cheng, Zaidi; NEPC excluded / included), clusterProfiler GSEA (minGSSize 10, maxGSSize 3000, eps 0, seed), gseaplot2. Env: rMAST.
# nPermSimple = 10000 added (fgsea could not compute one Npp53Rb1 p value with the default); not used in the original Fig 9.1 run.
suppressPackageStartupMessages({library(clusterProfiler); library(enrichplot); library(ggplot2); library(cowplot); library(ggplotify)})
O<-dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(FALSE),value=TRUE)))); MODEL<-commandArgs(TRUE)[1]
FSU_g<-readLines(file.path(O,paste0(MODEL,"_FSU_human.txt"))); FSD_g<-readLines(file.path(O,paste0(MODEL,"_FSD_human.txt")))
RES<-list(); RED<-"#d62728"; BLUE<-"#1f5fa8"
one_gsea<-function(csv){dem<-read.csv(csv); gl<-dem$avg_log2FC; names(gl)<-dem$gene; gl<-sort(gl,decreasing=TRUE)
  t2g<-rbind(data.frame(term="FSU",gene=intersect(FSU_g,names(gl))),data.frame(term="FSD",gene=intersect(FSD_g,names(gl))))
  set.seed(1); eg<-GSEA(gl,TERM2GENE=t2g,minGSSize=10,maxGSSize=3000,pvalueCutoff=1,eps=0,seed=TRUE,nPermSimple=10000); RES[[basename(csv)]]<<-data.frame(ranking=basename(csv),eg@result[,c("ID","setSize","NES","pvalue","p.adjust")]); list(eg=eg,n=length(gl))}
panel<-function(eg,id,col,n,rlab){r<-eg@result; nes<-r$NES[r$ID==id]; fdr<-r$p.adjust[r$ID==id]
  p<-gseaplot2(eg,geneSetID=id,title=sprintf("%s  NES=%.2f  FDR=%.3g",id,nes,fdr),color=col,base_size=12,pvalue_table=FALSE,rel_heights=c(1.5,0.4,0.9))
  p[[1]]<-p[[1]]+geom_hline(yintercept=0,linetype="dashed",colour="grey45",linewidth=0.4)+theme(plot.title=element_text(size=12.5,face="bold"))
  p[[3]]<-p[[3]]+scale_x_continuous(expand=expansion(mult=0.01))+xlab(NULL)+
    annotate("text",x=1,y=-Inf,label="CRPC",hjust=0,vjust=2.2,fontface="bold",size=3.6,colour="#b3261e")+
    annotate("text",x=n,y=-Inf,label=rlab,hjust=1,vjust=2.2,fontface="bold",size=3.6,colour="#1f5fa8")+
    coord_cartesian(clip="off")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=margin(3,7,17,3))
  as.ggplot(p)}
pair<-function(csv,subtitle,rlab){g<-one_gsea(csv)
  row<-plot_grid(panel(g$eg,"FSU",RED,g$n,rlab),panel(g$eg,"FSD",BLUE,g$n,rlab),ncol=2)
  plot_grid(ggdraw()+draw_label(subtitle,fontface="bold",size=14,colour="#0f2747"),row,ncol=1,rel_heights=c(0.10,1))}
CG<-"/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/progression_genesets/cheng_gsea"
ZG<-"/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/progression_genesets/zaidi_gsea"
A<-pair(file.path(CG,"FindMarkers_CRPC_vs_TreatNaive_nonNE_pct02.csv"),"Cheng (Eur Urol 2022): NEPC excluded","Primary (Treatment-naive)")
B<-pair(file.path(CG,"FindMarkers_CRPC_vs_TreatNaive_pct02.csv"),      "Cheng (Eur Urol 2022): NEPC included","Primary (Treatment-naive)")
Cc<-pair(file.path(ZG,"FindMarkers_CRPC_vs_CSPC_pct02.csv"),           "Zaidi MSK (PNAS 2024): NEPC excluded","Primary (CSPC)")
D<-pair(file.path(ZG,"FindMarkers_Advanced_vs_Primary_withNE_pct02.csv"),"Zaidi MSK (PNAS 2024): NEPC included","Primary (CSPC)")
grid<-plot_grid(A,B,Cc,D,ncol=2,labels=c("a","b","c","d"),label_size=19)
ggsave(file.path(O,sprintf("gsea_progression_2x2_%s.png",MODEL)),grid,width=14,height=10,dpi=135,bg="white")
write.csv(do.call(rbind,RES),file.path(O,sprintf("gsea_progression_%s.csv",MODEL)),row.names=FALSE); print(do.call(rbind,RES)[,c("ID","setSize","NES","p.adjust")])

suppressPackageStartupMessages({library(Matrix);library(Seurat);library(ggplot2)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1);out<-normalizePath(args[1])
counts<-as(readMM(gzfile(file.path(out,'dotplot_counts_PRIVATE.mtx.gz'))),'CsparseMatrix')
genes<-read.csv(file.path(out,'dotplot_genes.csv'))$gene
meta<-read.csv(file.path(out,'final_metadata_PRIVATE.csv'),row.names=1,check.names=FALSE)
states<-read.csv(file.path(out,'state_order.csv'))$state
stopifnot(ncol(counts)==nrow(meta),nrow(counts)==length(genes),!any(c('GFP','mCherry')%in%genes),all(meta$total_counts>0))
dimnames(counts)<-list(genes,rownames(meta));meta$putative_state<-factor(meta$putative_state,levels=rev(states))
obj<-CreateSeuratObject(counts=counts,meta.data=meta,min.cells=0,min.features=0)
norm<-counts%*%Diagonal(x=10000/meta$total_counts);norm@x<-log1p(norm@x);dimnames(norm)<-dimnames(counts)
LayerData(obj,assay='RNA',layer='data')<-norm;Idents(obj)<-'putative_state'
panels<-read.csv(file.path(out,'dotplot_panels.csv'),check.names=FALSE)
stopifnot(identical(panels$gene,genes),!anyDuplicated(genes))
 selected<-genes
 p<-DotPlot(obj,features=selected,cols=c('lightgrey','#08306b'),dot.scale=3.3,scale=TRUE,cluster.idents=FALSE,col.min=-2.5,col.max=2.5,scale.by='radius')+RotatedAxis()+theme_bw(base_size=12)+theme(axis.text.x=element_text(size=10,angle=90,hjust=1,vjust=.5,face='italic'),axis.text.y=element_text(size=10),legend.position='bottom',legend.box='horizontal',legend.margin=margin(t=8))+labs(x='Gene',y=NULL)+guides(size=guide_legend(title='Percent expressed',nrow=1,order=1),colour=guide_colourbar(title='Scaled expression',direction='horizontal',barwidth=grid::unit(4,'cm'),barheight=grid::unit(.35,'cm'),order=2))
 p<-p+scale_radius(range=c(0,3.3),limits=c(0,100),breaks=c(25,50,75))+scale_colour_gradient(low='lightgrey',high='#08306b',limits=c(-2.5,2.5),breaks=c(-2,-1,0,1,2))
if(length(unique(panels$panel))>1){
 p<-p+theme(axis.text.x=element_text(size=9,angle=45,hjust=1,vjust=1,face='italic'),plot.margin=margin(t=100,r=45,b=10,l=10))+
   coord_cartesian(ylim=c(.5,length(states)+.5),clip='off')
 for(group in unique(panels$panel)){
  ix<-which(panels$panel==group)
  p<-p+annotate('segment',x=min(ix)-.35,xend=max(ix)+.35,y=length(states)+.67,yend=length(states)+.67,linewidth=.5)+
   annotate('text',x=min(ix)-.2,y=length(states)+.85,label=group,angle=45,hjust=0,vjust=0,size=3.2,fontface='bold')
 }
}
ggsave(file.path(out,'final_marker_dotplot.png'),p,width=if(length(genes)>80)16 else 14,height=max(7.5,3.5+.33*length(states)),dpi=180,limitsize=FALSE)
write.csv(p$data,file.path(out,'dotplot_aggregate_values.csv'),row.names=FALSE)
cat('Actual Seurat DotPlot complete:',length(states),'states,',length(genes),'genes\n')

# Figure 2.2: existing joint Harmony theta1 UMAP, split by treatment, not refitted.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
run <- normalizePath(args[1], mustWork=TRUE)
coords <- read.csv(file.path(run,'UMAP_Harmony_theta1.csv'),row.names=1,check.names=FALSE)
meta <- read.csv(file.path(run,'review/marker_first_all_cell_annotations.csv'),row.names=1,check.names=FALSE)
stopifnot(nrow(coords)==56686L, !anyDuplicated(rownames(coords)), setequal(rownames(coords),rownames(meta)))
meta <- meta[rownames(coords),,drop=FALSE]
stopifnot(length(unique(meta$sample_id))==22L,
          !any(meta$sample_id %in% c('JZ257','JZ275','JZ278','JZ281')),
          min(table(meta$sample_id))>=400L, all(is.finite(as.matrix(coords))))
meta$treatment[meta$treatment=='Castration+Enzalutamide'] <- 'Castration+ENZA'
treatments <- c('Intact','Castration','Castration+ENZA')
stopifnot(setequal(unique(meta$treatment),treatments))
colours <- c(MJZ005='#1f77b4',MJZ008='#ff7f0e',MJZ025='#2ca02c',
             MJZ026='#d62728',MJZ027='#9467bd',MJZ028='#8c564b')
stopifnot(all(meta$library %in% names(colours)))
counts <- data.frame(treatment=treatments,cells=sapply(treatments,function(t)sum(meta$treatment==t)),
                    mice=sapply(treatments,function(t)length(unique(meta$sample_id[meta$treatment==t]))))
stopifnot(identical(as.integer(counts$cells),c(37413L,5420L,13853L)),
          identical(as.integer(counts$mice),c(8L,6L,8L)))
write.csv(counts,file.path(run,'harmony_theta1_treatment_counts.csv'),row.names=FALSE)
set.seed(19)
draw_order <- sample.int(nrow(coords))
xlim <- range(coords[[1]]);ylim <- range(coords[[2]])
png(file.path(run,'harmony_theta1_by_treatment.png'),width=3300,height=1400,res=180,pointsize=16,type='cairo')
par(mfrow=c(1,3),cex=1,mar=c(2.5,2.7,3.5,.5),oma=c(4,0,4,0),family='sans',las=1,mgp=c(1.2,.3,0))
for(t in treatments){
  plot(coords[[1]][draw_order],coords[[2]][draw_order],type='n',xlim=xlim,ylim=ylim,
       asp=1,xlab='UMAP_1',ylab='UMAP_2',axes=FALSE,bty='n',cex.lab=1.15,
       main='')
  mtext(t,side=3,line=1.8,cex=1.4,font=2)
  points(coords[[1]][draw_order],coords[[2]][draw_order],pch=16,cex=.22,col='#edf0f3')
  ix <- draw_order[meta$treatment[draw_order]==t]
  points(coords[[1]][ix],coords[[2]][ix],pch=16,cex=.48,col=colours[meta$library[ix]])
  j <- match(t,counts$treatment)
  mtext(sprintf('%s cells | %d mice',format(counts$cells[j],big.mark=',',trim=TRUE),counts$mice[j]),side=3,line=.3,cex=1.05)
}
mtext('Harmony theta1 | joint embedding split by treatment',side=3,outer=TRUE,line=2,cex=1.5,font=2)
mtext('22 retained mice | 56,686 cells | batch correction: library/run',side=3,outer=TRUE,line=.4,cex=1.05)
par(oma=c(0,0,0,0),fig=c(0,1,0,.11),new=TRUE,mar=c(0,0,0,0))
plot.new();legend(x=.5,y=.9,xjust=.5,yjust=1,legend=names(colours),col=colours,pch=16,pt.cex=1.3,horiz=TRUE,bty='n',cex=1.1,title='Colour: library / run')
dev.off()
print(counts)

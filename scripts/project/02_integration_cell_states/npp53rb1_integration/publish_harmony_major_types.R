# Descriptive marker review of the completed sweep; no new model fitting.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
run <- normalizePath(args[1],mustWork=TRUE)
r <- file.path(run,'harmony_major_types')
o <- read.csv(file.path(r,'resolution_cell_annotations.csv'),row.names=1,check.names=FALSE)
xy <- read.csv(file.path(run,'UMAP_Harmony_theta1.csv'),row.names=1,check.names=FALSE)
stopifnot(nrow(o)==56686, !anyDuplicated(rownames(o)), setequal(rownames(o),rownames(xy)))
o <- o[rownames(xy),,drop=FALSE]
stopifnot(length(unique(o$sample_id))==22, min(table(o$sample_id))>=400,
          !any(o$sample_id %in% c('JZ257','JZ275','JZ278','JZ281')))
cl <- as.character(o$leiden_r0p3)
stopifnot(setequal(unique(cl),as.character(0:12)))
labels <- c('Mesenchymal','Mesenchymal','Basal epithelial','Mesenchymal',
            'Myeloid','Myeloid','Fibroblast-like','Endothelial-like',
            'Luminal epithelial','Basal epithelial','Basal epithelial',
            'Basal epithelial','Mixed epithelial')
modifiers <- c('Cycling-enriched','Keratin-positive subset; heterogeneous',
 'Small cluster; mixed mesenchymal signal','Heterogeneous; includes stress/unresolved cells',
 'None','Cycling-enriched','None','Minor immune/mesenchymal signal',
 'Basal keratin co-expression','QC/stress-enriched; retain for review',
 'None','Cycling-enriched','Luminal/basal programs; do not force subtype')
reasons <- c('Vim/Zeb1/Col1a1/Fn1; Mki67/Top2a',
 'Vim/Zeb1/Col1a1; Krt8/Krt18-positive subset',
 'Krt5/Krt14/Trp63/Epcam; minority Vim',
 'Vim/Zeb1/Col1a1/Col3a1; Dcn',
 'Lyz2/Tyrobp/Ptprc; low epithelial markers',
 'Lyz2/Tyrobp/Ptprc; increased Mki67',
 'Col1a1/Dcn/Lum/Pdgfra',
 'Pecam1/Cdh5/Emcn; minority Lyz2',
 'Epcam/Krt8/Krt18/Krt19; Krt14 also detected',
 'Krt5/Krt14; lower gene counts, higher mitochondrial fraction',
 'Krt5/Krt14/Trp63/Epcam',
 'Krt5/Krt14/Trp63; Mki67/Top2a',
 'Krt8/Krt18/Epcam plus Krt5/Krt14 and minority Vim')
for (field in c('primary_cluster','major_type','putative_tumor')) {
  o[[paste0('source_fastMNN_',field)]] <- o[[field]]
  o[[field]] <- NULL
}
o$primary_cluster <- as.integer(cl)
o$marker_first_major_type <- labels[as.integer(cl)+1]
o$review_modifier <- modifiers[as.integer(cl)+1]
o$reporter_positive <- o$raw_GFP>0 | o$raw_mCherry>0
o$reporter_ge2 <- o$raw_GFP>=2 | o$raw_mCherry>=2
stopifnot(sum(o$reporter_positive)==39812, sum(o$reporter_ge2)==36530)
o$all_reporter_positive_interest <- o$reporter_positive
o$GFP_only <- o$raw_GFP>0 & o$raw_mCherry==0
o$mCherry_only <- o$raw_mCherry>0 & o$raw_GFP==0
o$double_positive <- o$raw_GFP>0 & o$raw_mCherry>0
write.csv(o,file.path(r,'marker_first_all_cell_annotations.csv'))
write.csv(o[o$reporter_positive,,drop=FALSE],file.path(r,'ALL_reporter_positive_cells_for_review.csv'))
write.csv(o[o$reporter_ge2,,drop=FALSE],file.path(r,'ALL_reporter_ge2_cells_for_review.csv'))
summarize <- function(keys,name) {
  counts <- data.frame(cells=rep(1,nrow(o)),reporter_positive=as.integer(o$reporter_positive),
    reporter_ge2=as.integer(o$reporter_ge2),GFP_only=as.integer(o$GFP_only),
    mCherry_only=as.integer(o$mCherry_only),double_positive=as.integer(o$double_positive))
  z <- aggregate(counts,by=o[keys],FUN=sum)
  z$reporter_percent <- 100*z$reporter_positive/z$cells
  stopifnot(sum(z$cells)==56686,sum(z$reporter_positive)==39812)
  write.csv(z,file.path(r,name),row.names=FALSE)
}
summarize('marker_first_major_type','major_types_then_reporters.csv')
summarize(c('primary_cluster','marker_first_major_type'),'clusters_then_reporters.csv')
summarize('sample_id','sample_reporters.csv')
summarize(c('sample_id','marker_first_major_type'),'sample_type_reporters.csv')
evidence <- do.call(rbind,lapply(0:12,function(k){
  z <- o[o$primary_cluster==k,,drop=FALSE]; dominant <- sort(table(z$sample_id),decreasing=TRUE)
  data.frame(cluster=k,cells=nrow(z),major_type=labels[k+1],modifier=modifiers[k+1],
    marker_reason=reasons[k+1],median_genes=median(z$n_genes_by_counts),
    median_mt=median(z$pct_counts_mt),dominant_mouse=names(dominant)[1],
    dominant_mouse_fraction=as.numeric(dominant[1])/nrow(z),reporter_positive=sum(z$reporter_positive))
}))
write.csv(evidence,file.path(r,'reviewed_cluster_evidence.csv'),row.names=FALSE)
palette <- c('#1f77b4','#70a9d3','#ff7f0e','#e5ab66','#24943a','#91c77d',
             '#c42d35','#e58a92','#7953a4','#bba0d1','#8c564b','#c3948b','#c253ad')
set.seed(19); ord <- sample.int(nrow(o))
baseplot <- function(cols,main='') {
  plot(xy[[1]],xy[[2]],type='n',axes=FALSE,bty='n',asp=1,xlab='UMAP_1',ylab='UMAP_2',main=main)
  points(xy[[1]][ord],xy[[2]][ord],pch=16,cex=.40,col=cols[ord])
}
png(file.path(r,'cluster_ids.png'),width=2800,height=2050,res=180,pointsize=18,type='cairo')
layout(matrix(c(1,2),1),widths=c(3.3,2));par(mar=c(2.7,2.7,5,.2),mgp=c(1.2,.3,0))
baseplot(palette[o$primary_cluster+1])
for(k in 0:12){
  pts <- as.matrix(xy[o$primary_cluster==k,]); center <- apply(pts,2,median)
  pos <- pts[which.min(rowSums(sweep(pts,2,center)^2)),]
  points(pos[1],pos[2],pch=21,bg='white',col=palette[k+1],cex=2.1,lwd=1.5)
  text(pos[1],pos[2],labels=k,cex=1.05,font=2)
}
mtext('Harmony theta1 | Leiden resolution 0.3',3,line=2.8,cex=1.15,font=2)
mtext('56,686 cells | 22 mice | 13 clusters',3,line=1.1,cex=.95)
par(mar=c(2,0,5,.5));plot.new()
legend('topleft',legend=sprintf('%d  %s\n     %s cells',0:12,labels,format(evidence$cells,big.mark=',',trim=TRUE)),
       col=palette,pch=16,pt.cex=1.3,cex=.84,bty='n',y.intersp=1.38,title='Cluster ID / preliminary identity')
dev.off()
means <- read.csv(file.path(r,'leiden_r0p3_marker_means.csv'),row.names=1,check.names=FALSE)
fracs <- read.csv(file.path(r,'leiden_r0p3_marker_fractions.csv'),row.names=1,check.names=FALSE)
genes <- c('Krt5','Krt14','Trp63','Krt8','Krt18','Krt19','Epcam','Vim','Zeb1',
 'Col1a1','Dcn','Lum','Pdgfra','Ptprc','Lyz2','Tyrobp','Cd3e','Ms4a1','Pecam1',
 'Cdh5','Emcn','Mki67','Top2a','GFP','mCherry')
m <- as.matrix(means[as.character(0:12),genes]);f <- as.matrix(fracs[as.character(0:12),genes])
scaled <- sweep(m,2,apply(m,2,max),'/');scaled[!is.finite(scaled)] <- 0
png(file.path(r,'marker_dotplot.png'),width=3000,height=1850,res=180,pointsize=17,type='cairo')
par(mar=c(8,15,4,2),mgp=c(3,.7,0),tcl=-.3)
plot(NA,xlim=c(.5,length(genes)+.5),ylim=c(.4,13.6),axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i')
box(bty='l',lwd=1.2)
axis(1,at=seq_along(genes),labels=genes,las=2,tick=TRUE,lwd=0,lwd.ticks=1.2)
axis(2,at=13:1,labels=paste(0:12,labels,sep=' | '),las=2,tick=TRUE,cex.axis=.85,lwd=0,lwd.ticks=1.2)
mtext('Marker gene',side=1,line=4.7,cex=1.05,font=2)
mtext('Cluster ID | preliminary cell type',side=2,line=12.7,cex=1.05,font=2)
ramp <- colorRampPalette(c('#eff4fa','#145792'))(101)
for(i in 1:13)points(seq_along(genes),rep(14-i,length(genes)),pch=21,
  cex=2.5*sqrt(f[i,]),bg=ramp[1+round(100*scaled[i,])],col='#a9b2bc',lwd=.5)
title('Harmony theta1 / Leiden 0.3: canonical-marker evidence',cex.main=1.2)
mtext('Dot area: fraction expressing | Blue: mean log-normalized expression / gene maximum | Reporters overlaid last',1,line=6.5,cex=.78)
dev.off()
types <- unique(labels); typecols <- setNames(c('#c87f2a','#426bb2','#38924e','#98784f','#c04e64','#9264b5','#767676'),types)
tab <- prop.table(table(o$sample_id,factor(o$marker_first_major_type,levels=types)),1)
png(file.path(r,'composition_by_mouse.png'),width=3000,height=1400,res=180,pointsize=16,type='cairo')
layout(matrix(c(1,2),1),widths=c(3.3,1));par(mar=c(5,4,3,.5))
barplot(t(tab),col=typecols[types],border=NA,las=2,ylab='Fraction per mouse',main='Harmony theta1 / Leiden 0.3 | preliminary major types')
par(mar=c(5,0,3,.5));plot.new()
legend('topleft',legend=types,fill=typecols[types],bty='n',cex=.9)
dev.off()
png(file.path(r,'reporters_QC.png'),width=3300,height=2100,res=180,pointsize=17,type='cairo')
par(mfrow=c(2,3),mar=c(2.6,2.8,3.2,.8),mgp=c(1.2,.3,0))
heat <- colorRampPalette(c('#e2e5e9','#efad48','#b02027'))(101)
for(field in c('raw_GFP','raw_mCherry','reporter_positive','n_genes_by_counts','pct_counts_mt','primary_cluster')){
  if(field=='reporter_positive') {cols <- ifelse(o[[field]],'#be3446','#d8dde4'); title <- 'GFP or mCherry >0: 39,812'}
  else if(field=='primary_cluster') {cols <- palette[o$primary_cluster+1];title <- 'Harmony cluster IDs: 0-12'}
  else {v <- o[[field]];if(grepl('raw_',field))v<-log1p(v)
    cap<-unname(quantile(v,.99));cols<-heat[1+round(100*pmin(v/cap,1))]
    title<-paste0(sub('raw_','log1p ',field),' (0-',round(cap,1),')')}
  baseplot(cols,title)
}
dev.off()
print(evidence);print(table(o$marker_first_major_type));cat('Harmony publication artifacts complete\n')

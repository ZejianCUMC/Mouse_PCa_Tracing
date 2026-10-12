suppressPackageStartupMessages({library(Seurat);library(Matrix);library(harmony)})
set.seed(19);out<-file.path(Sys.getenv('RB1_RUN'),'tumor_umap_20260921')
m<-as(readMM(file.path(out,'tumor_hvg_counts.mtx')),'CsparseMatrix')
md<-read.csv(file.path(out,'metadata.csv'),row.names=1,check.names=FALSE)
rownames(m)<-readLines(file.path(out,'hvg.txt'));colnames(m)<-rownames(md)
so<-CreateSeuratObject(m,meta.data=md)
# Normalize by full-transcriptome counts, not by the selected-HVG library size.
ln<-m %*% Diagonal(x=1e4/md$total_counts);ln@x<-log1p(ln@x)
LayerData(so,assay='RNA',layer='data')<-ln
VariableFeatures(so)<-rownames(m)
so<-ScaleData(so,features=rownames(m),verbose=FALSE)
so<-RunPCA(so,features=rownames(m),npcs=30,verbose=FALSE)
write.csv(Embeddings(so,'pca'),file.path(out,'latent_uncorrected.csv'))
configs<-list(Harmony_run=list(vars='run_id',theta=1),
              Harmony_sample=list(vars='sample_id',theta=1),
              Harmony_run_sample=list(vars=c('run_id','sample_id'),theta=c(1,.25)))
for(name in names(configs)){
  cfg<-configs[[name]];set.seed(19)
  fit<-RunHarmony(so,group.by.vars=cfg$vars,theta=cfg$theta,verbose=TRUE)
  write.csv(Embeddings(fit,'harmony'),file.path(out,paste0('latent_',name,'.csv')))
  cat('COMPLETED',name,'\n')
}
writeLines(capture.output(sessionInfo()),file.path(out,'harmony_sessionInfo.txt'))

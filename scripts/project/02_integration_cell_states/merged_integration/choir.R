suppressPackageStartupMessages({library(Matrix);library(Seurat);library(CHOIR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
# Tested one-token parser fix from Npp53Rb1; never edit installed package.
original<-getFromNamespace('pruneTree','CHOIR')
old<-'P\\d*_L\\d*_';new<-'^P\\d*_L-?\\d*_';changes<-0L
replace_literal<-function(node,from,to){
 if(is.character(node)&&length(node)==1L&&identical(node,from)){changes<<-changes+1L;return(to)}
 if(is.call(node)||is.pairlist(node)||is.expression(node)) for(i in seq_along(node))
  if(!identical(node[[i]],quote(expr=))) node[i]<-list(replace_literal(node[[i]],from,to))
 node
}
patched<-original;body(patched)<-replace_literal(body(original),old,new);stopifnot(changes==1L)
changes<-0L;restored<-replace_literal(body(patched),new,old);stopifnot(changes==1L,identical(restored,body(original)))
stopifnot(identical(as.numeric(sub(new,'',c('P0_L0_32','P0_L-1_2','P0_L0000_7','P12_L45_19'))),c(32,2,7,19)))
if(args[1]=='--test-only'){cat('CHOIR parser regression PASS; package unchanged\n');quit(status=0)}
out<-normalizePath(args[1])
# CHOIR inputs are exported before the independent Leiden sweep starts.
plan_file<-file.path(out,'PREPARED.json')
if(!file.exists(plan_file)) plan_file<-file.path(dirname(out),'INTEGRATION_SELECTION.json')
plan<-fromJSON(plan_file)
stopifnot(!file.exists(file.path(out,'CHOIR_COMPLETE.json')))
runner<-getFromNamespace('CHOIR','CHOIR');e<-new.env(parent=environment(runner));e$pruneTree<-patched;environment(runner)<-e
mat<-as(readMM(gzfile(file.path(out,'hvg_counts_PRIVATE.mtx.gz'))),'CsparseMatrix')
rownames(mat)<-readLines(file.path(out,'features.txt'));colnames(mat)<-readLines(file.path(out,'barcodes_PRIVATE.txt'))
meta<-read.csv(file.path(out,'metadata_PRIVATE.csv'),row.names=1,check.names=FALSE)
latent<-as.matrix(read.csv(file.path(out,'scvi_latent_PRIVATE.csv'),row.names=1))
stopifnot(ncol(mat)==plan$cells,nrow(mat)==3000,ncol(latent)==30,identical(colnames(mat),rownames(meta)),identical(rownames(latent),rownames(meta)),all(is.finite(latent)))
seu<-CreateSeuratObject(counts=mat,meta.data=meta,min.cells=0,min.features=0)
norm<-mat%*%Diagonal(x=10000/meta$total_counts);norm@x<-log1p(norm@x);dimnames(norm)<-dimnames(mat)
LayerData(seu,assay='RNA',layer='data')<-norm;VariableFeatures(seu)<-rownames(mat)
colnames(latent)<-paste0('scVI_',seq_len(ncol(latent)))
seu[['scVI']]<-CreateDimReducObject(embeddings=latent,key='scVI_',assay='RNA')
RNGkind("L'Ecuyer-CMRG");set.seed(42);warnings<-character()
cat('CHOIR START',plan$cells,format(Sys.time()),'\n');flush.console()
seu<-withCallingHandlers(runner(object=seu,alpha=.05,p_adjust='bonferroni',feature_set='var',reduction=latent,var_features=rownames(mat),subtree_reductions=FALSE,normalization_method='none',batch_correction_method='none',use_assay='RNA',use_slot='data',neighbor_params=list(k.param=15),n_iterations=100,n_trees=50,n_cores=8,random_seed=42,verbose=TRUE),warning=function(w){
 warnings<<-c(warnings,conditionMessage(w));if(grepl('NAs introduced by coercion',conditionMessage(w),fixed=TRUE)) stop('Invalid CHOIR label coercion')
})
labels<-as.character(seu$CHOIR_clusters_0.05);groups<-seu@misc$CHOIR$clusters
final<-groups[['CHOIR_clusters_0.05']];steps<-groups[['stepwise_clusters_0.05']]
stopifnot(length(labels)==plan$cells,!anyNA(labels),nrow(final)==plan$cells,!anyNA(final),!any(grepl('_NA|_NaN|_Inf|-Inf',final$Record_cluster_label)))
stopifnot(!any(vapply(steps[-1],function(x)any(grepl('_NA|_NaN|_Inf|-Inf',x)),logical(1))))
write.csv(data.frame(cell=colnames(seu),CHOIR_cluster=labels),file.path(out,'choir_labels_PRIVATE.csv'),row.names=FALSE)
saveRDS(seu,file.path(out,'choir_HVG_object_PRIVATE.RDS'))
write_json(list(cells=plan$cells,clusters=length(unique(labels)),seed=42,RNG_kind=RNGkind(),cores=8,alpha=.05,p_adjust='bonferroni',iterations=100,trees=50,normalized_by='Original full-library counts',reference='Raw unpooled CHOIR labels',project_local_parser_fix=TRUE,invalid_stepwise_labels=0,installed_package_modified=FALSE,version=as.character(packageVersion('CHOIR')),warnings=warnings,completed=as.character(Sys.time())),file.path(out,'CHOIR_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
cat('CHOIR COMPLETE',length(unique(labels)),format(Sys.time()),'\n')

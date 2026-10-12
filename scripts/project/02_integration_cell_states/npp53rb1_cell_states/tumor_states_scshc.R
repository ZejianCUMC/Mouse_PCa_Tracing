# Genuine sc-SHC testClusters on the requested raw-CHOIR AMI candidates.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args) %in% c(1,2))
primary_seed <- if(length(args)==2) as.integer(args[2]) else 42L
out <- normalizePath(args[1]);.libPaths(c(file.path(out,'R_library'),file.path(dirname(out),'R_library'),.libPaths()))
suppressPackageStartupMessages({library(Matrix);library(scSHC);library(jsonlite)})
stopifnot(file.exists(file.path(out,'RANKED.json')), !file.exists(file.path(out,'scSHC_COMPLETE.json')))
mat <- as(readMM(gzfile(file.path(out,'hvg_counts_PRIVATE.mtx.gz'))),'CsparseMatrix')
rownames(mat)<-readLines(file.path(out,'features.txt'));colnames(mat)<-readLines(file.path(out,'barcodes_PRIVATE.txt'))
meta<-read.csv(file.path(out,'metadata_PRIVATE.csv'),row.names=1,check.names=FALSE)
assign<-read.csv(file.path(out,'resolution_assignments_PRIVATE.csv'),row.names=1,check.names=FALSE)
candidates<-read.csv(file.path(out,'scSHC_candidates.csv'))
expected_cells <- fromJSON(file.path(out,'RANKED.json'))$cells
stopifnot(ncol(mat)==expected_cells,identical(colnames(mat),rownames(meta)),identical(rownames(meta),rownames(assign)))
stopifnot(all(mat@x>=0),all(mat@x==round(mat@x)),length(unique(meta$run_id))==6)
# Progress logging only; the upstream statistic and null model are unmodified.
trace('test_split',where=asNamespace('scSHC'),print=FALSE,tracer=quote({
  cat('  Testing split:',length(ids1),'vs',length(ids2),'eligible cells;',
      length(unique(batch[c(ids1,ids2)])),'runs; alpha=',alpha_level,format(Sys.time()),'\n')
  flush.console()
}))
results<-list()
for(i in seq_len(nrow(candidates))){
  r<-candidates$resolution[i];tag<-sprintf('r%.1f',r)
  original<-as.character(assign[[sprintf('leiden_r%.1f_seed%d',r,primary_seed)]])
  cat('START sc-SHC',tag,'input clusters',length(unique(original)),format(Sys.time()),'\n');flush.console()
  begin<-Sys.time();RNGkind("L'Ecuyer-CMRG");set.seed(primary_seed)
  tested<-testClusters(data=mat,cluster_ids=original,batch=as.character(meta$run_id),
                        var.genes=rownames(mat),alpha=.05,num_PCs=30,parallel=TRUE,cores=4)
  stopifnot(length(tested[[1]])==expected_cells,!anyNA(tested[[1]]))
  saveRDS(tested,file.path(out,paste0('scSHC_',tag,'_PRIVATE.RDS')))
  write.csv(data.frame(cell=colnames(mat),leiden_cluster=original,scSHC_cluster=as.character(tested[[1]])),
            file.path(out,paste0('scSHC_',tag,'_labels_PRIVATE.csv')),row.names=FALSE)
  write.csv(as.data.frame.matrix(table(original,tested[[1]])),file.path(out,paste0('scSHC_',tag,'_merge_table.csv')))
  capture.output(print(tested[[2]]),file=file.path(out,paste0('scSHC_',tag,'_tree.txt')))
  results[[i]]<-data.frame(resolution=r,input_clusters=length(unique(original)),defended_clusters=length(unique(tested[[1]])),
                           retained_fraction=length(unique(tested[[1]]))/length(unique(original)),
                           elapsed_minutes=as.numeric(difftime(Sys.time(),begin,units='mins')))
  write.csv(do.call(rbind,results),file.path(out,'scSHC_summary.csv'),row.names=FALSE)
  cat('COMPLETE sc-SHC',tag,'\n');print(results[[i]]);flush.console()
}
write_json(list(cells=expected_cells,candidates=candidates$resolution,seed=primary_seed,RNG_kind=RNGkind(),alpha=.05,batch='run_id',PCs=30,cores=4,
                version=as.character(packageVersion('scSHC')),completed=as.character(Sys.time()),
                caveat='Splits require >20 cells from each side within an eligible run; untestable batch-specific/rare splits may merge. Merging is not proof of biological identity.'),
           file.path(out,'scSHC_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)

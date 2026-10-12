suppressPackageStartupMessages({library(Seurat);library(Matrix);library(loupeR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
out<-normalizePath(args[1]);inputs<-file.path(out,'loupe_inputs')
stopifnot(file.exists(file.path(out,'FINAL_PROJECTIONS_COMPLETE.json')))
manifest<-fromJSON(file.path(inputs,'COUNTS_COMPLETE.json'))
cat('Loading all-gene corrected final count matrix\n')
mat<-Read10X_h5(file.path(inputs,'final_counts.h5'),use.names=TRUE,unique.features=FALSE)
features<-read.csv(file.path(inputs,'features.csv'),check.names=FALSE)
meta<-read.csv(file.path(inputs,'metadata_PRIVATE.csv'),row.names=1,check.names=FALSE)
checks<-read.csv(file.path(inputs,'cell_count_checks_PRIVATE.csv'),check.names=FALSE)
stopifnot(nrow(mat)==33702,nrow(mat)==manifest$genes,ncol(mat)==40141,identical(rownames(mat),features$symbol),
          identical(colnames(mat),rownames(meta)),identical(checks$cell,colnames(mat)),
          all(Matrix::colSums(mat)==checks$sum_corrected_counts),all(mat@x==round(mat@x)))
lev<-unique(meta$putative_state);lev<-lev[order(as.integer(sub('c([0-9]+)_.*','\\1',lev)))]
stopifnot(length(lev)==16,!any(grepl('^c(17|18|19)_',lev)))
clusters<-list(Putative_cell_state=factor(meta$putative_state,levels=lev),
              Treatment=factor(meta$treatment,levels=c('Intact','Castration','Castration+ENZA')),
              Run_ID=factor(meta$run_id),Sample_ID=factor(meta$sample_id),Hash_ID=factor(meta$hashtag),
              Original_Leiden_r1_1=factor(meta$original_r1_1_cluster))
projections<-list()
for(nm in c('UMAP','tSNE','Diffmap')){
 fname<-switch(nm,UMAP='umap_PRIVATE.csv',tSNE='tsne_PRIVATE.csv',Diffmap='diffmap_PRIVATE.csv')
 df<-read.csv(file.path(out,fname),row.names=1,check.names=FALSE)
 stopifnot(identical(rownames(df),colnames(mat)),ncol(df)==2,all(is.finite(as.matrix(df))))
 projections[[nm]]<-as.matrix(df)
}
stopifnot(validate_count_mat(mat,feature_ids=features$gene_id)$success,
          validate_clusters(clusters,ncol(mat))$success,validate_projections(projections,ncol(mat))$success)
cat('All barcodes, count totals, state categories and three projections validated\n')
# Keep the exact louper input for an independent HDF5 audit. No hidden features,
# count normalization or coordinate changes occur during export.
h5file<-file.path(inputs,'louper_input_PRIVATE.h5')
ok<-loupeR:::create_hdf5(mat,clusters,projections,h5file,feature_ids=features$gene_id,seurat_obj_version=NULL)
stopifnot(ok$success)
target<-file.path(out,'Npp53Rb1_final_tumor_states_c1_c16.cloupe')
stopifnot(!file.exists(target))
result<-loupeR:::louper_create_cloupe(h5file,output_dir=out,output_name='Npp53Rb1_final_tumor_states_c1_c16',force=FALSE)
stopifnot(result$success,file.exists(target),file.info(target)$size>1e6)
write_json(list(cells=40141,genes=nrow(mat),states=16,projections=names(projections),
           categories=names(clusters),loupeR_version=as.character(packageVersion('loupeR')),
           counts_total=sum(mat),bytes=file.info(target)$size,count_checks=TRUE,
           projection_barcode_order_checked=TRUE,feature_ids_preserved=TRUE),
           file.path(out,'CLOUPE_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PASS: final cloupe built with UMAP, tSNE and Diffmap; c1-c16 only\n')

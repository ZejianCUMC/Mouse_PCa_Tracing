# CHOIR reference on the fixed tumour scVI latent; no refitting or filtering.
suppressPackageStartupMessages({library(Matrix); library(Seurat); library(CHOIR); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
out <- normalizePath(args[1])
stopifnot(file.exists(file.path(out,'PREPARED.json')), !file.exists(file.path(out,'CHOIR_COMPLETE.json')))
mat <- as(readMM(gzfile(file.path(out,'hvg_counts_PRIVATE.mtx.gz'))),'CsparseMatrix')
rownames(mat) <- readLines(file.path(out,'features.txt'))
colnames(mat) <- readLines(file.path(out,'barcodes_PRIVATE.txt'))
meta <- read.csv(file.path(out,'metadata_PRIVATE.csv'),row.names=1,check.names=FALSE)
latent <- as.matrix(read.csv(file.path(out,'scvi_latent_PRIVATE.csv'),row.names=1,check.names=FALSE))
expected_cells <- fromJSON(file.path(out,'PREPARED.json'))$cells
stopifnot(ncol(mat)==expected_cells,nrow(mat)==3000,ncol(latent)==30,
          identical(colnames(mat),rownames(meta)),identical(rownames(latent),rownames(meta)))
stopifnot(all(is.finite(latent)),all(meta$total_counts>0),all(colSums(mat)<=meta$total_counts))
seu <- CreateSeuratObject(counts=mat,meta.data=meta,min.cells=0,min.features=0)
# Preserve full-library normalization denominators, not HVG-only library sizes.
norm <- mat %*% Diagonal(x=10000/meta$total_counts)
norm@x <- log1p(norm@x)
dimnames(norm) <- dimnames(mat)
LayerData(seu,assay='RNA',layer='data') <- norm
VariableFeatures(seu) <- rownames(mat)
colnames(latent) <- paste0('scVI_',seq_len(ncol(latent)))
seu[['scVI']] <- CreateDimReducObject(embeddings=latent,key='scVI_',assay='RNA')
cat('CHOIR start:',format(Sys.time()),'; cells:',ncol(seu),'; genes:',nrow(seu),'\n')
RNGkind("L'Ecuyer-CMRG");set.seed(42)
seu <- CHOIR(object=seu,alpha=0.05,p_adjust='bonferroni',feature_set='var',
             reduction=latent,var_features=rownames(mat),subtree_reductions=FALSE,
             normalization_method='none',batch_correction_method='none',
             use_assay='RNA',use_slot='data',neighbor_params=list(k.param=15),
             n_iterations=100,n_trees=50,n_cores=8,random_seed=42,verbose=TRUE)
labels <- as.character(seu$CHOIR_clusters_0.05)
stopifnot(length(labels)==expected_cells,!anyNA(labels),identical(colnames(seu),rownames(meta)))
write.csv(data.frame(cell=colnames(seu),CHOIR_cluster=labels),file.path(out,'choir_labels_PRIVATE.csv'),row.names=FALSE)
write.csv(as.data.frame(table(labels)),file.path(out,'choir_cluster_sizes.csv'),row.names=FALSE)
saveRDS(seu,file.path(out,'choir_HVG_object_PRIVATE.RDS'))
write_json(list(cells=ncol(seu),clusters=length(unique(labels)),alpha=.05,p_adjust='bonferroni',
                iterations=100,trees=50,seed=42,RNG_kind=RNGkind(),cores=8,subtree_reductions=FALSE,
                normalized_by='Original full-library UMI totals',
                anchor='Raw CHOIR clusters; no marker-based pooling',
                caveat='CHOIR tests use normalized expression, not batch-corrected counts; sample/run specificity must be reviewed',
                version=as.character(packageVersion('CHOIR')),completed=as.character(Sys.time())),
           file.path(out,'CHOIR_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
cat('CHOIR COMPLETE:',length(unique(labels)),'raw clusters\n')

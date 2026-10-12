args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
out <- normalizePath(args[1]);fix <- file.path(out,'choir_label_fix_seed42')
stopifnot(file.exists(file.path(fix,'CHOIR_LABEL_FIX_VERIFIED.json')))
suppressPackageStartupMessages({library(Seurat);library(jsonlite)})
old <- readRDS(file.path(out,'choir_HVG_object_PRIVATE.RDS'))
corrected <- readRDS(file.path(fix,'choir_HVG_object_PRIVATE.RDS'))
same <- identical(old@misc$CHOIR$clusters$full_tree,corrected@misc$CHOIR$clusters$full_tree)
stopifnot(same,identical(colnames(old),colnames(corrected)),identical(rownames(old),rownames(corrected)))
write_json(list(same_initial_tree=same,cells=ncol(old),genes=nrow(old),
                original_object_md5=unname(tools::md5sum(file.path(out,'choir_HVG_object_PRIVATE.RDS'))),
                corrected_object_md5=unname(tools::md5sum(file.path(fix,'choir_HVG_object_PRIVATE.RDS'))),
                checked=as.character(Sys.time()),
                caveat='Initial tree identity is checked; this is not an independent full-pipeline bitwise replay.'),
           file.path(fix,'TREE_COMPARISON.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PASS: initial CHOIR tree, cell order and gene order are identical before the pruning parser correction.\n')

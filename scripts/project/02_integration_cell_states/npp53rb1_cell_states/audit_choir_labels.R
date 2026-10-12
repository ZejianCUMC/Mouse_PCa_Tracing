args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
out <- normalizePath(args[1])
suppressPackageStartupMessages({library(Seurat);library(jsonlite)})
path <- file.path(out,'choir_HVG_object_PRIVATE.RDS')
x <- readRDS(path)
groups <- x@misc$CHOIR$clusters
final <- groups[['CHOIR_clusters_0.05']]
steps <- groups[['stepwise_clusters_0.05']]
csv_path <- file.path(out,'choir_labels_PRIVATE.csv')
exported <- read.csv(csv_path,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(identical(exported$cell,final$CellID),
          identical(as.character(exported$CHOIR_cluster),as.character(final$CHOIR_clusters_0.05)))
bad <- function(v)is.na(v)|grepl('_NA|_NaN|_Inf|-Inf',v)
step_bad <- vapply(steps[-1],function(v)sum(bad(v)),integer(1))
affected <- Reduce(`|`,lapply(steps[-1],bad))
final_bad <- bad(final$Record_cluster_label)|is.na(final$CHOIR_clusters_0.05)
record <- list(valid=!any(affected)&!any(final_bad),cells=nrow(final),
               final_invalid_cells=sum(final_bad),ever_invalid_cells=sum(affected),
               invalid_cells_by_step=as.list(step_bad[step_bad>0]),
               object_md5=unname(tools::md5sum(path)),labels_csv_md5=unname(tools::md5sum(csv_path)),
               exported_labels_match_object=TRUE,checked=as.character(Sys.time()))
write_json(record,file.path(out,'CHOIR_LABEL_AUDIT.json'),pretty=TRUE,auto_unbox=TRUE)
write.csv(data.frame(cell=steps$CellID[affected]),file.path(out,'choir_invalid_internal_labels_PRIVATE.csv'),row.names=FALSE)
print(record)

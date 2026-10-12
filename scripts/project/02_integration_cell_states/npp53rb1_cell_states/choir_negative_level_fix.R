# Project-local, single-token parser correction; installed CHOIR is unchanged.
suppressPackageStartupMessages({library(CHOIR);library(jsonlite)})
original_prune <- getFromNamespace('pruneTree','CHOIR')
old_pattern <- 'P\\d*_L\\d*_'
new_pattern <- '^P\\d*_L-?\\d*_'
test_labels <- c('P0_L0_32','P0_L-1_2','P0_L0000_7','P12_L45_19')
stopifnot(anyNA(suppressWarnings(as.numeric(sub(old_pattern,'',test_labels)))))
stopifnot(identical(as.numeric(sub(new_pattern,'',test_labels)),c(32,2,7,19)))
replacement_count <- 0L
replace_pattern <- function(node) {
  if(is.character(node) && length(node)==1L && identical(node,old_pattern)) {
    replacement_count <<- replacement_count+1L
    return(new_pattern)
  }
  if(is.call(node) || is.pairlist(node) || is.expression(node)) {
    for(i in seq_along(node)) {
      # Missing function arguments are not evaluated by this traversal.
      if(!identical(node[[i]],quote(expr=))) node[i] <- list(replace_pattern(node[[i]]))
    }
  }
  node
}
patched_prune <- original_prune
body(patched_prune) <- replace_pattern(body(original_prune))
stopifnot(replacement_count==1L)
saved_patterns <- c(old_pattern,new_pattern)
old_pattern <- saved_patterns[2];new_pattern <- saved_patterns[1]
replacement_count <- 0L
restored_body <- replace_pattern(body(patched_prune))
stopifnot(replacement_count==1L,identical(restored_body,body(original_prune)))
old_pattern <- saved_patterns[1];new_pattern <- saved_patterns[2]
cat('PASS: reproduced negative-level coercion; corrected parser handles signed and zero-padded levels; exactly one literal changed.\n')
if('--test-only' %in% commandArgs(trailingOnly=TRUE)) quit(status=0)

wrapper <- getFromNamespace('CHOIR','CHOIR')
wrapper_environment <- new.env(parent=environment(wrapper))
wrapper_environment$pruneTree <- patched_prune
environment(wrapper) <- wrapper_environment
CHOIR <- wrapper
script_path <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)))
project <- dirname(dirname(script_path))
warning_log <- character()
withCallingHandlers(source(file.path(project,'scripts/tumor_states_choir.R'),local=.GlobalEnv),warning=function(w) {
  warning_log <<- c(warning_log,conditionMessage(w))
  if(grepl('NAs introduced by coercion',conditionMessage(w),fixed=TRUE))
    stop('Corrected CHOIR still produced a coercion warning; refuse acceptance')
})
# The original runner leaves the completed Seurat object in `seu` and output in `out`.
groups <- seu@misc$CHOIR$clusters
final <- groups[['CHOIR_clusters_0.05']]
steps <- groups[['stepwise_clusters_0.05']]
stopifnot(nrow(final)==42806,!anyNA(final),
          !any(grepl('_NA|_NaN|_Inf|-Inf',final$Record_cluster_label)))
stopifnot(!any(vapply(steps[-1],function(x)any(grepl('_NA|_NaN|_Inf|-Inf',x)),logical(1))))
record <- fromJSON(file.path(out,'CHOIR_COMPLETE.json'))
record$project_local_patch <- 'pruneTree label parser accepts negative L levels; one regex literal changed'
record$original_pattern <- old_pattern
record$corrected_pattern <- new_pattern
record$upstream_commit <- as.character(packageDescription('CHOIR')$RemoteSha)
record$warnings <- warning_log
write_json(record,file.path(out,'CHOIR_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
write_json(list(cells=42806,clusters=length(unique(final$CHOIR_clusters_0.05)),
                stepwise_invalid_labels=0,final_invalid_labels=0,replacements=replacement_count,
                installed_package_modified=FALSE,completed=as.character(Sys.time())),
           file.path(out,'CHOIR_LABEL_FIX_VERIFIED.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PASS: corrected full CHOIR rerun; no malformed labels in any saved pruning step.\n')

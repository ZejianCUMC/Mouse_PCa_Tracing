#!/usr/bin/env Rscript
# Filter cohort_cache.rds to protein-coding genes only (Gencode v44, n=20,020 symbols).
# Saves as cohort_cache_pc.rds; reports per-cohort n_genes_before/after.
suppressPackageStartupMessages(library(data.table))
plog <- function(...) cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n", sep="")

NPP53 <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
CACHE_IN  <- file.path(NPP53, "analysis/fitness_signatures_survival/cohort_cache.rds")
CACHE_OUT <- file.path(NPP53, "analysis/fitness_signatures_survival/cohort_cache_pc.rds")
PCFILE    <- file.path(NPP53, "analysis/fitness_signatures_survival/protein_coding_human.txt")

plog("Loading cache + PC gene list ...")
COH <- readRDS(CACHE_IN)
PC  <- unique(readLines(PCFILE)); PC <- PC[nchar(PC)>0]
plog("cache cohorts: ", paste(names(COH), collapse=", "))
plog("PC gene list: ", length(PC), " symbols")

for (nm in names(COH)) {
  co <- COH[[nm]]
  n_before <- nrow(co$mat)
  keep <- rownames(co$mat) %in% PC
  co$mat <- co$mat[keep, , drop=FALSE]
  n_after <- nrow(co$mat)
  plog(sprintf("%-24s  %6d → %6d  (%.1f%% retained; %d samples)",
               nm, n_before, n_after, 100*n_after/n_before, ncol(co$mat)))
  COH[[nm]] <- co
}

plog("Saving filtered cache → ", CACHE_OUT)
saveRDS(COH, CACHE_OUT)
plog("Done. Size = ",
     round(file.info(CACHE_OUT)$size/1e6, 1), " MB")

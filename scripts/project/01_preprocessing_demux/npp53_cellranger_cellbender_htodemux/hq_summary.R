suppressPackageStartupMessages({library(Seurat); library(dplyr)})
ROOT <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/seurat_objs"
cohorts <- c("MJZ019","MJZ021","MJZ022","MJZ023","MJZ024")
rows <- list()
for (co in cohorts) {
  f <- file.path(ROOT, paste0(co,"_per_sample_DF"), paste0(co,"_seurat_DF_clean.RDS"))
  if (!file.exists(f)) { cat("MISSING:", f, "\n"); next }
  md <- readRDS(f)@meta.data
  rows[[co]] <- md %>% group_by(sample_id, treatment, cohort_id) %>%
                summarise(n_hq_cells=n(), .groups="drop")
}
summ <- bind_rows(rows) %>% mutate(pass_ge1000 = n_hq_cells >= 1000) %>% arrange(cohort_id, sample_id)
write.csv(summ, file.path(ROOT,"hq_cell_summary.csv"), row.names=FALSE)
cat("\n=== HIGH-QUALITY CELLS PER SAMPLE (HTODemux-singlet + DoubletFinder, cellbender counts) ===\n")
print(as.data.frame(summ))
cat("\nExcluded (<1000):", paste(summ$sample_id[!summ$pass_ge1000], collapse=", "), "\n")
cat("PASS:", sum(summ$pass_ge1000), "/", nrow(summ), "samples\n")

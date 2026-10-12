#!/usr/bin/env Rscript
# GSEA of Fate 2 (c3+c4+c5) vs all other CARTA states using clusterProfiler,
# ranked by Seurat FindMarkers Wilcoxon DE stats (proper per-cell test, not
# pseudobulk mean-diff).
suppressPackageStartupMessages({
  library(Seurat)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(presto)   # fast Wilcoxon for Seurat
})

H5AD <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/data/npp53_merge/npp53_edited.h5ad"
OUT  <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12/biology_review"

message("[1/6] Reading h5ad → SingleCellExperiment ...")
sce <- readH5AD(H5AD)
message("     dims: ", paste(dim(sce), collapse=" x "))
message("[2/6] Converting SCE → Seurat ...")
message("     assays available: ", paste(assayNames(sce), collapse=", "))
counts_mat <- assay(sce, "logcounts")   # h5ad X was stored as logcounts by zellkonverter
seu <- CreateSeuratObject(counts = counts_mat, meta.data = as.data.frame(colData(sce)))
# Mark the assay as already-normalized (data slot), avoid re-log
seu[["RNA"]]$data <- counts_mat
STATES <- c('Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR','Luminal-Persistent',
            'Luminal-Squamous','Basal-like','Partial-EMT','Mes-Quiescent','Mes-Collagen',
            'Mes-Stemness','Mes-myCAF')
seu <- subset(seu, subset = putative_state %in% STATES)
Idents(seu) <- seu$putative_state
message("     kept ", ncol(seu), " cells across 12 states")

message("[3/6] DE analysis Fate 2 (c3+c4+c5) vs all others (Wilcoxon via presto) ...")
FATE2 <- c('Luminal-ARdep','Luminal-CR','Luminal-Persistent')
seu$grp <- ifelse(seu$putative_state %in% FATE2, "Fate2", "Other")
Idents(seu) <- seu$grp
de <- FindMarkers(seu, ident.1 = "Fate2", ident.2 = "Other",
                  logfc.threshold = 0, min.pct = 0.01, test.use = "wilcox",
                  min.cells.group = 5)
de$gene <- rownames(de)
message("     n genes tested: ", nrow(de))
# Rank statistic: signed -log10(p) * sign(avg_log2FC)
de$rank_stat <- -log10(pmax(de$p_val, 1e-300)) * sign(de$avg_log2FC)
de <- de[order(-de$rank_stat), ]
write_tsv(de, file.path(OUT, "de_fate2_vs_other_wilcox.tsv"))
message("     saved DE: de_fate2_vs_other_wilcox.tsv")

message("[4/6] Loading Hallmark + KEGG + Reactome from msigdbr (mouse) ...")
# msigdbr v25 API: use msigdbr with species + collection
hallmark <- msigdbr(species = "Mus musculus", collection = "H")
kegg     <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
if (nrow(kegg) == 0) kegg <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_MEDICUS")
reactome <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
message("     Hallmark: ", length(unique(hallmark$gs_name)), " sets")
message("     KEGG:     ", length(unique(kegg$gs_name)),     " sets")
message("     Reactome: ", length(unique(reactome$gs_name)), " sets")

message("[5/6] Running clusterProfiler::GSEA on Wilcoxon-ranked list ...")
gene_list <- setNames(de$rank_stat, de$gene)

run_gsea <- function(t2g, label) {
  t2g_df <- t2g[, c("gs_name","gene_symbol")]
  colnames(t2g_df) <- c("term","gene")
  res <- GSEA(gene_list, TERM2GENE = t2g_df,
              pvalueCutoff = 1.0, minGSSize = 10, maxGSSize = 500,
              seed = TRUE, verbose = FALSE, eps = 0)
  if (nrow(res@result) == 0) {
    message("     [", label, "] no enrichment")
    return(NULL)
  }
  df <- as.data.frame(res@result) %>% arrange(desc(NES))
  df$db <- label
  write_tsv(df, file.path(OUT, paste0("gsea_fate2_", label, "_result.tsv")))
  message("     [", label, "] n sets: ", nrow(df),
          "; n FDR<0.25: ", sum(df$p.adjust < 0.25))
  df
}

hall_res <- run_gsea(hallmark, "Hallmark")
kegg_res <- run_gsea(kegg,     "KEGG")
reac_res <- run_gsea(reactome, "Reactome")

message("[6/6] Printing top 15 UP and DOWN per database ...")
show_top <- function(res, label) {
  if (is.null(res)) return(invisible())
  cat("\n=== ", label, " — TOP 15 UP (Fate 2 enriched) ===\n", sep="")
  up <- res %>% arrange(desc(NES)) %>% head(15)
  for (i in seq_len(nrow(up))) {
    marker <- if (up$p.adjust[i] < 0.05) " * " else if (up$p.adjust[i] < 0.25) " ~ " else "   "
    cat(sprintf("%s NES=%+.2f  padj=%.4f  %s\n",
                marker, up$NES[i], up$p.adjust[i],
                substr(up$ID[i], 1, 90)))
  }
  cat("\n=== ", label, " — TOP 10 DOWN (Other enriched) ===\n", sep="")
  dn <- res %>% arrange(NES) %>% head(10)
  for (i in seq_len(nrow(dn))) {
    marker <- if (dn$p.adjust[i] < 0.05) " * " else if (dn$p.adjust[i] < 0.25) " ~ " else "   "
    cat(sprintf("%s NES=%+.2f  padj=%.4f  %s\n",
                marker, dn$NES[i], dn$p.adjust[i],
                substr(dn$ID[i], 1, 90)))
  }
}
show_top(hall_res, "Hallmark")
show_top(kegg_res, "KEGG")
show_top(reac_res, "Reactome")

message("\nDONE. Outputs at ", OUT, ":")
message("  de_fate2_vs_other_wilcox.tsv")
message("  gsea_fate2_Hallmark_result.tsv")
message("  gsea_fate2_KEGG_result.tsv")
message("  gsea_fate2_Reactome_result.tsv")

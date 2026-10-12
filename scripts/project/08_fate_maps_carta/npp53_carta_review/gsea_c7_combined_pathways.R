#!/usr/bin/env Rscript
# c7 Basal-like vs all others — clusterProfiler GSEA using the COMBINED
# 75-pathway library (50 Hallmark + 25 PCa/CRPC custom, proper ortholog
# mouse conversion). Companion to the pure-Hallmark result.
suppressPackageStartupMessages({
  library(Seurat); library(zellkonverter); library(SingleCellExperiment)
  library(clusterProfiler); library(dplyr); library(readr); library(presto)
})
H5AD <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/data/npp53_merge/npp53_edited.h5ad"
GMT  <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/aucell_arhi/pathway_genesets_mouse_via_ortholog.gmt"
OUT  <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12/biology_review"

STATES <- c('Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR','Luminal-Persistent',
            'Luminal-Squamous','Basal-like','Partial-EMT','Mes-Quiescent','Mes-Collagen',
            'Mes-Stemness','Mes-myCAF')

message("[1/4] Read h5ad + build Seurat obj")
sce <- readH5AD(H5AD)
counts_mat <- assay(sce, "logcounts")
seu <- CreateSeuratObject(counts = counts_mat, meta.data = as.data.frame(colData(sce)))
seu[["RNA"]]$data <- counts_mat
seu <- subset(seu, subset = putative_state %in% STATES)
seu$grp <- ifelse(seu$putative_state == "Basal-like", "c7", "Other")
Idents(seu) <- seu$grp
message("     kept ", ncol(seu), " cells (c7 = ", sum(seu$grp == "c7"), ")")

# Reuse DE if it already exists to save time
de_file <- file.path(OUT, "de_c7_vs_other_wilcox.tsv")
if (file.exists(de_file)) {
  message("[2/4] Reusing existing DE table")
  de <- read_tsv(de_file, show_col_types = FALSE)
} else {
  message("[2/4] Wilcoxon DE c7 vs Other")
  de <- FindMarkers(seu, ident.1 = "c7", ident.2 = "Other",
                    logfc.threshold = 0, min.pct = 0.01, test.use = "wilcox",
                    min.cells.group = 5)
  de$gene <- rownames(de)
  de$rank_stat <- -log10(pmax(de$p_val, 1e-300)) * sign(de$avg_log2FC)
  de <- de[order(-de$rank_stat), ]
  write_tsv(de, de_file)
}

message("[3/4] Load combined mouse GMT (Hallmark + PCa/CRPC)")
t2g_df <- read.gmt(GMT)
cat("     n pathways in GMT: ", length(unique(t2g_df$term)), "\n")

message("[4/4] clusterProfiler GSEA")
gene_list <- setNames(de$rank_stat, de$gene)
res <- GSEA(gene_list, TERM2GENE = t2g_df,
            pvalueCutoff = 1.0, minGSSize = 5, maxGSSize = 500,
            seed = TRUE, verbose = FALSE, eps = 0)
df <- as.data.frame(res@result) %>% arrange(desc(NES))
write_tsv(df, file.path(OUT, "gsea_c7_combined_result.tsv"))
message("     n sets: ", nrow(df),
        "; n FDR<0.05: ", sum(df$p.adjust < 0.05),
        "; n FDR<0.25: ", sum(df$p.adjust < 0.25))

# Print top 20 UP, top 15 DOWN
cat("\n=== c7 (COMBINED) — TOP 20 UP ===\n")
up <- df %>% arrange(desc(NES)) %>% head(20)
for (i in seq_len(nrow(up))) {
  marker <- if (up$p.adjust[i] < 0.05) " * " else if (up$p.adjust[i] < 0.25) " ~ " else "   "
  cat(sprintf("%s NES=%+.2f  padj=%.4f  n=%3d  %s\n",
              marker, up$NES[i], up$p.adjust[i], up$setSize[i], substr(up$ID[i], 1, 55)))
}
cat("\n=== c7 (COMBINED) — TOP 15 DOWN ===\n")
dn <- df %>% arrange(NES) %>% head(15)
for (i in seq_len(nrow(dn))) {
  marker <- if (dn$p.adjust[i] < 0.05) " * " else if (dn$p.adjust[i] < 0.25) " ~ " else "   "
  cat(sprintf("%s NES=%+.2f  padj=%.4f  n=%3d  %s\n",
              marker, dn$NES[i], dn$p.adjust[i], dn$setSize[i], substr(dn$ID[i], 1, 55)))
}
message("\nDONE")

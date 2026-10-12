#!/usr/bin/env Rscript
# GSEA of c7 Basal-like vs ALL OTHER 11 CARTA states — Wilcoxon-DE ranked,
# clusterProfiler on mouse Hallmark. Companion to the Fate 1/2/3 analyses.
suppressPackageStartupMessages({
  library(Seurat); library(zellkonverter); library(SingleCellExperiment)
  library(clusterProfiler); library(msigdbr); library(dplyr); library(readr); library(presto)
})
H5AD <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/data/npp53_merge/npp53_edited.h5ad"
OUT  <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12/biology_review"
STATES <- c('Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR','Luminal-Persistent',
            'Luminal-Squamous','Basal-like','Partial-EMT','Mes-Quiescent','Mes-Collagen',
            'Mes-Stemness','Mes-myCAF')
message("[1/4] read h5ad + subset"); sce <- readH5AD(H5AD)
counts_mat <- assay(sce, "logcounts")
seu <- CreateSeuratObject(counts = counts_mat, meta.data = as.data.frame(colData(sce)))
seu[["RNA"]]$data <- counts_mat
seu <- subset(seu, subset = putative_state %in% STATES)
seu$grp <- ifelse(seu$putative_state == "Basal-like", "c7", "Other")
Idents(seu) <- seu$grp
message("     kept ", ncol(seu), " cells (c7 = ", sum(seu$grp == "c7"), ")")
message("[2/4] Wilcoxon DE c7 vs Other")
de <- FindMarkers(seu, ident.1 = "c7", ident.2 = "Other",
                  logfc.threshold = 0, min.pct = 0.01, test.use = "wilcox",
                  min.cells.group = 5)
de$gene <- rownames(de)
de$rank_stat <- -log10(pmax(de$p_val, 1e-300)) * sign(de$avg_log2FC)
de <- de[order(-de$rank_stat), ]
write_tsv(de, file.path(OUT, "de_c7_vs_other_wilcox.tsv"))
message("[3/4] clusterProfiler GSEA (Hallmark)")
hallmark <- msigdbr(species = "Mus musculus", collection = "H")
t2g <- hallmark[, c("gs_name","gene_symbol")]; colnames(t2g) <- c("term","gene")
gene_list <- setNames(de$rank_stat, de$gene)
res <- GSEA(gene_list, TERM2GENE = t2g, pvalueCutoff = 1.0, minGSSize = 10, maxGSSize = 500,
            seed = TRUE, verbose = FALSE, eps = 0)
df <- as.data.frame(res@result) %>% arrange(desc(NES))
write_tsv(df, file.path(OUT, "gsea_c7_Hallmark_result.tsv"))
message("     n sets: ", nrow(df),
        "; n FDR<0.05: ", sum(df$p.adjust < 0.05),
        "; n FDR<0.25: ", sum(df$p.adjust < 0.25))
cat("\n=== c7 — TOP 15 UP ===\n")
up <- df %>% arrange(desc(NES)) %>% head(15)
for (i in seq_len(nrow(up))) {
  marker <- if (up$p.adjust[i] < 0.05) " * " else if (up$p.adjust[i] < 0.25) " ~ " else "   "
  cat(sprintf("%s NES=%+.2f  padj=%.4f  %s\n", marker, up$NES[i], up$p.adjust[i], substr(up$ID[i], 1, 60)))
}
cat("\n=== c7 — TOP 10 DOWN ===\n")
dn <- df %>% arrange(NES) %>% head(10)
for (i in seq_len(nrow(dn))) {
  marker <- if (dn$p.adjust[i] < 0.05) " * " else if (dn$p.adjust[i] < 0.25) " ~ " else "   "
  cat(sprintf("%s NES=%+.2f  padj=%.4f  %s\n", marker, dn$NES[i], dn$p.adjust[i], substr(dn$ID[i], 1, 60)))
}
message("\nDONE")

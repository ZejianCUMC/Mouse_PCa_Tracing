# GO BP + Hallmark enrichment per FINAL module (same settings as 02_go_enrich.R), for the module annotation table. Env: rMAST.
suppressMessages({library(clusterProfiler); library(org.Mm.eg.db); library(msigdbr); library(data.table)})
F <- fread("final_modules.csv"); U <- fread("universe_genes.csv")
H <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))[, .(gs_name, gene_symbol)]
out <- list(); hout <- list()
for (i in seq_len(nrow(F))) {
  g <- F$genotype[i]; m <- F$final_module[i]; genes <- strsplit(F$genes[i], " ")[[1]]; uni <- U[genotype == g]$gene
  e <- enrichGO(genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", universe = uni, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
  if (!is.null(e) && nrow(as.data.frame(e)) > 0) out[[i]] <- data.table(genotype = g, final_module = m, n_genes = length(genes), as.data.frame(e))
  h <- enricher(genes, TERM2GENE = H, universe = uni, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10)
  if (!is.null(h) && nrow(as.data.frame(h)) > 0) hout[[i]] <- data.table(genotype = g, final_module = m, n_genes = length(genes), as.data.frame(h))
  cat(g, m, length(genes), if (!is.null(e)) nrow(as.data.frame(e)) else 0, if (!is.null(h)) nrow(as.data.frame(h)) else 0, "\n")
}
fwrite(rbindlist(out), "go_bp_final_modules.csv"); fwrite(rbindlist(hout), "hallmark_final_modules.csv")

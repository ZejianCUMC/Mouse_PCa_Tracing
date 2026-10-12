# GO BP + Hallmark enrichment per Hotspot module (universe = genes in that genotype's meta-analysis). Env: rMAST. Deterministic.
suppressMessages({library(clusterProfiler); library(org.Mm.eg.db); library(msigdbr); library(data.table)})
G <- fread("hotspot_module_genes.csv"); U <- fread("universe_genes.csv")
H <- as.data.table(msigdbr(species = "Mus musculus", category = "H"))[, .(gs_name, gene_symbol)]
out <- list(); hout <- list()
for (g in unique(G$genotype)) for (m in unique(G[genotype == g]$hotspot_module)) {
  genes <- G[genotype == g & hotspot_module == m]$gene; uni <- U[genotype == g]$gene
  e <- enrichGO(genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", universe = uni, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
  if (!is.null(e) && nrow(as.data.frame(e)) > 0) out[[paste(g, m)]] <- data.table(genotype = g, hotspot_module = m, n_genes = length(genes), as.data.frame(e))
  h <- enricher(genes, TERM2GENE = H, universe = uni, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10)
  if (!is.null(h) && nrow(as.data.frame(h)) > 0) hout[[paste(g, m)]] <- data.table(genotype = g, hotspot_module = m, n_genes = length(genes), as.data.frame(h))
  cat(g, m, length(genes), if (!is.null(e)) nrow(as.data.frame(e)) else 0, if (!is.null(h)) nrow(as.data.frame(h)) else 0, "\n")
}
fwrite(rbindlist(out), "go_bp_per_module.csv"); fwrite(rbindlist(hout), "hallmark_per_module.csv")
cat(as.character(packageVersion("clusterProfiler")), as.character(packageVersion("org.Mm.eg.db")), as.character(packageVersion("msigdbr")), "\n")

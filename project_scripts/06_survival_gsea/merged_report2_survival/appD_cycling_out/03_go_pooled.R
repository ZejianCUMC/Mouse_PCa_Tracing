# GO BP + Hallmark enrichment of the top pooled candidates (mouse genes; universe = genes in the pooled meta-analysis), as fig45_function/04_go_final.R. Env: rMAST.
suppressMessages({library(clusterProfiler); library(org.Mm.eg.db); library(msigdbr); library(data.table); library(jsonlite)})
G <- fromJSON("pooled_top_mouse.json"); U <- fread(cmd="zcat ../../../metaNorm_20261011_all_cycout_bcr/meta_all.csv.gz")$gene
H <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))[, .(gs_name, gene_symbol)]; out <- list()
for (c in names(G)) {
  e <- enrichGO(G[[c]], OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", universe = U, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
  h <- enricher(G[[c]], TERM2GENE = H, universe = U, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10)
  for (x in list(list(e, "GO"), list(h, "Hallmark"))) if (!is.null(x[[1]]) && nrow(as.data.frame(x[[1]])) > 0) out[[length(out) + 1]] <- data.table(cand = c, db = x[[2]], as.data.frame(x[[1]])[, c("ID", "Description", "Count", "p.adjust")]) }
fwrite(rbindlist(out), "pooled_top_enrichment.csv")

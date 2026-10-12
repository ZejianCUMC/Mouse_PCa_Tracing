#!/usr/bin/env Rscript
# Export ggtree fan-layout tree coords (parent/node/x/y/isTip/label/angle) for JZ201/JZ205/JZ206
# Adapted from 23-ancestorOT/report/export_tree_coords.R
suppressPackageStartupMessages({library(ggtree); library(ape)})
BASE <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/videos_out"
NWK  <- file.path(BASE, "nwk_pt"); ANN  <- file.path(BASE, "anno")
OUT  <- file.path(BASE, "anim");   dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
for(t in c("JZ201","JZ202","JZ203","JZ205","JZ206")){
  tree <- ape::read.tree(file.path(NWK, sprintf("%s_cassiopeia_pt0.nwk", t)))
  ann  <- read.csv(file.path(ANN, sprintf("%s_rings.csv", t)), stringsAsFactors=FALSE)
  rownames(ann) <- ann$TipLabel
  common <- intersect(tree$tip.label, ann$TipLabel)
  tree   <- ape::keep.tip(tree, common)
  p <- ggtree(tree, layout="fan", open.angle=8)
  d <- p$data[, c("parent","node","x","y","isTip","label","angle")]
  d$CellType   <- NA_character_
  d$Pseudotime <- NA_real_
  tipmask <- d$isTip
  d$CellType[tipmask]   <- ann[d$label[tipmask], "CellType"]
  d$Pseudotime[tipmask] <- ann[d$label[tipmask], "Pseudotime"]
  write.csv(d, file.path(OUT, sprintf("%s_treecoords.csv", t)), row.names=FALSE)
  cat(sprintf("  %s: %d nodes (%d tips), x=[%.3f,%.3f]\n",
              t, nrow(d), sum(d$isTip), min(d$x, na.rm=TRUE), max(d$x, na.rm=TRUE)))
}
cat("DONE\n")

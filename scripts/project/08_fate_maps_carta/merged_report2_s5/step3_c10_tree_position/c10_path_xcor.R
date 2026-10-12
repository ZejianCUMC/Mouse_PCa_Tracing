#!/usr/bin/env Rscript
# =============================================================================
# Step 3b: where do c10 (cycling) cells sit on each tree? PATH phylogenetic cross-correlation
# =============================================================================
# Description : Same PATH calls and tree preparation as Npp53 evolution_features_v2/npp53_PATH_run.R
#               (unit branch lengths if missing, multi2di for polytomies, zero edges -> 1e-3,
#               catMat + one_node_tree_dist + xcor), on the merged 10-state labels of ALL tree-matched
#               cells of the 47 current trees. Reports the full z-score matrix per tumor and the c10 row.
# Input       : tumors.tsv, labels_all/{tumor}.tsv (from export_labels_all.py)
# Output      : path_z/{tumor}_zscore.csv, c10_path_summary.tsv
# Conda env   : rMAST (/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript)
# Key deps    : PATH, ape, Matrix
# Validated   : 2026-10-11
# =============================================================================
suppressMessages({ library(ape); library(Matrix); library(PATH) })
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
C10 <- "c10_Mes-Proliferation"; MIN_TIPS <- 30
tum <- read.delim(file.path(here, "tumors.tsv"), stringsAsFactors = FALSE)
dir.create(file.path(here, "path_z"), showWarnings = FALSE)
rows <- list()
for (i in seq_len(nrow(tum))) {
  t_id <- tum$tumor[i]
  tr <- read.tree(tum$tree[i])
  lab <- read.delim(file.path(here, "labels_all", paste0(t_id, ".tsv")), stringsAsFactors = FALSE)
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1.0, nrow(tr$edge))
  keep <- tr$tip.label %in% lab$leaf
  if (!all(keep)) tr <- drop.tip(tr, tr$tip.label[!keep])
  if (!is.binary(tr)) { tr <- multi2di(tr, random = FALSE) }
  tr$edge.length[tr$edge.length == 0] <- 1e-3
  st <- setNames(lab$new_state, lab$leaf)[tr$tip.label]
  present <- sort(unique(st))
  n_c10 <- sum(st == C10)
  if (length(tr$tip.label) < MIN_TIPS || !(C10 %in% present) || n_c10 < 2) {
    rows[[t_id]] <- data.frame(tumor = t_id, model = tum$model[i], treatment = tum$treatment[i],
                               n_tips = length(tr$tip.label), n_c10 = n_c10, c10_self_z = NA,
                               top_other_state = NA, top_other_z = NA, n_other_z_gt2 = NA, other_z = NA)
    next
  }
  X <- catMat(factor(st, levels = present)); W <- one_node_tree_dist(tr)
  xc <- xcor(X, W)
  z <- matrix(as.matrix(xc$Z.score), length(present), length(present), dimnames = list(present, present))
  write.csv(z, file.path(here, "path_z", paste0(t_id, "_zscore.csv")))
  zr <- z[C10, setdiff(present, C10)]
  zr <- sort(zr, decreasing = TRUE)
  rows[[t_id]] <- data.frame(tumor = t_id, model = tum$model[i], treatment = tum$treatment[i],
                             n_tips = length(tr$tip.label), n_c10 = n_c10, c10_self_z = round(z[C10, C10], 2),
                             top_other_state = names(zr)[1], top_other_z = round(zr[1], 2),
                             n_other_z_gt2 = sum(zr > 2),
                             other_z = paste(sprintf("%s:%.1f", sub("_.*", "", names(zr)), zr), collapse = " "))
  cat(sprintf("[%s] tips=%d c10=%d self_z=%.2f top=%s (%.2f)\n", t_id, length(tr$tip.label), n_c10, z[C10, C10], names(zr)[1], zr[1]))
}
out <- do.call(rbind, rows)
write.table(out, file.path(here, "c10_path_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

#!/usr/bin/env Rscript
# =============================================================================
# Npp53 PATH runs — EDIT-DISTANCE branch lengths (all 25 castration-cohort samples)
# =============================================================================
# Description : Same as npp53_PATH_run.R, but reads trees whose branch lengths
#               have been re-annotated by Camin-Sokal ancestral reconstruction
#               (each edge length = Hamming distance between parent/child
#               character states, zero pseudocounted to 1e-3).  See
#               annotate_trees_editdist.py for the tree annotation step.
# Input       : PATH_out/trees_editdist/{sample}_tree_greedy_editdist.nwk (25)
#               + evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv
# Output      : PATH_out/editdist/{sample}_PATH_{phylocorr,zscore,pvalue,rate_Q,
#                                                 transition_P,transition_Pt}.csv
#               + editdist/PATH_per_sample_scalars.csv
#               + editdist/PATH_heritability_long.csv
# Conda env   : rMAST (/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript)
# Key deps    : PATH, ape, Matrix, expm
# Validated   : 2026-09-08
# Parameters  : See DEFAULTS list below
# =============================================================================
suppressMessages({
  library(ape); library(Matrix); library(PATH)
})

DEFAULTS <- list(
  proj      = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing",
  tree_dir  = "1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/PATH_out/trees_editdist",
  tree_suf  = "_tree_greedy_editdist.nwk",
  scfp_csv  = "1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv",   # v10 = RECOMPUTED UMAP+scVI+diffmap on 18,241 no-cyc-no-ion cells
  outdir    = "1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/PATH_out/editdist",
  samples   = c("JZ136","JZ201","JZ202","JZ203","JZ204","JZ205","JZ206","JZ207",
                "JZ208","JZ209","JZ210","JZ211","JZ212","JZ213","JZ214","JZ224",
                "JZ225","JZ226","JZ227","JZ228","JZ229","JZ230","JZ231","JZ232","JZ234"),
  # putative_state levels observed in sc_fitness_plasticity (auto-detect at runtime)
  state_col = "putative_state"
)

run_sample <- function(sample, sc_fp, p, all_states) {
  cat(sprintf("[PATH] %s ...\n", sample))
  treef <- file.path(p$proj, p$tree_dir, paste0(sample, p$tree_suf))
  if (!file.exists(treef)) { cat("  no tree, skip\n"); return(NULL) }
  tr <- read.tree(treef)
  # Edit-distance trees carry Hamming-based branch lengths (with 1e-3 pseudocount
  # for zero-Hamming edges).  We still resolve polytomies for PATH.
  if (is.null(tr$edge.length)) stop(sprintf("expected edit-distance branch lengths on %s", treef))
  if (!is.binary(tr)) {
    cat("  polytomies found, breaking with multi2di\n")
    tr <- multi2di(tr, random = FALSE)
    # multi2di adds zero-length edges for resolved multichotomies — pseudocount
    tr$edge.length[tr$edge.length == 0] <- 1e-3
  }
  sub <- sc_fp[sc_fp$sample == sample, ]
  lab_map <- setNames(sub[[p$state_col]], sub$cellBC)
  fit_map <- setNames(sub$fitness,        sub$cellBC)
  has <- tr$tip.label %in% names(lab_map) & tr$tip.label %in% names(fit_map)
  if (!all(has)) tr <- drop.tip(tr, tr$tip.label[!has])
  states_raw  <- unname(lab_map[tr$tip.label])
  fitness_raw <- as.numeric(unname(fit_map[tr$tip.label]))
  present <- intersect(all_states, unique(states_raw))
  cat(sprintf("  tips=%d, states_present=%d\n", length(tr$tip.label), length(present)))
  if (length(tr$tip.label) < 30 || length(present) < 2) {
    cat("  too small, skip\n"); return(NULL)
  }
  tr$states <- factor(states_raw, levels = present)
  X  <- catMat(tr$states)
  W  <- one_node_tree_dist(tr)
  # discrete state phylo-corr
  xc <- xcor(X, W)
  nm <- present
  phy_cor <- matrix(as.matrix(xc$phy_cor), length(nm), length(nm), dimnames = list(nm, nm))
  zsc     <- matrix(as.matrix(xc$Z.score), length(nm), length(nm), dimnames = list(nm, nm))
  pval    <- matrix(as.matrix(xc$one.sided.pvalue), length(nm), length(nm), dimnames = list(nm, nm))
  # continuous fitness phylo-corr — single-column xcor
  fx <- matrix(fitness_raw, ncol = 1); colnames(fx) <- "fitness"
  fit_xc <- xcor(fx, W)
  fitness_phy_cor <- as.numeric(fit_xc$phy_cor[1,1])
  fitness_z       <- as.numeric(fit_xc$Z.score[1,1])
  fitness_p       <- as.numeric(fit_xc$one.sided.pvalue[1,1])
  # transitions
  Pinf <- tryCatch(PATH_inf(tree = tr, cell_states = "states", nstates = length(present)),
                   error = function(e) { message("   PATH_inf failed: ", conditionMessage(e)); NULL })
  # save
  outdir <- p$outdir
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  write.csv(phy_cor, file.path(outdir, paste0(sample, "_PATH_phylocorr.csv")))
  write.csv(zsc,     file.path(outdir, paste0(sample, "_PATH_zscore.csv")))
  write.csv(pval,    file.path(outdir, paste0(sample, "_PATH_pvalue.csv")))
  if (!is.null(Pinf)) {
    if (!is.null(Pinf$P))  { m <- as.matrix(Pinf$P);  dimnames(m) <- list(nm,nm); write.csv(m, file.path(outdir, paste0(sample, "_PATH_transition_P.csv"))) }
    if (!is.null(Pinf$Pt)) { m <- as.matrix(Pinf$Pt); dimnames(m) <- list(nm,nm); write.csv(m, file.path(outdir, paste0(sample, "_PATH_transition_Pt.csv"))) }
    if (!is.null(Pinf$Q))  { m <- as.matrix(Pinf$Q);  dimnames(m) <- list(nm,nm); write.csv(m, file.path(outdir, paste0(sample, "_PATH_rate_Q.csv"))) }
  }
  # per-mouse scalars
  herit_per_state <- diag(phy_cor)
  z_per_state     <- diag(zsc)
  n_z2 <- sum(z_per_state > 2, na.rm=TRUE)
  list(
    sample     = sample,
    n_tips     = length(tr$tip.label),
    n_states   = length(present),
    mean_herit = mean(herit_per_state, na.rm=TRUE),
    max_herit  = max(herit_per_state, na.rm=TRUE),
    mean_z_state = mean(z_per_state, na.rm=TRUE),
    max_z_state  = max(z_per_state, na.rm=TRUE),
    n_states_z2  = n_z2,
    fitness_phy_cor = fitness_phy_cor,
    fitness_z       = fitness_z,
    fitness_p       = fitness_p,
    herit_per_state = herit_per_state,
    z_per_state     = z_per_state
  )
}

main <- function(p = DEFAULTS) {
  outdir <- p$outdir
  dir.create(file.path(p$proj, outdir), showWarnings = FALSE, recursive = TRUE)
  p$outdir <- file.path(p$proj, outdir)
  sc_fp <- read.csv(file.path(p$proj, p$scfp_csv), stringsAsFactors = FALSE)
  # Zejian 2026-09-01: exclude Cycling + Ionocyte from PATH heritability analysis
  EXCLUDE <- c('Cycling', 'Ionocyte')
  n_before <- nrow(sc_fp)
  sc_fp <- sc_fp[!(sc_fp[[p$state_col]] %in% EXCLUDE), ]
  cat(sprintf("  excluded %s: %d -> %d cells retained\n", paste(EXCLUDE, collapse=", "), n_before, nrow(sc_fp)))
  all_states <- sort(unique(as.character(sc_fp[[p$state_col]])))
  cat(sprintf("[main] %d samples, %d discrete states in sc_fp\n", length(p$samples), length(all_states)))
  scalar_rows <- list(); herit_rows <- list()
  for (s in p$samples) {
    r <- run_sample(s, sc_fp, p, all_states)
    if (is.null(r)) next
    treatment <- unique(sc_fp$treatment[sc_fp$sample == s])[1]
    scalar_rows[[s]] <- data.frame(
      sample = s, treatment = treatment, n_tips = r$n_tips, n_states = r$n_states,
      mean_herit = r$mean_herit, max_herit = r$max_herit,
      mean_z_state = r$mean_z_state, max_z_state = r$max_z_state,
      n_states_z2 = r$n_states_z2,
      fitness_phy_cor = r$fitness_phy_cor, fitness_z = r$fitness_z, fitness_p = r$fitness_p
    )
    herit_rows[[s]] <- data.frame(
      sample = s, treatment = treatment,
      state = names(r$herit_per_state),
      heritability = as.numeric(r$herit_per_state),
      z = as.numeric(r$z_per_state)
    )
  }
  scalar_df <- do.call(rbind, scalar_rows)
  herit_df  <- do.call(rbind, herit_rows)
  write.csv(scalar_df, file.path(p$outdir, "PATH_per_sample_scalars.csv"), row.names = FALSE)
  write.csv(herit_df,  file.path(p$outdir, "PATH_heritability_long.csv"),  row.names = FALSE)
  cat("\n=== SCALAR SUMMARY ===\n")
  print(scalar_df)
  cat(sprintf("\nsaved: %s/PATH_per_sample_scalars.csv\n", p$outdir))
  cat(sprintf("saved: %s/PATH_heritability_long.csv\n", p$outdir))
}

if (sys.nframe() == 0) main()

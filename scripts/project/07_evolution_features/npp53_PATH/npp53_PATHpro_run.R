#!/usr/bin/env Rscript
# =============================================================================
# Npp53 PATHpro runs — per-state proliferation rates (γ) for all 25 samples
# =============================================================================
# Description : Compute PATH's PATHpro-inferred per-state proliferation rates
#               (γ) under a state-dependent birth-death CTMC on each sample's
#               tree. `total_prolif=1` normalization makes γ comparable across
#               samples with different state compositions (sum(γ * freq) = 1).
# Input       : cassiopeia_trees_v8_rebuilt/greedy/{sample}_tree_greedy.nwk (25)
#               + evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv
# Output      : evolution_features_v2/PATH_out/{sample}_PATHpro_{Q,P,gamma}.csv
#               + PATHpro_all_samples_gamma.csv (long: sample × state × γ)
#               + PATHpro_run_log.csv (per-sample status + timing)
# Conda env   : rMAST (/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript)
# Key deps    : PATH, ape
# Validated   : 2026-09-08
# =============================================================================
suppressMessages({ library(ape); library(PATH) })

DEFAULTS <- list(
  proj      = "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing",
  tree_dir  = "1-mouse_based/Npp53_tracing/lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy",
  tree_suf  = "_tree_greedy.nwk",
  scfp_csv  = "1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv",
  outdir    = "1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/PATH_out",
  samples   = c("JZ136","JZ201","JZ202","JZ203","JZ204","JZ205","JZ206","JZ207",
                "JZ208","JZ209","JZ210","JZ211","JZ212","JZ213","JZ214","JZ224",
                "JZ225","JZ226","JZ227","JZ228","JZ229","JZ230","JZ231","JZ232","JZ234"),
  state_col = "putative_state",
  # ≥20 cells per state per sample (CARTA convention; excludes very-rare noise states)
  min_cells_per_state = 20,
  # Per-sample soft wall-clock cap (CPU-secs); PATHpro is unconstrained beyond this.
  # setTimeLimit only fires between top-level R evaluations, so we log & report slow samples.
  cpu_cap_sec = 1200
)

run_sample <- function(sample, sc_fp, p) {
  t0 <- Sys.time()
  treef <- file.path(p$proj, p$tree_dir, paste0(sample, p$tree_suf))
  if (!file.exists(treef)) {
    return(list(sample=sample, status="no_tree_file",
                elapsed_sec=NA, tips=NA, states=NA, err=NA))
  }
  tr <- read.tree(treef)
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1.0, nrow(tr$edge))
  if (!is.binary(tr)) tr <- multi2di(tr, random = FALSE)
  tr$edge.length[tr$edge.length == 0] <- 1e-3
  sub <- sc_fp[sc_fp$sample == sample, ]
  lab_map <- setNames(sub[[p$state_col]], sub$cellBC)
  has <- tr$tip.label %in% names(lab_map)
  if (!all(has)) tr <- drop.tip(tr, tr$tip.label[!has])
  if (length(tr$tip.label) == 0) {
    return(list(sample=sample, status="no_tip_match",
                elapsed_sec=NA, tips=0, states=0, err=NA))
  }
  lab <- unname(lab_map[tr$tip.label])
  tab <- table(lab)
  keep_states <- names(tab)[tab >= p$min_cells_per_state]
  if (length(keep_states) < 2) {
    return(list(sample=sample, status=sprintf("too_few_states(only_%d_>=%d)",
                length(keep_states), p$min_cells_per_state),
                elapsed_sec=NA, tips=length(tr$tip.label), states=length(keep_states),
                err=NA))
  }
  keep_cells <- lab %in% keep_states
  tr  <- drop.tip(tr, tr$tip.label[!keep_cells])
  lab <- unname(lab_map[tr$tip.label])
  n_tips   <- length(tr$tip.label)
  n_states <- length(unique(lab))
  cat(sprintf("[PATHpro] %s tips=%d states=%d ...", sample, n_tips, n_states))
  setTimeLimit(cpu = p$cpu_cap_sec, elapsed = p$cpu_cap_sec, transient = TRUE)
  res <- tryCatch(
    PATHpro(tree=tr, cell_states=lab, depth=1, total_prolif=1),
    error = function(e) list(err = paste0(conditionMessage(e),
                                          " | trace: ",
                                          paste(sys.calls(), collapse=" << ")))
  )
  setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units="secs"))
  if (!is.null(res$err)) {
    cat(sprintf(" FAILED [%.1fs] %s\n", elapsed, substr(res$err, 1, 160)))
    return(list(sample=sample, status="PATHpro_error",
                elapsed_sec=round(elapsed,2), tips=n_tips, states=n_states,
                err=res$err))
  }
  cat(sprintf(" ok [%.1fs] γ.range=[%.3f,%.3f]\n",
              elapsed, min(res$gamma), max(res$gamma)))
  # attach state names to Q/P if PATHpro didn't (it uses default 1..K ordering
  # matching the levels of factor(lab))
  state_levels <- levels(factor(lab))
  if (is.null(dimnames(res$Q)))
    dimnames(res$Q) <- list(state_levels, state_levels)
  if (is.null(dimnames(res$P)))
    dimnames(res$P) <- list(state_levels, state_levels)
  gamma_named <- setNames(as.numeric(res$gamma), state_levels)
  # save per-sample matrices
  dir.create(p$outdir, showWarnings=FALSE, recursive=TRUE)
  write.csv(as.matrix(res$Q), file.path(p$outdir, paste0(sample, "_PATHpro_Q.csv")))
  write.csv(as.matrix(res$P), file.path(p$outdir, paste0(sample, "_PATHpro_P.csv")))
  write.csv(data.frame(state=state_levels, gamma=gamma_named, row.names=NULL),
            file.path(p$outdir, paste0(sample, "_PATHpro_gamma.csv")),
            row.names=FALSE)
  list(sample=sample, status="ok",
       elapsed_sec=round(elapsed, 2), tips=n_tips, states=n_states, err=NA,
       gamma_named=gamma_named)
}

main <- function(p = DEFAULTS) {
  p$outdir <- file.path(p$proj, p$outdir)
  dir.create(p$outdir, showWarnings=FALSE, recursive=TRUE)
  sc_fp <- read.csv(file.path(p$proj, p$scfp_csv), stringsAsFactors=FALSE)
  cat(sprintf("[main] %d samples, sc_fp = %d cells (v10)\n",
              length(p$samples), nrow(sc_fp)))
  log_rows <- list(); gamma_long <- list()
  for (s in p$samples) {
    r <- run_sample(s, sc_fp, p)
    treatment <- if (!is.null(sc_fp$treatment[sc_fp$sample == s])) {
      as.character(unique(sc_fp$treatment[sc_fp$sample == s]))[1]
    } else NA_character_
    log_rows[[s]] <- data.frame(
      sample=s, treatment=treatment, status=r$status,
      elapsed_sec=r$elapsed_sec, tips=r$tips, states=r$states,
      err=if (is.na(r$err)) "" else substr(r$err, 1, 400),
      stringsAsFactors=FALSE
    )
    if (r$status == "ok") {
      gamma_long[[s]] <- data.frame(
        sample=s, treatment=treatment,
        state=names(r$gamma_named),
        gamma=as.numeric(r$gamma_named),
        stringsAsFactors=FALSE
      )
    }
  }
  log_df   <- do.call(rbind, log_rows)
  gamma_df <- do.call(rbind, gamma_long)
  write.csv(log_df,   file.path(p$outdir, "PATHpro_run_log.csv"), row.names=FALSE)
  write.csv(gamma_df, file.path(p$outdir, "PATHpro_all_samples_gamma.csv"),
            row.names=FALSE)
  cat("\n=== RUN LOG ===\n"); print(log_df)
  cat(sprintf("\nsaved: %s/PATHpro_run_log.csv\n", p$outdir))
  cat(sprintf("saved: %s/PATHpro_all_samples_gamma.csv (%d rows)\n",
              p$outdir, nrow(gamma_df)))
}

if (sys.nframe() == 0) main()

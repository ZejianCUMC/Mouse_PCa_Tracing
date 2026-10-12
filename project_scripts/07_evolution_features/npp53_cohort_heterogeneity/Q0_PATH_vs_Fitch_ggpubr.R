#!/usr/bin/env Rscript
# Q0 · Fitch vs PATH — 2-panel comparison figure.
# Left panel: Fitch featurization (parent-child transitions, 144-d, Pearson corr dist).
# Right panel: PATH featurization (Moran's-I Z matrix, 144-d, Pearson corr dist).
# ggpubr violin + inner white boxplot (width=0.1), NO dots, similarity on y (higher = more similar).
# Wilcoxon bracket (number-only) between within / between per user preference.
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(ggsignif); library(dplyr); library(patchwork); library(tidyr)
})

M   <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
OUT <- file.path(M, "analysis/evolution_features_v2/lineage_6Q")

# ── Fitch pairs (already computed) ────────────────────────────────────────────
d_fitch <- as.matrix(read.csv(file.path(OUT, "Q0_pairwise_phylo_distance_corr.csv"),
                              row.names=1, check.names=FALSE))
s_fitch <- 1 - d_fitch   # Pearson r similarity

# ── PATH pairs (from long CSV — already tagged) ──────────────────────────────
path_long <- read.csv(file.path(OUT, "Q0_PATH_pair_long.csv"))
path_long <- path_long %>% filter(metric == "pearson_corr_dist")

# ── Fitch long-format ─────────────────────────────────────────────────────────
SCFP  <- read.csv(file.path(M, "analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv"))
sample_arm <- SCFP %>% group_by(sample) %>% summarise(treatment=first(treatment), .groups="drop")
ARM_SHORT  <- c("Intact"="Intact","Castration"="Cas","Castration + Enzalutamide"="Cas+Enz")
sample_arm$arm <- ARM_SHORT[sample_arm$treatment]

build_fitch_long <- function(){
  samples <- rownames(s_fitch); n <- length(samples)
  out <- list()
  for (i in seq_len(n-1)) for (j in (i+1):n){
    if (is.na(s_fitch[i,j])) next
    ai <- sample_arm$arm[match(samples[i], sample_arm$sample)]
    aj <- sample_arm$arm[match(samples[j], sample_arm$sample)]
    kind <- if (ai == aj) "within-arm" else "between-arm"
    out[[length(out)+1]] <- data.frame(method="Fitch", kind=kind, similarity=s_fitch[i,j])
  }
  bind_rows(out)
}
fitch_long <- build_fitch_long()
path_long_v <- data.frame(method="PATH", kind=path_long$kind, similarity=path_long$similarity)
long <- bind_rows(fitch_long, path_long_v)
long$kind <- factor(long$kind, levels=c("within-arm","between-arm"))

# ── PERMANOVA stats ───────────────────────────────────────────────────────────
fitch_stats <- read.csv(file.path(OUT, "Q0_permanova_pvalues.csv"))
path_stats  <- read.csv(file.path(OUT, "Q0_PATH_permanova_pvalues.csv"))
fitch_row <- fitch_stats[fitch_stats$metric == "pearson_corr_dist", ]
path_row  <- path_stats [path_stats$metric  == "pearson_corr_dist", ]

# ── per-method panel builder ──────────────────────────────────────────────────
make_panel <- function(method_name, stats_row){
  df <- long %>% filter(method == method_name)
  n_w <- sum(df$kind == "within-arm"); n_b <- sum(df$kind == "between-arm")
  # Use CSV mean-DISTANCE ratio (spec convention: <1 = within-arm tighter than between-arm)
  dist_ratio <- stats_row$within_over_between
  subtitle_txt <- sprintf(
    "PERMANOVA F=%.2f · p=%s\nmean(within dist) / mean(between dist) = %.2f  (<1 = within tighter)",
    stats_row$F, formatC(stats_row$permanova_p, format="g", digits=2), dist_ratio)
  title_txt <- if (method_name == "Fitch")
    "Fitch parent→child transition counts (144-d)"
  else
    "PATH Moran-I autocorrelation Z-matrix (144-d)"
  ggviolin(df, x="kind", y="similarity", fill="kind",
           palette=c("within-arm"="#4a7fbf","between-arm"="#c44e52"),
           trim=FALSE, alpha=0.45, width=0.80) +
    geom_boxplot(width=0.1, fill="white", color="grey20",
                 outlier.shape=NA, linewidth=0.5) +
    geom_signif(comparisons=list(c("within-arm","between-arm")), test="wilcox.test",
                map_signif_level=function(p) formatC(p, format="g", digits=2),
                textsize=3.6, tip_length=0.01, vjust=-0.15) +
    labs(title=title_txt, subtitle=subtitle_txt,
         y="phylogenetic similarity (Pearson r)", x=NULL) +
    scale_x_discrete(labels=c("within-arm"=sprintf("within-arm\nn=%d pairs", n_w),
                              "between-arm"=sprintf("between-arm\nn=%d pairs", n_b))) +
    theme_bw(base_size=11) +
    theme(legend.position="none",
          plot.title=element_text(size=11.5, face="bold", color="#0f2747"),
          plot.subtitle=element_text(size=10, color="#555"),
          axis.text.x=element_text(size=10),
          panel.grid.minor=element_blank())
}

fitch_panel <- make_panel("Fitch", fitch_row)
path_panel  <- make_panel("PATH",  path_row)

# ── concordance verdict ──────────────────────────────────────────────────────
fitch_sig <- fitch_row$permanova_p < 0.01 && fitch_row$within_over_between < 1.0
path_sig  <- path_row$permanova_p  < 0.01 && path_row$within_over_between  < 1.0
verdict <- if (fitch_sig == path_sig) "AGREE" else "DISAGREE"

combined <- (fitch_panel | path_panel) +
  plot_annotation(
    title="Q0 cohort-heterogeneity - Fitch vs PATH phylogenetic featurization",
    subtitle=sprintf(
      "Concordance verdict: %s (Fitch sig=%s, PATH sig=%s at p<0.01 & W/B<1)",
      verdict, fitch_sig, path_sig),
    theme=theme(
      plot.title=element_text(size=13, face="bold", color="#0f2747", hjust=0.5),
      plot.subtitle=element_text(size=11, color="#0f2747", hjust=0.5)))

ggsave(file.path(OUT, "Q0_PATH_comparison.pdf"), combined,
       width=11.5, height=5.8, device=cairo_pdf)
ggsave(file.path(OUT, "Q0_PATH_comparison.png"), combined,
       width=11.5, height=5.8, dpi=170)
message("[saved] Q0_PATH_comparison.png / .pdf")
message(sprintf("  Fitch: F=%.2f p=%s W/B=%.2f", fitch_row$F,
                formatC(fitch_row$permanova_p, format="g", digits=2),
                fitch_row$within_over_between))
message(sprintf("  PATH : F=%.2f p=%s W/B=%.2f", path_row$F,
                formatC(path_row$permanova_p, format="g", digits=2),
                path_row$within_over_between))
message(sprintf("  Verdict: %s", verdict))

#!/usr/bin/env Rscript
# Q0 · within-arm vs between-arm phylogenetic-distance — ggpubr boxplot version.
# No dots (per user), Wilcoxon p as bracket (number only).
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(ggsignif); library(dplyr); library(patchwork); library(tidyr)
})

M   <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
OUT <- file.path(M, "analysis/evolution_features_v2/lineage_6Q")

d_corr <- as.matrix(read.csv(file.path(OUT, "Q0_pairwise_phylo_distance_corr.csv"), row.names=1, check.names=FALSE))
d_bray <- as.matrix(read.csv(file.path(OUT, "Q0_pairwise_phylo_distance_bray.csv"), row.names=1, check.names=FALSE))
perm   <- read.csv(file.path(OUT, "Q0_permanova_pvalues.csv"))
# convert distances -> similarities (higher = more similar) per user request
s_corr <- 1 - d_corr    # Pearson r directly, range -1 .. +1
s_bray <- 1 - d_bray    # 1 - Bray-Curtis, range 0 .. 1

SCFP  <- read.csv(file.path(M, "analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv"))
sample_arm <- SCFP %>% group_by(sample) %>% summarise(treatment=first(treatment), .groups="drop")
ARM_SHORT <- c("Intact"="Intact","Castration"="Cas","Castration + Enzalutamide"="Cas+Enz")
sample_arm$arm <- ARM_SHORT[sample_arm$treatment]

# build long-format pair table (higher = more similar per user request)
build_long <- function(S, metric_name){
  samples <- rownames(S); n <- length(samples)
  out <- data.frame()
  for (i in seq_len(n-1)) for (j in (i+1):n){
    if (is.na(S[i,j])) next
    ai <- sample_arm$arm[match(samples[i], sample_arm$sample)]
    aj <- sample_arm$arm[match(samples[j], sample_arm$sample)]
    kind <- if (ai == aj) "within-arm" else "between-arm"
    out <- rbind(out, data.frame(metric=metric_name, kind=kind, similarity=S[i,j]))
  }
  out
}
long <- rbind(build_long(s_corr, "Pearson r (similarity)"),
              build_long(s_bray, "1 - Bray-Curtis (similarity)"))
long$kind <- factor(long$kind, levels=c("within-arm","between-arm"))

# per-metric panel — boxplot only, no dots
make_panel <- function(metric_name){
  df <- long %>% filter(metric == metric_name)
  perm_row <- perm[match(ifelse(grepl("Pearson", metric_name),"pearson_corr_dist","bray_curtis"), perm$metric), ]
  n_w <- sum(df$kind == "within-arm"); n_b <- sum(df$kind == "between-arm")
  # ratio computed on similarity (within/between; >1 means within more similar)
  med_w <- median(df$similarity[df$kind=="within-arm"], na.rm=TRUE)
  med_b <- median(df$similarity[df$kind=="between-arm"], na.rm=TRUE)
  sim_ratio <- med_w / med_b
  ggviolin(df, x="kind", y="similarity", fill="kind",
           palette=c("within-arm"="#4a7fbf","between-arm"="#c44e52"),
           trim=FALSE, alpha=0.45, width=0.80) +
    # inner boxplot (width=0.1, white fill) per user
    geom_boxplot(width=0.1, fill="white", color="grey20",
                 outlier.shape=NA, linewidth=0.5) +
    geom_signif(comparisons=list(c("within-arm","between-arm")), test="wilcox.test",
                map_signif_level=function(p) formatC(p, format="g", digits=2),
                textsize=3.6, tip_length=0.01, vjust=-0.15) +
    labs(title=metric_name,
         subtitle=sprintf("PERMANOVA F=%.2f · p=%s  |  within/between similarity ratio = %.2f",
                          perm_row$F, formatC(perm_row$permanova_p, format="g", digits=2),
                          sim_ratio),
         y="phylogenetic similarity", x=NULL) +
    scale_x_discrete(labels=c("within-arm"=sprintf("within-arm\nn=%d pairs", n_w),
                              "between-arm"=sprintf("between-arm\nn=%d pairs", n_b))) +
    theme_bw(base_size=11) +
    theme(legend.position="none",
          plot.title=element_text(size=12, face="bold", color="#0f2747"),
          plot.subtitle=element_text(size=10, color="#555"),
          axis.text.x=element_text(size=10),
          panel.grid.minor=element_blank())
}

panels <- lapply(c("Pearson r (similarity)","1 - Bray-Curtis (similarity)"), make_panel)
combined <- wrap_plots(panels, nrow=1) +
  plot_annotation(
    title="Cohort-heterogeneity — within-arm vs between-arm phylogenetic distance",
    theme=theme(plot.title=element_text(size=13, face="bold", color="#0f2747", hjust=0.5)))

ggsave(file.path(OUT, "Q0_within_vs_between.pdf"), combined,
       width=11.5, height=5.5, device=cairo_pdf)
ggsave(file.path(OUT, "Q0_within_vs_between.png"), combined,
       width=11.5, height=5.5, dpi=170)
message("[saved] Q0_within_vs_between.png (ggpubr boxplot, no dots)")

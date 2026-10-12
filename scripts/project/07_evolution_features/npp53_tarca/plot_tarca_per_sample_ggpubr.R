#!/usr/bin/env Rscript
# TarCA EffN per sample — ggpubr violin + inner boxplot (width=0.1, fill=white).
# One violin per sample; y = EffN across states; colored by treatment arm (same palette as Fig 6.1B).
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(dplyr); library(tidyr)
})

M <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
OUT <- file.path(M, "analysis/evolution_features_v2/tarca_npp53")
long <- read.csv(file.path(OUT, "tarca_all_samples_long.csv"), stringsAsFactors=FALSE)

# Arm palette — same as Fig 6.1B (Intact/Cas/Cas+Enz)
ARM_COLORS <- c("Intact"="#4a7fbf", "Cas"="#e0955f", "Cas+Enz"="#c44e52")
ARM_ORDER  <- c("Intact","Cas","Cas+Enz")

# Handle Inf (perfectly monophyletic clade, CR=0) — clip to visualization ceiling
LOG_CEIL <- 1e4
d <- long %>%
  filter(!is.na(EffN)) %>%
  mutate(EffN_viz = ifelse(is.infinite(EffN), LOG_CEIL, EffN),
         arm = factor(arm, levels=ARM_ORDER))

# Order samples on x by mean EffN (ascending) — computed on FINITE values only
# (drop Inf/NaN before mean so the order isn't dominated by monophyletic clades hitting the viz ceiling)
sample_order <- d %>%
  group_by(sample) %>%
  summarise(mean_effN = mean(EffN[is.finite(EffN)], na.rm=TRUE), .groups='drop') %>%
  arrange(mean_effN) %>%
  pull(sample)
d$sample <- factor(d$sample, levels=sample_order)

# Sample count per violin annotation
n_per_sample <- d %>% group_by(sample) %>% summarise(n=n(), .groups='drop')

g <- ggboxplot(d, x="sample", y="EffN_viz",
               color="arm", palette=ARM_COLORS, fill="white",
               outlier.shape=NA, size=0.75, width=0.7) +
  # visualization ceiling reference line (Inf clipping)
  geom_hline(yintercept=LOG_CEIL, color="#c44e52", linetype="dotted", linewidth=0.4, alpha=0.7) +
  scale_y_log10(labels=scales::comma,
                breaks=c(1, 10, 100, 1000, 10000),
                limits=c(0.5, LOG_CEIL*2)) +
  labs(y="Effective number of progenitors",
       x=NULL, color="Treatment") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        axis.title.y=element_text(size=12),
        legend.position="top",
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

ggsave(file.path(OUT, "tarca_per_sample_violin.pdf"), g,
       width=14, height=5.5, device=cairo_pdf)
ggsave(file.path(OUT, "tarca_per_sample_violin.png"), g,
       width=14, height=5.5, dpi=170)
message("[saved] tarca_per_sample_violin.png")

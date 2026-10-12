#!/usr/bin/env Rscript
# 2-class TarCA (Epi/Mes) — aggregate + per-arm boxplot (line-color by arm).
suppressPackageStartupMessages({library(ggplot2); library(ggpubr); library(dplyr); library(tidyr)})

M    <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
RES  <- file.path(M, "analysis/evolution_features_v2/tarca_npp53/results_2class")
OUT  <- file.path(M, "analysis/evolution_features_v2/tarca_npp53")
SCFP <- file.path(M, "analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv")

scfp <- read.csv(SCFP, stringsAsFactors=FALSE)
sample_arm <- scfp %>% group_by(sample) %>% summarise(treatment=first(treatment), .groups="drop")
arm_short <- c("Intact"="Intact","Castration"="Cas","Castration + Enzalutamide"="Cas+Enz")

# aggregate — one row per (sample, class)
rows <- list()
for (d in list.files(RES, full.names=TRUE)) {
  f <- file.path(d, "tarca_EffN.csv")
  if (!file.exists(f)) next
  s <- basename(d)
  df <- read.csv(f, stringsAsFactors=FALSE)
  df$sample <- s
  df$arm <- arm_short[sample_arm$treatment[match(s, sample_arm$sample)]]
  rows[[length(rows)+1]] <- df
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(OUT, "tarca_2class_all_samples_long.csv"), row.names=FALSE)

cat(sprintf("\n=== 2-class TarCA EffN per arm × class ===\n"))
for (cls in c("Epithelial","Mesenchymal")) {
  cat(sprintf("\n[%s]\n", cls))
  for (a in c("Intact","Cas","Cas+Enz")) {
    v <- tab$EffN[tab$TipAnn==cls & tab$arm==a]
    vf <- v[is.finite(v)]
    cat(sprintf("  %-8s  n=%d  n_finite=%d  n_inf=%d  n_nan=%d  median=%.1f  Q1-Q3=%.1f-%.1f\n",
                a, length(v), length(vf), sum(is.infinite(v)), sum(is.na(v)),
                if (length(vf)) median(vf) else NA,
                if (length(vf)) quantile(vf,0.25) else NA,
                if (length(vf)) quantile(vf,0.75) else NA))
  }
}

# --- boxplot: per class × arm, line color by arm, fill white ---
ARM_COLORS <- c("Intact"="#4a7fbf", "Cas"="#e0955f", "Cas+Enz"="#c44e52")
ARM_ORDER  <- c("Intact","Cas","Cas+Enz")

LOG_CEIL <- 1e4
d <- tab %>%
  filter(!is.na(EffN)) %>%
  mutate(EffN_viz = ifelse(is.infinite(EffN), LOG_CEIL, EffN),
         arm = factor(arm, levels=ARM_ORDER),
         TipAnn = factor(TipAnn, levels=c("Epithelial","Mesenchymal")))

g <- ggboxplot(d, x="arm", y="EffN_viz", color="arm",
               palette=ARM_COLORS, fill="white",
               outlier.shape=NA, size=0.85, width=0.6) +
  geom_jitter(aes(color=arm), width=0.15, size=1.8, alpha=0.85, show.legend=FALSE) +
  facet_wrap(~TipAnn, nrow=1) +
  geom_hline(yintercept=LOG_CEIL, color="#c44e52", linetype="dotted", linewidth=0.4, alpha=0.7) +
  scale_y_log10(labels=scales::comma,
                breaks=c(1, 10, 100, 1000, 10000),
                limits=c(0.5, LOG_CEIL*2)) +
  labs(y="Effective number of progenitors", x=NULL, color="Treatment") +
  theme_bw(base_size=12) +
  theme(strip.background=element_rect(fill="#f7f9fb"),
        strip.text=element_text(size=12, color="#0f2747", face="bold"),
        axis.title.y=element_text(size=12),
        legend.position="top",
        panel.grid.minor=element_blank())

ggsave(file.path(OUT, "tarca_2class_per_arm.pdf"), g, width=9, height=5, device=cairo_pdf)
ggsave(file.path(OUT, "tarca_2class_per_arm.png"), g, width=9, height=5, dpi=170)
message("[saved] tarca_2class_per_arm.png")

#!/usr/bin/env Rscript
# 2-class TarCA — paired Epi vs Mes for tumors that have BOTH valid values.
# Rule: Mes n_tips ≥ 20 AND finite EffN; Epi always finite (min tip count 49).
# n=12 tumors: 7 Intact + 2 Cas + 3 Cas+Enz.
# Test: paired Wilcoxon signed-rank (Epi vs Mes within tumor).
suppressPackageStartupMessages({library(ggplot2); library(ggpubr); library(ggsignif); library(dplyr); library(tidyr)})

M   <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
OUT <- file.path(M, "analysis/evolution_features_v2/tarca_npp53")
tab <- read.csv(file.path(OUT, "tarca_2class_all_samples_long.csv"), stringsAsFactors=FALSE)

ARM_COLORS <- c("Intact"="#4a7fbf", "Cas"="#e0955f", "Cas+Enz"="#c44e52")
MIN_TIPS_MES <- 20

# Keep only tumors where BOTH classes are valid
tab$class <- tab$TipAnn
paired_samples <- tab %>%
  filter(is.finite(EffN)) %>%
  filter(!(class=="Mesenchymal" & Total < MIN_TIPS_MES)) %>%
  group_by(sample) %>%
  filter(n() == 2) %>%      # keep only samples with both Epi and Mes rows surviving
  ungroup() %>%
  distinct(sample) %>%
  pull(sample)

d <- tab %>%
  filter(sample %in% paired_samples, is.finite(EffN)) %>%
  mutate(class = factor(class, levels=c("Epithelial","Mesenchymal")),
         arm   = factor(arm,   levels=names(ARM_COLORS)))

cat(sprintf("[paired set] n=%d tumors (Intact %d, Cas %d, Cas+Enz %d)\n",
            length(paired_samples),
            sum(d$arm[d$class=='Epithelial']=='Intact'),
            sum(d$arm[d$class=='Epithelial']=='Cas'),
            sum(d$arm[d$class=='Epithelial']=='Cas+Enz')))

# Paired Wilcoxon signed-rank on log-EffN
w <- wide <- d %>% select(sample, arm, class, EffN) %>%
  pivot_wider(names_from=class, values_from=EffN)
wt <- wilcox.test(log(w$Epithelial), log(w$Mesenchymal), paired=TRUE)
cat(sprintf("Paired Wilcoxon (log EffN, Epi vs Mes, n=%d): V=%.1f, p=%.4g\n",
            nrow(w), wt$statistic, wt$p.value))
cat(sprintf("Median EffN: Epi=%.2f  Mes=%.2f  |  log2FC (Mes/Epi) = %.2f\n",
            median(w$Epithelial), median(w$Mesenchymal),
            log2(median(w$Mesenchymal)/median(w$Epithelial))))

# --- Plot ---
p_label <- sprintf("paired Wilcoxon p = %.3g", wt$p.value)

g <- ggplot(d, aes(x=class, y=EffN)) +
  geom_boxplot(width=0.55, outlier.shape=NA,
               fill="white", color="#333", linewidth=0.7) +
  geom_line(aes(group=sample, color=arm),
            linewidth=0.5, alpha=0.6) +
  geom_point(aes(color=arm), size=2.8, alpha=0.9) +
  geom_signif(comparisons=list(c("Epithelial","Mesenchymal")),
              test="wilcox.test", test.args=list(paired=TRUE),
              map_signif_level=function(p) formatC(p, format="g", digits=2),
              tip_length=0.01, textsize=3.6, y_position=log10(180),
              size=0.5) +
  scale_color_manual(values=ARM_COLORS, name="Treatment") +
  scale_y_log10(labels=scales::comma,
                breaks=c(1, 3, 10, 30, 100),
                limits=c(0.9, 300)) +
  labs(x=NULL, y="Effective number of progenitors",
       subtitle=sprintf("Paired within tumor · n=%d tumors with Epi + Mes both n_tips ≥ %d",
                        length(paired_samples), MIN_TIPS_MES)) +
  theme_bw(base_size=12) +
  theme(axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12, face="bold"),
        legend.position="right",
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        plot.subtitle=element_text(size=9, color="grey35"))

ggsave(file.path(OUT, "tarca_2class_paired_epi_vs_mes.pdf"), g, width=6, height=5.5, device=cairo_pdf)
ggsave(file.path(OUT, "tarca_2class_paired_epi_vs_mes.png"), g, width=6, height=5.5, dpi=170)
message("[saved] tarca_2class_paired_epi_vs_mes.png")

# Print the paired table
cat("\n--- Paired EffN table ---\n")
print(as.data.frame(w %>% arrange(arm, sample) %>%
                    mutate(Epi = round(Epithelial,2),
                           Mes = round(Mesenchymal,2),
                           ratio_MesOverEpi = round(Mesenchymal/Epithelial,2)) %>%
                    select(sample, arm, Epi, Mes, ratio_MesOverEpi)),
      row.names=FALSE)

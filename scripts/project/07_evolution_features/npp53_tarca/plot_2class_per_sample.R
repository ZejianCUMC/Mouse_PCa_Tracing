#!/usr/bin/env Rscript
# 2-class TarCA (Epi/Mes) — per-sample paired view.
# 25 tumors on x, grouped by treatment arm.
# Each sample has 2 points connected by a line: Epithelial (open square) vs Mesenchymal (filled circle).
# Triangle glyph marks Inf-clipped values (perfectly monophyletic clade — CR estimated as 0).
suppressPackageStartupMessages({library(ggplot2); library(ggpubr); library(dplyr); library(tidyr)})

M   <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
OUT <- file.path(M, "analysis/evolution_features_v2/tarca_npp53")
tab <- read.csv(file.path(OUT, "tarca_2class_all_samples_long.csv"), stringsAsFactors=FALSE)

ARM_COLORS <- c("Intact"="#4a7fbf", "Cas"="#e0955f", "Cas+Enz"="#c44e52")
ARM_ORDER  <- c("Intact","Cas","Cas+Enz")
LOG_CEIL   <- 1e4

d <- tab %>%
  mutate(EffN_viz = ifelse(is.infinite(EffN), LOG_CEIL, EffN),
         censored = is.infinite(EffN),
         arm = factor(arm, levels=ARM_ORDER),
         TipAnn = factor(TipAnn, levels=c("Epithelial","Mesenchymal")))

# sample order: by arm, then by Epithelial EffN within arm (Inf → last)
epi_order <- d %>%
  filter(TipAnn=="Epithelial") %>%
  mutate(sort_key = ifelse(is.infinite(EffN), Inf, EffN)) %>%
  arrange(arm, sort_key) %>%
  pull(sample)
# tumors that have no Epi row (shouldn't happen) — append at end grouped by arm
all_samps <- unique(d$sample)
missing   <- setdiff(all_samps, epi_order)
sample_order <- c(epi_order, missing)
d$sample <- factor(d$sample, levels=sample_order)

# arm banding (background rectangles per arm block)
arm_of <- d %>% select(sample, arm) %>% distinct() %>% arrange(sample)
band <- arm_of %>%
  mutate(x = as.integer(sample)) %>%
  group_by(arm) %>%
  summarise(xmin=min(x)-0.5, xmax=max(x)+0.5, .groups="drop") %>%
  mutate(fill = ARM_COLORS[as.character(arm)])

# paired line (Epi ↔ Mes) — one line per sample; drop where either side is NA
paired <- d %>%
  select(sample, arm, TipAnn, EffN_viz, EffN) %>%
  filter(!is.na(EffN)) %>%
  group_by(sample) %>%
  filter(n() == 2) %>%
  ungroup()

# NaN Mes samples (only 1 tip in Mes) — no Mes point drawn
n_mes_missing <- sum(d$TipAnn=="Mesenchymal" & is.na(d$EffN))
n_mes_inf     <- sum(d$TipAnn=="Mesenchymal" & is.infinite(d$EffN))
n_epi_inf     <- sum(d$TipAnn=="Epithelial"  & is.infinite(d$EffN))
message(sprintf("[info] Mes: %d NaN (single-tip), %d Inf (monophyletic)  |  Epi: %d Inf",
                n_mes_missing, n_mes_inf, n_epi_inf))

d_draw <- d %>% filter(!is.na(EffN))

g <- ggplot() +
  # arm banding
  geom_rect(data=band,
            aes(xmin=xmin, xmax=xmax, ymin=0.5, ymax=LOG_CEIL*2, fill=fill),
            alpha=0.06, inherit.aes=FALSE, show.legend=FALSE) +
  scale_fill_identity() +
  # paired connector line, one per sample, colored by arm
  geom_line(data=paired,
            aes(x=sample, y=EffN_viz, group=sample, color=arm),
            linewidth=0.6, alpha=0.55, show.legend=FALSE) +
  # finite points — Epi = open square, Mes = filled circle
  geom_point(data=d_draw %>% filter(!censored),
             aes(x=sample, y=EffN_viz, color=arm, shape=TipAnn, fill=TipAnn),
             size=3.2, stroke=1.0) +
  # censored (Inf) points — triangle at the ceiling
  geom_point(data=d_draw %>% filter(censored),
             aes(x=sample, y=EffN_viz, color=arm),
             shape=17, size=3.2, alpha=0.9) +
  # ceiling reference
  geom_hline(yintercept=LOG_CEIL, color="#c44e52", linetype="dotted", linewidth=0.4, alpha=0.7) +
  scale_shape_manual(values=c("Epithelial"=22, "Mesenchymal"=21),
                     name="Cell state") +
  scale_fill_manual(values=c("Epithelial"="white", "Mesenchymal"="grey40"),
                    name="Cell state", guide="none") +
  scale_color_manual(values=ARM_COLORS, name="Treatment") +
  scale_y_log10(labels=scales::comma,
                breaks=c(1, 10, 100, 1000, 10000),
                limits=c(0.5, LOG_CEIL*2)) +
  labs(y="Effective number of progenitors", x=NULL,
       caption=sprintf("▲ = Inf-clipped (perfectly monophyletic clade; CR≈0). %d Mes samples with only 1 tip omitted.",
                       n_mes_missing)) +
  guides(color=guide_legend(order=1, override.aes=list(shape=16, size=3)),
         shape=guide_legend(order=2, override.aes=list(color="black", size=3))) +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        axis.title.y=element_text(size=12),
        legend.position="top",
        legend.box="horizontal",
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        plot.caption=element_text(size=8, color="grey40", hjust=0))

ggsave(file.path(OUT, "tarca_2class_per_sample.pdf"), g, width=14, height=5.5, device=cairo_pdf)
ggsave(file.path(OUT, "tarca_2class_per_sample.png"), g, width=14, height=5.5, dpi=170)
message("[saved] tarca_2class_per_sample.png")

# Also print the full per-sample table sorted by arm then Epi EffN, for the record
d_wide <- d %>%
  select(sample, arm, TipAnn, EffN) %>%
  mutate(EffN_disp = ifelse(is.infinite(EffN), "Inf",
                            ifelse(is.na(EffN), "NaN", sprintf("%.2f", EffN)))) %>%
  select(-EffN) %>%
  pivot_wider(names_from=TipAnn, values_from=EffN_disp) %>%
  arrange(match(arm, ARM_ORDER), sample)
print(as.data.frame(d_wide), row.names=FALSE)
write.csv(d_wide, file.path(OUT, "tarca_2class_per_sample_table.csv"), row.names=FALSE)

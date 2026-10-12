#!/usr/bin/env Rscript
# 2-class TarCA — per-class filtering (Mes n_tips ≥ 20; Epi all 25 kept).
# Two plots:
#   (A) per-sample paired view — all 25 Epi + 12 Mes; missing Mes shown as small grey "n.d." marker
#   (B) per-arm boxplot — same filter, arm-level stats
suppressPackageStartupMessages({library(ggplot2); library(ggpubr); library(dplyr); library(tidyr)})

M   <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
OUT <- file.path(M, "analysis/evolution_features_v2/tarca_npp53")
tab <- read.csv(file.path(OUT, "tarca_2class_all_samples_long.csv"), stringsAsFactors=FALSE)

ARM_COLORS <- c("Intact"="#4a7fbf", "Cas"="#e0955f", "Cas+Enz"="#c44e52")
ARM_ORDER  <- c("Intact","Cas","Cas+Enz")
MIN_TIPS_MES <- 20    # per-class threshold for Mes only

# --- per-class filtering ---
# Add sc_fitness_plasticity cell counts (n_tips) — but Total in tarca csv already = n_tips
tab$n_tips <- tab$Total
tab$class  <- tab$TipAnn
# Filter Mes rows failing threshold (drop from analysis, but their Epi twin stays)
tab_kept <- tab %>%
  filter(!(class == "Mesenchymal" & (n_tips < MIN_TIPS_MES | is.na(n_tips))))
# Also drop residual Inf (JZ210 Mes n_tips=12 already gone; safety catch)
tab_kept <- tab_kept %>% filter(is.finite(EffN))
cat(sprintf("[filter] Epi kept: %d/25 tumors; Mes kept: %d/25 tumors (n_tips ≥ %d, all finite)\n",
            sum(tab_kept$class=='Epithelial'), sum(tab_kept$class=='Mesenchymal'), MIN_TIPS_MES))

# Include tumors with 0 Mes cells (JZ226, JZ232) — they weren't in tarca output; add Epi rows only.
# (Already present in tarca since Epi is always reported.)

# ── PLOT A — per-sample paired view ─────────────────────────────────
d <- tab_kept %>%
  mutate(arm = factor(arm, levels=ARM_ORDER),
         class = factor(class, levels=c("Epithelial","Mesenchymal")))
# sample order: by arm, then by Epi EffN
epi_order <- d %>% filter(class=="Epithelial") %>% arrange(arm, EffN) %>% pull(sample)
missing   <- setdiff(unique(tab$sample), epi_order)
sample_order <- c(epi_order, missing)
d$sample <- factor(d$sample, levels=sample_order)

# arm banding
arm_of <- d %>% select(sample, arm) %>% distinct() %>% arrange(sample)
band <- arm_of %>%
  mutate(x=as.integer(sample)) %>%
  group_by(arm) %>%
  summarise(xmin=min(x)-0.5, xmax=max(x)+0.5, .groups="drop") %>%
  mutate(fill=ARM_COLORS[as.character(arm)])

# paired connector (only samples with BOTH Epi and Mes kept)
paired <- d %>% group_by(sample) %>% filter(n()==2) %>% ungroup()

# "n.d." marker for Mes-missing samples — placed at y=0.55 (near bottom)
mes_missing <- tab %>%
  filter(class=="Mesenchymal") %>%
  mutate(missing = is.na(EffN) | !is.finite(EffN) | n_tips < MIN_TIPS_MES) %>%
  filter(missing) %>% pull(sample)
# also add the 2 zero-Mes samples that had no Mes row at all
zero_mes <- setdiff(unique(tab$sample[tab$class=="Epithelial"]),
                    unique(tab$sample[tab$class=="Mesenchymal"]))
mes_missing_all <- unique(c(mes_missing, zero_mes))
nd <- data.frame(sample = factor(mes_missing_all, levels=sample_order),
                 arm    = arm_of$arm[match(mes_missing_all, arm_of$sample)],
                 EffN   = 0.6)

g <- ggplot() +
  geom_rect(data=band,
            aes(xmin=xmin, xmax=xmax, ymin=0.4, ymax=2e4, fill=fill),
            alpha=0.06, inherit.aes=FALSE, show.legend=FALSE) +
  scale_fill_identity() +
  geom_line(data=paired,
            aes(x=sample, y=EffN, group=sample, color=arm),
            linewidth=0.6, alpha=0.55, show.legend=FALSE) +
  geom_point(data=d,
             aes(x=sample, y=EffN, color=arm, shape=class, fill=class),
             size=3.2, stroke=1.0) +
  geom_point(data=nd,
             aes(x=sample, y=EffN),
             shape=4, color="grey55", size=2.4, stroke=0.9) +
  annotate("text", x=1, y=0.44, label="× = Mes n<20 (not shown)",
           hjust=0, color="grey45", size=3) +
  scale_shape_manual(values=c("Epithelial"=22, "Mesenchymal"=21), name="Cell state") +
  scale_fill_manual(values=c("Epithelial"="white", "Mesenchymal"="grey40"),
                    name="Cell state", guide="none") +
  scale_color_manual(values=ARM_COLORS, name="Treatment") +
  scale_y_log10(labels=scales::comma,
                breaks=c(1, 10, 100, 1000),
                limits=c(0.4, 2e2)) +
  labs(y="Effective number of progenitors", x=NULL,
       subtitle=sprintf("Per-class filter: Epi all 25 tumors; Mes n_tips ≥ %d (12 tumors kept)", MIN_TIPS_MES)) +
  guides(color=guide_legend(order=1, override.aes=list(shape=16, size=3)),
         shape=guide_legend(order=2, override.aes=list(color="black", size=3))) +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        axis.title.y=element_text(size=12),
        legend.position="top", legend.box="horizontal",
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        plot.subtitle=element_text(size=9, color="grey35"))

ggsave(file.path(OUT, "tarca_2class_per_sample_filtered.pdf"), g, width=14, height=5.5, device=cairo_pdf)
ggsave(file.path(OUT, "tarca_2class_per_sample_filtered.png"), g, width=14, height=5.5, dpi=170)
message("[saved] tarca_2class_per_sample_filtered.png")

# ── PLOT B — per-arm boxplot (Epi all 25 | Mes n≥20) ────────────────
h <- ggboxplot(d, x="arm", y="EffN", color="arm",
               palette=ARM_COLORS, fill="white",
               outlier.shape=NA, size=0.85, width=0.6) +
  geom_jitter(aes(color=arm), width=0.15, size=1.8, alpha=0.85, show.legend=FALSE) +
  facet_wrap(~class, nrow=1) +
  scale_y_log10(labels=scales::comma,
                breaks=c(1, 10, 100, 1000),
                limits=c(0.5, 200)) +
  labs(y="Effective number of progenitors", x=NULL, color="Treatment",
       subtitle=sprintf("Per-class filter: Epi all 25 tumors; Mes n_tips ≥ %d (12 tumors kept)", MIN_TIPS_MES)) +
  theme_bw(base_size=12) +
  theme(strip.background=element_rect(fill="#f7f9fb"),
        strip.text=element_text(size=12, color="#0f2747", face="bold"),
        axis.title.y=element_text(size=12),
        legend.position="top",
        panel.grid.minor=element_blank(),
        plot.subtitle=element_text(size=9, color="grey35"))

ggsave(file.path(OUT, "tarca_2class_per_arm_filtered.pdf"), h, width=9, height=5, device=cairo_pdf)
ggsave(file.path(OUT, "tarca_2class_per_arm_filtered.png"), h, width=9, height=5, dpi=170)
message("[saved] tarca_2class_per_arm_filtered.png")

# ── per-arm summary numbers (filtered) ──
cat("\n=== Per-arm EffN (filtered) ===\n")
for (cls in c("Epithelial","Mesenchymal")) {
  cat(sprintf("\n[%s]\n", cls))
  for (a in ARM_ORDER) {
    v <- d$EffN[d$class==cls & d$arm==a]
    cat(sprintf("  %-8s  n=%d  median=%.2f  IQR=%.2f-%.2f\n",
                a, length(v),
                if(length(v)) median(v) else NA,
                if(length(v)) quantile(v,0.25) else NA,
                if(length(v)) quantile(v,0.75) else NA))
  }
}

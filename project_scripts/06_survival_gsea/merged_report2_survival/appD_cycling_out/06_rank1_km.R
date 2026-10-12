#!/usr/bin/env Rscript
# Appendix D (cycling cells excluded): KM curves (median split, all five cohorts) for the rank-1 merged candidate, same style and scoring as ../01_km_all_cohorts.R. Env: rMAST.
suppressPackageStartupMessages({library(singscore); library(survival); library(survminer); library(jsonlite); library(ggplot2); library(cowplot); library(data.table)})
O <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); S <- dirname(dirname(dirname(O)))
C <- readRDS(file.path(S, "cohorts_corrected.rds")); T <- fread(file.path(O, "pooled_top.csv")); G <- fromJSON(file.path(O, "pooled_top_human.json"))[[paste0("all|", T$cand[1])]]
CO <- c("SU2C-mCRPC-union"="SU2C mCRPC · overall survival", "WCDT-JCO"="WCDT mCRPC · overall survival", "TCGA-PRAD"="TCGA-PRAD primary · progression-free interval",
        "DKFZ-primary"="DKFZ primary · biochemical recurrence", "MSKCC-primary"="MSKCC primary · biochemical recurrence")
kmp <- function(n) { co <- C[[n]]; rk <- rankGenes(co$mat); s <- setNames(simpleScore(rk, upSet=intersect(G, rownames(rk)), centerScore=FALSE)$TotalScore, colnames(rk)); id <- co$map[names(s)]; ok <- !is.na(id)
  x <- tapply(s[ok], id[ok], mean); d <- merge(data.frame(id=names(x), score=as.numeric(x)), co$cli, by="id"); d <- d[complete.cases(d) & d$time > 0, ]; d$time_mo <- d$time / 30.44
  d$grp <- factor(ifelse(d$score > median(d$score), "high", "low"), levels=c("low", "high")); p <- 1 - pchisq(survdiff(Surv(time, event) ~ grp, data=d)$chisq, 1); HR <- exp(coef(coxph(Surv(time, event) ~ grp, data=d)))[1]
  cap <- min(120, ceiling(max(d$time_mo) / 12) * 12)
  pp <- ggsurvplot(survfit(Surv(time_mo, event) ~ grp, data=d), data=d, pval=sprintf("p = %.2g", p), pval.size=4.5, pval.coord=c(3, .12), conf.int=TRUE, risk.table=TRUE, risk.table.height=.3,
                   risk.table.title="At risk", risk.table.y.text=FALSE, palette=c("#4a7fbf", "#c44e52"), legend.title="score", legend.labs=c("low", "high"), xlab="Time (months)", ylab="Event-free probability",
                   xlim=c(0, cap), break.time.by=if (cap > 72) 24 else 12, title=sprintf("%s · median%s\n(HR=%.2f · p=%.2g · n=%d, %d events)", CO[[n]], if (p < .05 && HR > 1) " *" else "", HR, p, nrow(d), sum(d$event)),
                   ggtheme=theme_bw(base_size=10) + theme(plot.title=element_text(hjust=0, size=10, face="bold", color="#0f2747"), legend.position="top", panel.grid.minor=element_blank()),
                   font.legend=c(9, "plain"), font.tickslab=c(9), font.x=c(10), font.y=c(10))
  plot_grid(pp$plot, pp$table + theme(plot.title=element_text(size=9)), ncol=1, rel_heights=c(.72, .28)) }
pl <- lapply(names(CO), kmp)
out <- plot_grid(ggdraw() + draw_label(sprintf("Rank-1 merged candidate (no cycling cells) (%d human genes) · median split of the module score", length(G)), x=.01, hjust=0, size=13, fontface="bold", colour="#0f2747"), plot_grid(plotlist=pl, ncol=3), ncol=1, rel_heights=c(.05, .95))
ggsave(file.path(O, "figC3_rank1_km.png"), out, width=16, height=11.5, dpi=105, bg="white")

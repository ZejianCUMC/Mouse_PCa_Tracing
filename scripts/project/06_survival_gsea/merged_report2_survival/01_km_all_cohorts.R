#!/usr/bin/env Rscript
# Report-2 §4.8 rework (Zejian 2026-10-11): median split as the default, more prostate-cancer cohorts, one WCDT cohort only, mean split
# moved to Appendix A. Same scoring and KM style as report2_figures_20261010/03_km_forest.R (singscore simpleScore uncentred on
# cohorts_corrected.rds, patient mean; survminer KM with CI band and number-at-risk table; HR from Cox on the groups, log-rank p).
# Cohorts: SU2C-mCRPC-union (OS), WCDT-JCO (OS, EU-updated), TCGA-PRAD (PFI), DKFZ-primary (BCR), MSKCC-primary (BCR).
# Outputs: fig_km_<module>_median.png, fig_km_<module>_mean.png, km_cohorts.csv (all splits + continuous per SD). Env: rMAST.
suppressPackageStartupMessages({library(singscore); library(survival); library(survminer); library(jsonlite); library(ggplot2); library(cowplot); library(data.table)})
O <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); S <- dirname(dirname(O))
C <- readRDS(file.path(S, "cohorts_corrected.rds"))
MOD <- list("Npp53 risk module"=unique(fromJSON(file.path(S, "metaNorm_20261010_Npp53_bcr/candidates_human.json"))[["Q03389"]]),
            "Npp53Rb1 risk module"=unique(fromJSON(file.path(S, "metaNorm_20261010_Npp53Rb1_bcr/candidates_human.json"))[["Q03123"]]))
FN <- c("Npp53 risk module"="Npp53", "Npp53Rb1 risk module"="Npp53Rb1")
CO <- c("SU2C-mCRPC-union"="SU2C mCRPC · overall survival", "WCDT-JCO"="WCDT mCRPC · overall survival", "TCGA-PRAD"="TCGA-PRAD primary · progression-free interval",
        "DKFZ-primary"="DKFZ primary · biochemical recurrence", "MSKCC-primary"="MSKCC primary · biochemical recurrence")
P <- lapply(setNames(names(CO), names(CO)), function(n) list(rk=rankGenes(C[[n]]$mat), map=C[[n]]$map, cli=C[[n]]$cli))
dat <- function(n, g) { p <- P[[n]]; s <- setNames(simpleScore(p$rk, upSet=intersect(g, rownames(p$rk)), centerScore=FALSE)$TotalScore, colnames(p$rk)); id <- p$map[names(s)]; ok <- !is.na(id)
  x <- tapply(s[ok], id[ok], mean); d <- merge(data.frame(id=names(x), score=as.numeric(x)), p$cli, by="id"); d <- d[complete.cases(d) & d$time > 0, ]; d$time_mo <- d$time / 30.44; d }
grp <- function(x, sp) { q3 <- quantile(x, c(1/3, 2/3)); q4 <- quantile(x, c(.25, .75))
  switch(sp, median=ifelse(x > median(x), "high", "low"), mean=ifelse(x > mean(x), "high", "low"), tertile=ifelse(x <= q3[1], "low", ifelse(x > q3[2], "high", NA)),
         q1q4=ifelse(x <= q4[1], "low", ifelse(x >= q4[2], "high", NA))) }
SPN <- c(median="median", mean="mean", tertile="T3 vs T1", q1q4="Q4 vs Q1")
kmp <- function(d, sp, ttl0) { d$grp <- factor(grp(d$score, sp), levels=c("low", "high")); d <- d[!is.na(d$grp), ]
  p <- 1 - pchisq(survdiff(Surv(time, event) ~ grp, data=d)$chisq, 1); HR <- exp(coef(coxph(Surv(time, event) ~ grp, data=d)))[1]
  ttl <- sprintf("%s · %s%s\n(HR=%.2f · p=%.2g · n=%d, %d events)", ttl0, SPN[sp], if (p < .05 && HR > 1) " *" else "", HR, p, nrow(d), sum(d$event))
  cap <- min(120, ceiling(max(d$time_mo) / 12) * 12)
  pp <- ggsurvplot(survfit(Surv(time_mo, event) ~ grp, data=d), data=d, pval=sprintf("p = %.2g", p), pval.size=4.5, pval.coord=c(3, .12), conf.int=TRUE, risk.table=TRUE,
                   risk.table.height=.3, risk.table.title="At risk", risk.table.y.text=FALSE, palette=c("#4a7fbf", "#c44e52"), legend.title="score", legend.labs=c("low", "high"),
                   xlab="Time (months)", ylab="Event-free probability", xlim=c(0, cap), break.time.by=if (cap > 72) 24 else 12,
                   ggtheme=theme_bw(base_size=10) + theme(plot.title=element_text(hjust=0, size=10, face="bold", color="#0f2747"), legend.position="top", panel.grid.minor=element_blank()),
                   title=ttl, font.legend=c(9, "plain"), font.tickslab=c(9), font.x=c(10), font.y=c(10))
  plot_grid(pp$plot, pp$table + theme(plot.title=element_text(size=9)), ncol=1, rel_heights=c(.72, .28)) }
res <- list()
for (m in names(MOD)) {
  D <- lapply(names(CO), function(n) dat(n, MOD[[m]])); names(D) <- names(CO)
  for (sp in c("median", "mean")) {
    pl <- lapply(names(CO), function(n) kmp(D[[n]], sp, CO[[n]]))
    ttl <- sprintf("%s (%d human genes) · %s split of the module score", m, length(MOD[[m]]), sp)
    out <- plot_grid(ggdraw() + draw_label(ttl, x=.01, hjust=0, size=13, fontface="bold", colour="#0f2747"), plot_grid(plotlist=pl, ncol=3), ncol=1, rel_heights=c(.05, .95))
    ggsave(file.path(O, sprintf("fig_km_%s_%s.png", FN[[m]], sp)), out, width=16, height=11.5, dpi=105, bg="white") }
  for (n in names(CO)) { d <- D[[n]]
    for (sp in names(SPN)) { dd <- d; dd$grp <- factor(grp(dd$score, sp), levels=c("low", "high")); dd <- dd[!is.na(dd$grp), ]
      p <- 1 - pchisq(survdiff(Surv(time, event) ~ grp, data=dd)$chisq, 1); f <- coxph(Surv(time, event) ~ grp, data=dd)
      res[[length(res) + 1]] <- data.table(module=m, cohort=n, split=sp, HR=exp(coef(f)), p=p, n=nrow(dd), events=sum(dd$event)) }
    f <- coxph(Surv(time, event) ~ scale(score), data=d); res[[length(res) + 1]] <- data.table(module=m, cohort=n, split="cont", HR=exp(coef(f)), p=summary(f)$coefficients[1, 5], n=nrow(d), events=sum(d$event)) } }
R <- rbindlist(res); fwrite(R, file.path(O, "km_cohorts.csv")); print(R[, .(module=substr(module, 1, 9), cohort, split, HR=round(HR, 2), p=signif(p, 2), n, events)])

#!/usr/bin/env Rscript
# opt5 step 26 — median-split KM curves, identical layout for every module (Zejian 2026-10-10: same split and same cohorts for Npp53
# and Npp53Rb1). Usage: 26_km_median.R <out_png> <label>=<arm>:<cand> [...]
# Scoring exactly as km_scan.R: singscore simpleScore (uncentred) per sample on cohorts_corrected.rds, patient mean, high = above the
# cohort median; HR from Cox on the groups, p from the log-rank test. Panels per module: SU2C-mCRPC-union (OS), WCDT-JCO (OS).
# Env: rMAST.
suppressPackageStartupMessages({library(singscore); library(survival); library(jsonlite)})
S <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); a <- commandArgs(TRUE)
C <- readRDS(file.path(S, "cohorts_corrected.rds")); CO <- c("SU2C-mCRPC-union", "WCDT-JCO"); LAB <- c("SU2C mCRPC (OS)", "WCDT JCO-only (OS)")
P <- lapply(setNames(CO, CO), function(n) list(rk=rankGenes(C[[n]]$mat), map=C[[n]]$map, cli=C[[n]]$cli))
mods <- a[-1]; png(a[1], width=1100, height=430 * length(mods), res=110); par(mfrow=c(length(mods), 2), mar=c(4.2, 4.2, 4.2, 1))
for (m in mods) {
  lab <- sub("=.*", "", m); ac <- strsplit(sub(".*=", "", m), ":")[[1]]
  g <- unique(fromJSON(file.path(S, ac[1], "candidates_human.json"))[[ac[2]]]); meta <- read.csv(file.path(S, ac[1], "candidates_meta.csv"))
  nm <- meta$n_mouse[meta$cand == ac[2]]
  for (i in seq_along(CO)) { p <- P[[CO[i]]]; gg <- intersect(g, rownames(p$rk))
    sc <- setNames(simpleScore(p$rk, upSet=gg, centerScore=FALSE)$TotalScore, colnames(p$rk)); id <- p$map[names(sc)]; ok <- !is.na(id)
    x <- tapply(sc[ok], id[ok], mean); d <- merge(data.frame(id=names(x), score=as.numeric(x)), p$cli, by="id"); d <- d[complete.cases(d) & d$time > 0, ]
    d$grp <- factor(ifelse(d$score > median(d$score), "high", "low"), levels=c("low", "high"))
    hr <- exp(coef(coxph(Surv(time, event) ~ grp, data=d)))[1]; pv <- 1 - pchisq(survdiff(Surv(time, event) ~ grp, data=d)$chisq, 1)
    f <- survfit(Surv(time / 30.44, event) ~ grp, data=d)   # cohort times are days; drawn in months (HR / p unaffected)
    plot(f, col=c("#2c7bb6", "#d7191c"), lwd=2, mark.time=TRUE, xlab="months", ylab="overall survival", las=1)
    title(main=sprintf("%s  |  %s\n%d mouse genes; median split: HR %.2f, log-rank p %.3f", lab, LAB[i], nm, hr, pv), cex.main=0.85, font.main=1)
    legend("topright", sprintf("%s (n=%d)", c("low", "high"), as.integer(table(d$grp))), col=c("#2c7bb6", "#d7191c"), lwd=2, bty="n", cex=0.8) } }
invisible(dev.off()); cat("wrote", a[1], "\n")

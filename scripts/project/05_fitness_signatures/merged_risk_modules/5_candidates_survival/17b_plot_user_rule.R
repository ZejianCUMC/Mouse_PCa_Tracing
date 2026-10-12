#!/usr/bin/env Rscript
# opt5 step 17b — figures for the FSU modules passing Zejian's survival rule (usage: 17b_plot_user_rule.R <arm>).
# Per module, one row of 4 panels: KM curves in SU2C-mCRPC-union, WCDT-merged (ComBat-seq on raw counts -> TPM) and WCDT-JCO, using
# the split under which the module passes (median / mean: high vs low; tertile: T3 vs T1; q1q4: Q4 vs Q1; cont: median split drawn, with
# the continuous Cox HR per SD annotated); log-rank p and Cox HR on each panel. 4th panel: mouse per-tumour module score (z) by
# treatment and genotype (17a), with the OLS treatment coefficients of step 16. Scores = singscore, as in the scan.
# Output <arm>/figs/user_rule_<cand>.png (one per module) and <arm>/figs/user_rule_all.pdf. Env: rMAST.
suppressPackageStartupMessages({library(singscore); library(survival); library(data.table); library(maxstat)})
S <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); A <- file.path(S, commandArgs(TRUE)[1]); dir.create(file.path(A, "figs"), showWarnings=FALSE)
C <- readRDS(file.path(S, "cohorts_corrected.rds")); R <- fread(file.path(A, "user_rule_passers.csv")); TS <- fread(file.path(A, "user_rule_tumour_scores.csv"))
H <- jsonlite::fromJSON(file.path(A, "candidates_human.json"))
CO <- c("SU2C-mCRPC-union"="SU2C mCRPC (OS)", "WCDT-merged"="WCDT merged, ComBat-seq (OS)", "WCDT-JCO"="WCDT JCO-only (OS)")
P <- lapply(names(CO), function(n) list(rk=rankGenes(C[[n]]$mat), map=C[[n]]$map, cli=C[[n]]$cli)); names(P) <- names(CO)
frame <- function(p, g) { x <- setNames(simpleScore(p$rk, upSet=intersect(g, rownames(p$rk)), centerScore=FALSE)$TotalScore, colnames(p$rk)); id <- p$map[names(x)]; ok <- !is.na(id)
  v <- tapply(x[ok], id[ok], mean); d <- merge(data.frame(id=names(v), score=as.numeric(v)), p$cli, by="id"); d <- d[complete.cases(d) & d$time > 0, ]; d$time_m <- d$time / 30.44; d }
grp <- function(x, sp, d=NULL) switch(sp, median=, cont=ifelse(x > median(x), "high", "low"), mean=ifelse(x > mean(x), "high", "low"),
  optimal={ ms <- maxstat.test(Surv(time, event) ~ score, data=d, smethod="LogRank", pmethod="Lau92", minprop=.2, maxprop=.8); ifelse(x > ms$estimate, "high", "low") },
  tertile={ q <- quantile(x, c(1/3, 2/3)); ifelse(x <= q[1], "low", ifelse(x > q[2], "high", NA)) }, q1q4={ q <- quantile(x, c(.25, .75)); ifelse(x <= q[1], "low", ifelse(x >= q[2], "high", NA)) })
lab <- c(median="median split", mean="mean split", optimal="maximally selected cut (maxstat, Lau92 p in scan)", tertile="top vs bottom tertile", q1q4="top vs bottom quartile", cont="median split shown; HR per SD from continuous score")
panel <- function(d, sp, title) { d$g <- factor(grp(d$score, sp, d), levels=c("low", "high")); dd <- d[!is.na(d$g), ]
  f <- survfit(Surv(time_m, event) ~ g, data=dd); lp <- 1 - pchisq(survdiff(Surv(time_m, event) ~ g, data=dd)$chisq, 1); hr <- exp(coef(coxph(Surv(time_m, event) ~ g, data=dd)))
  plot(f, col=c("#3b75af", "#c0392b"), lwd=2, xlab="months", ylab="overall survival", mark.time=TRUE); title(main=title, line=2.5, cex.main=.95)
  txt <- c(sprintf("high vs low: HR %.2f, log-rank p %.3f", hr, lp), "")
  if (sp == "cont") { s <- summary(coxph(Surv(time_m, event) ~ scale(score), data=d))$coefficients; txt <- c(sprintf("continuous score: HR %.2f per SD, p %.3f", exp(s[1, 1]), s[1, 5]), sprintf("(median split drawn: HR %.2f, p %.3f)", hr, lp)) }
  legend("topright", c(sprintf("low (n=%d)", sum(dd$g == "low")), sprintf("high (n=%d)", sum(dd$g == "high"))), col=c("#3b75af", "#c0392b"), lwd=2, bty="n", cex=.85)
  mtext(txt[1], side=3, line=1.2, cex=.72); if (nzchar(txt[2])) mtext(txt[2], side=3, line=.25, cex=.68) }
mouse <- function(c, r) { d <- TS[cand == c]; d[, x := as.integer(factor(treatment, levels=c("Intact", "Castration", "Castration+ENZA"))) + ifelse(model == "Npp53", -.15, .15)]
  set.seed(42); d[, x := x + runif(.N, -.05, .05)]
  plot(d$x, d$z, pch=ifelse(d$model == "Npp53", 16, 17), col=ifelse(d$model == "Npp53", "#2c7fb8", "#d95f0e"), xaxt="n", xlim=c(.5, 3.5), xlab="", ylab="module score (z, per tumour)",
       main=""); title(main="Mouse tumours: score by treatment", line=2.5, cex.main=.95)
  axis(1, at=1:3, labels=c("Intact", "Castration", "Cas+ENZA"), cex.axis=.85)
  for (t in 1:3) for (m in c("Npp53", "Npp53Rb1")) { v <- d[x %between% c(t - .4, t + .4) & model == m, z]; if (length(v)) segments(t + ifelse(m == "Npp53", -.15, .15) - .12, median(v), t + ifelse(m == "Npp53", -.15, .15) + .12, median(v), lwd=2.5) }
  legend("topleft", c("Npp53", "Npp53Rb1"), pch=c(16, 17), col=c("#2c7fb8", "#d95f0e"), bty="n", cex=.8)
  mtext(sprintf("vs Intact, all tumours: Cas %s; Cas+ENZA %s", r$Cas_all, r$CasENZA_all), side=3, line=1.2, cex=.68)
  mtext(sprintf("within %s: Cas %s; Cas+ENZA %s", r$genotype, r$Cas_own, r$CasENZA_own), side=3, line=.25, cex=.68) }
pdf(file.path(A, "figs", "user_rule_all.pdf"), width=19, height=4.8)
for (i in seq_len(nrow(R))) { r <- R[i]; c <- r$cand; g <- H[[c]]
  for (dev in c("png", "pdf")) {
    if (dev == "png") png(file.path(A, "figs", sprintf("user_rule_%s.png", c)), width=2300, height=600, res=120) else dev.set(2)
    par(mfrow=c(1, 4), mar=c(4.2, 4.2, 4.3, 1), oma=c(0, 0, 3.2, 0))
    for (cn in names(CO)) panel(frame(P[[cn]], g), r$best_split, CO[[cn]])
    mouse(c, r)
    gl <- strsplit(r$genes, " ")[[1]]; gs <- if (length(gl) > 22) paste(c(gl[1:22], sprintf("... (+%d)", length(gl) - 22)), collapse=" ") else paste(gl, collapse=" ")
    mtext(sprintf("%s  |  %s FSU module, %d genes  |  passes under: %s (%s)", c, r$genotype, r$n_genes, r$splits_passing, lab[[r$best_split]]), outer=TRUE, line=1.7, cex=.95, font=2)
    mtext(gs, outer=TRUE, line=.4, cex=.72)
    if (dev == "png") invisible(dev.off()) } }
invisible(dev.off()); cat(nrow(R), "module figures written to", file.path(A, "figs"), "\n")

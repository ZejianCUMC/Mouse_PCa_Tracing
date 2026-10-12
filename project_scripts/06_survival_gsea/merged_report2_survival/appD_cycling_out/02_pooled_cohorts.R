#!/usr/bin/env Rscript
# Appendix D (cycling cells excluded) (pooled): median-split HR (p) of the top pooled candidates in all five cohorts, scored as in ../01_km_all_cohorts.R. Output pooled_top_cohorts.csv. Env: rMAST.
suppressPackageStartupMessages({library(singscore); library(survival); library(jsonlite); library(data.table)})
O <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); S <- dirname(dirname(dirname(O)))
C <- readRDS(file.path(S, "cohorts_corrected.rds")); G <- fromJSON(file.path(O, "pooled_top_human.json"))
CO <- c("SU2C-mCRPC-union", "WCDT-JCO", "TCGA-PRAD", "DKFZ-primary", "MSKCC-primary")
P <- lapply(setNames(CO, CO), function(n) list(rk=rankGenes(C[[n]]$mat), map=C[[n]]$map, cli=C[[n]]$cli)); res <- list()
for (k in names(G)) for (n in CO) { p <- P[[n]]; s <- setNames(simpleScore(p$rk, upSet=intersect(G[[k]], rownames(p$rk)), centerScore=FALSE)$TotalScore, colnames(p$rk))
  id <- p$map[names(s)]; ok <- !is.na(id); x <- tapply(s[ok], id[ok], mean); d <- merge(data.frame(id=names(x), score=as.numeric(x)), p$cli, by="id"); d <- d[complete.cases(d) & d$time > 0, ]
  d$grp <- factor(ifelse(d$score > median(d$score), "high", "low"), levels=c("low", "high"))
  res[[length(res) + 1]] <- data.table(set=k, cohort=n, HR=exp(coef(coxph(Surv(time, event) ~ grp, data=d)))[1], p=1 - pchisq(survdiff(Surv(time, event) ~ grp, data=d)$chisq, 1), n=nrow(d)) }
R <- rbindlist(res); fwrite(R, file.path(O, "pooled_top_cohorts.csv")); print(dcast(R[, .(set, cohort, v=sprintf("%.2f (%.3f)", HR, p))], set ~ cohort, value.var="v"))

#!/usr/bin/env Rscript
# opt5 survival scanner (usage: km_scan.R <candidates json> <out csv> [maxstat 0/1] [sel|val]).
# Cohorts from cohorts_corrected.rds (step 00, audited inputs):
#   sel = SU2C-mCRPC-union (OS), WCDT-merged (OS; ComBat-seq on raw counts -> TPM; primary WCDT), WCDT-JCO (OS; raw counts -> TPM)
#   val = TCGA-PRAD (PFI), DKFZ-primary (BCR), MSKCC-primary (BCR)   [Michigan dropped: z-scores only]
# Sample score = singscore simpleScore (uncentred; ranks within sample over the cohort's protein-coding genes); patient score = mean
# over that patient's samples. Candidate = gene array (up-set score) or {up, dn} (bidirectional). Splits: median, mean, tertile
# (T3 vs T1), q1q4 (log-rank p; HR from Cox on the groups); cont = Cox per SD (Wald p); optional optimal = maxstat Lau92 (minprop 0.2).
# HR is always high vs low score (direction-aware reading in 08). Deterministic. Env: rMAST.
suppressPackageStartupMessages({library(singscore); library(survival); library(maxstat); library(data.table); library(jsonlite); library(parallel)})
S <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); a <- commandArgs(TRUE)
MS <- length(a) > 2 && a[3] == "1"; SET <- if (length(a) > 3) a[4] else "sel"
C <- readRDS(file.path(S, "cohorts_corrected.rds"))
CO <- if (SET == "sel") c("SU2C-mCRPC-union", "WCDT-merged", "WCDT-JCO") else c("TCGA-PRAD", "DKFZ-primary", "MSKCC-primary")
P <- lapply(setNames(CO, CO), function(n) { co <- C[[n]]; list(rk=rankGenes(co$mat), map=co$map, cli=co$cli) })
agg <- function(p, sc) { id <- p$map[names(sc)]; ok <- !is.na(id); x <- tapply(sc[ok], id[ok], mean)
  d <- merge(data.frame(id=names(x), score=as.numeric(x)), p$cli, by="id"); d[complete.cases(d) & d$time > 0, ] }
cl <- function(g) unique(g[nchar(g) > 0]); cand <- lapply(fromJSON(a[1]), function(g) if (is.list(g)) lapply(g, cl) else cl(g))
lr <- function(d, g) { d$grp <- factor(g, levels=c("low", "high")); d <- d[!is.na(d$grp), ]; if (length(unique(d$grp)) < 2) return(c(NA, NA, NA))
  c(exp(coef(coxph(Surv(time, event) ~ grp, data=d)))[1], 1 - pchisq(survdiff(Surv(time, event) ~ grp, data=d)$chisq, 1), nrow(d)) }
one <- function(k) { out <- list()
  for (cn in names(P)) { rn <- rownames(P[[cn]]$rk); g <- cand[[k]]
    if (is.list(g)) { up <- intersect(g$up, rn); dn <- intersect(g$dn, rn); if (length(up) < 3 || length(dn) < 3) next; ng <- length(up) + length(dn)
      ss <- simpleScore(P[[cn]]$rk, upSet=up, downSet=dn, centerScore=FALSE)$TotalScore
    } else { p <- intersect(g, rn); if (length(p) < 3) next; ng <- length(p); ss <- simpleScore(P[[cn]]$rk, upSet=p, centerScore=FALSE)$TotalScore }
    d <- agg(P[[cn]], setNames(ss, colnames(P[[cn]]$rk))); x <- d$score; q3 <- quantile(x, c(1/3, 2/3)); q4 <- quantile(x, c(.25, .75))
    G <- list(median=ifelse(x > median(x), "high", "low"), mean=ifelse(x > mean(x), "high", "low"),
              tertile=ifelse(x <= q3[1], "low", ifelse(x > q3[2], "high", NA)), q1q4=ifelse(x <= q4[1], "low", ifelse(x >= q4[2], "high", NA)))
    for (sp in names(G)) { r <- lr(d, G[[sp]]); out[[length(out) + 1]] <- data.frame(cand=k, cohort=cn, split=sp, HR=r[1], p=r[2], n=r[3], n_genes=ng) }
    d$zs <- as.numeric(scale(x)); cf <- summary(coxph(Surv(time, event) ~ zs, data=d))$coefficients
    out[[length(out) + 1]] <- data.frame(cand=k, cohort=cn, split="cont", HR=exp(cf[1, 1]), p=cf[1, 5], n=nrow(d), n_genes=ng)
    if (MS) { ms <- try(maxstat.test(Surv(time, event) ~ score, data=d, smethod="LogRank", pmethod="Lau92", minprop=.2, maxprop=.8), silent=TRUE)
      if (!inherits(ms, "try-error")) { r <- lr(d, ifelse(x > ms$estimate, "high", "low")); out[[length(out) + 1]] <- data.frame(cand=k, cohort=cn, split="optimal", HR=r[1], p=ms$p.value, n=r[3], n_genes=ng) } } }
  if (length(out)) do.call(rbind, out) else NULL }
t0 <- Sys.time(); res <- rbindlist(mclapply(names(cand), one, mc.cores=as.integer(Sys.getenv("NCORES", "12"))))
fwrite(res, a[2]); cat(length(cand), "candidates,", nrow(res), "rows (", SET, "),", format(Sys.time() - t0), "\n")

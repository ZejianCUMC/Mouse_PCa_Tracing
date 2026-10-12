#!/usr/bin/env Rscript
# opt5 step 14 — how much of each candidate is proliferation? (usage: 14_prolif_corr.R <arm>). For every non-bidirectional candidate:
# Spearman correlation of its patient-level singscore with the Hallmark E2F score in SU2C-mCRPC-union and WCDT-merged, and the Cox HR per
# SD adjusted for E2F (risk scale: FSD sign-flipped). Lets 15_summary stratify passers into proliferation-type (rho >= 0.7) and
# non-proliferation (rho < 0.3) modules. Output <arm>/prolif_corr.csv. Env: rMAST.
suppressPackageStartupMessages({library(singscore); library(survival); library(data.table); library(jsonlite); library(parallel)})
S <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); A <- file.path(S, commandArgs(TRUE)[1])
C <- readRDS(file.path(S, "cohorts_corrected.rds")); HS <- fromJSON(file.path(dirname(S), "opt3_cycout_20261009/hallmark_sets.json"))
H <- fromJSON(file.path(A, "candidates_human.json")); M <- fread(file.path(A, "candidates_meta.csv")); H <- H[M[direction != "BI", cand]]; dir_of <- setNames(M$direction, M$cand)
P <- lapply(c(SU2C="SU2C-mCRPC-union", WCDTm="WCDT-merged"), function(n) { p <- list(rk=rankGenes(C[[n]]$mat), map=C[[n]]$map, cli=C[[n]]$cli)
  e <- setNames(simpleScore(p$rk, upSet=intersect(HS[["HALLMARK_E2F_TARGETS"]], rownames(p$rk)), centerScore=FALSE)$TotalScore, colnames(p$rk))
  id <- p$map[names(e)]; ok <- !is.na(id); ev <- tapply(e[ok], id[ok], mean); p$d <- merge(data.frame(id=names(ev), E2F=as.numeric(scale(as.numeric(ev)))), p$cli, by="id"); p$d <- p$d[p$d$time > 0, ]; p })
one <- function(k) { r <- list(cand=k)
  for (cn in names(P)) { p <- P[[cn]]; g <- intersect(H[[k]], rownames(p$rk)); if (length(g) < 3) next
    x <- setNames(simpleScore(p$rk, upSet=g, centerScore=FALSE)$TotalScore, colnames(p$rk)); id <- p$map[names(x)]; ok <- !is.na(id); v <- tapply(x[ok], id[ok], mean)
    d <- merge(p$d, data.frame(id=names(v), s=as.numeric(v)), by="id"); d$risk <- as.numeric(scale(if (dir_of[[k]] == "FSD") -d$s else d$s))
    cf <- summary(coxph(Surv(time, event) ~ risk + E2F, data=d))$coefficients
    r[[paste0("rho_", cn)]] <- cor(d$s, d$E2F, method="spearman"); r[[paste0("HRadj_", cn)]] <- exp(cf["risk", "coef"]); r[[paste0("padj_", cn)]] <- cf["risk", "Pr(>|z|)"] }
  as.data.frame(r) }
R <- rbindlist(mclapply(names(H), one, mc.cores=as.integer(Sys.getenv("NCORES", "12"))), fill=TRUE); fwrite(R, file.path(A, "prolif_corr.csv")); cat(basename(A), nrow(R), "candidates\n")

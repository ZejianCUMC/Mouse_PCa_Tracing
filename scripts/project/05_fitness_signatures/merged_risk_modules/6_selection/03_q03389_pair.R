#!/usr/bin/env Rscript
# Matched pair from the identical setting (bcr arm, list FSU_top1000, Hotspot k30 / top 1000, create_modules min 30, adjacent pair):
# Npp53 Q03389 and Npp53Rb1 Q03123. Survival (median, mean), E2F rho / adjusted Cox, score concordance and gene overlap (incl. Q03340).
suppressPackageStartupMessages({library(singscore); library(survival); library(jsonlite); library(data.table)})
O <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); S <- dirname(O)
C <- readRDS(file.path(S, "cohorts_corrected.rds")); HS <- fromJSON(file.path(dirname(S), "opt3_cycout_20261009/hallmark_sets.json"))
np <- fromJSON(file.path(S, "metaNorm_20261010_Npp53_bcr/candidates_human.json")); rb <- fromJSON(file.path(S, "metaNorm_20261010_Npp53Rb1_bcr/candidates_human.json"))
G <- list(Q03389=unique(np[["Q03389"]]), Q03123=unique(rb[["Q03123"]]), Q03340=unique(np[["Q03340"]]))
jac <- function(a, b) sprintf("%d shared, Jaccard %.3f", length(intersect(a, b)), length(intersect(a, b)) / length(union(a, b)))
cat("Q03389 vs Q03123:", jac(G$Q03389, G$Q03123), "\nQ03389 vs Q03340:", jac(G$Q03389, G$Q03340), "\n")
cat("Q03389 E2F/G2M members:", length(intersect(G$Q03389, HS[["HALLMARK_E2F_TARGETS"]])), length(intersect(G$Q03389, HS[["HALLMARK_G2M_CHECKPOINT"]])), "\n")
out <- list()
for (n in c("SU2C-mCRPC-union", "WCDT-JCO", "WCDT-merged")) { co <- C[[n]]; rk <- rankGenes(co$mat); cli <- co$cli[complete.cases(co$cli) & co$cli$time > 0, ]
  sc <- function(g) { s <- setNames(simpleScore(rk, upSet=intersect(g, rownames(rk)), centerScore=FALSE)$TotalScore, colnames(rk)); id <- co$map[names(s)]; ok <- !is.na(id); tapply(s[ok], id[ok], mean) }
  a <- sc(G$Q03389); b <- sc(G$Q03123); e <- sc(HS[["HALLMARK_E2F_TARGETS"]]); ids <- intersect(intersect(names(a), names(e)), cli$id)
  d <- merge(data.frame(id=ids, s=as.numeric(a[ids]), E2F=as.numeric(e[ids])), cli, by="id"); d$risk <- as.numeric(scale(d$s)); d$E2Fz <- as.numeric(scale(d$E2F))
  cf <- summary(coxph(Surv(time, event) ~ risk + E2Fz, data=d))$coefficients
  out[[n]] <- data.frame(cohort=n, Q03389_rho_E2F=round(cor(d$s, d$E2F, method="spearman"), 3), Q03389_HRadj=round(exp(cf["risk", "coef"]), 3), Q03389_padj=round(cf["risk", "Pr(>|z|)"], 3),
                         spearman_Q03389_vs_Q03123=round(cor(a[ids], b[ids], method="spearman"), 3), same_median_group=sum((a[ids] > median(a[ids])) == (b[ids] > median(b[ids]))), n=length(ids)) }
R <- do.call(rbind, out); fwrite(R, file.path(O, "pair_q03389_q03123.csv")); print(R)
K <- fread(file.path(S, "metaNorm_20261010_Npp53_bcr/km_scan.csv"))[cand == "Q03389" & split %in% c("median", "mean")]; fwrite(K, file.path(O, "q03389_km.csv")); print(K[, .(cohort, split, HR=round(HR, 2), p=round(p, 3))])

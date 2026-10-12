#!/usr/bin/env Rscript
# Q03123 (Npp53Rb1, best candidate of Q03340's Hotspot run): all splits, E2F rho / E2F-adjusted Cox in SU2C, WCDT-JCO, WCDT-merged,
# score correlation with Q03340, human-ortholog Jaccard with Q03340 / Q02570 / Q02660. Env: rMAST.
suppressPackageStartupMessages({library(singscore); library(survival); library(jsonlite); library(data.table)})
O <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); S <- dirname(O)
C <- readRDS(file.path(S, "cohorts_corrected.rds")); HS <- fromJSON(file.path(dirname(S), "opt3_cycout_20261009/hallmark_sets.json"))
rb <- fromJSON(file.path(S, "metaNorm_20261010_Npp53Rb1_bcr/candidates_human.json")); np <- fromJSON(file.path(S, "metaNorm_20261010_Npp53_bcr/candidates_human.json"))
G <- list(Q03123=unique(rb[["Q03123"]]), Q03340=unique(np[["Q03340"]]), Q02570=unique(rb[["Q02570"]]), Q02660=unique(rb[["Q02660"]]))
J <- sapply(c("Q03340", "Q02570", "Q02660"), function(k) round(length(intersect(G$Q03123, G[[k]])) / length(union(G$Q03123, G[[k]])), 3))
sh <- sapply(c("Q03340", "Q02570", "Q02660"), function(k) length(intersect(G$Q03123, G[[k]])))
res <- list(); cc <- list()
for (n in c("SU2C-mCRPC-union", "WCDT-JCO", "WCDT-merged")) { co <- C[[n]]; rk <- rankGenes(co$mat); cli <- co$cli[complete.cases(co$cli) & co$cli$time > 0, ]
  sc <- function(g) { s <- setNames(simpleScore(rk, upSet=intersect(g, rownames(rk)), centerScore=FALSE)$TotalScore, colnames(rk)); id <- co$map[names(s)]; ok <- !is.na(id); tapply(s[ok], id[ok], mean) }
  a <- sc(G$Q03123); b <- sc(G$Q03340); e <- sc(HS[["HALLMARK_E2F_TARGETS"]]); ids <- intersect(intersect(names(a), names(e)), cli$id)
  d <- merge(data.frame(id=ids, s=as.numeric(a[ids]), E2F=as.numeric(e[ids])), cli, by="id"); d$risk <- as.numeric(scale(d$s)); d$E2Fz <- as.numeric(scale(d$E2F))
  cf <- summary(coxph(Surv(time, event) ~ risk + E2Fz, data=d))$coefficients
  res[[n]] <- data.frame(cohort=n, rho_E2F=round(cor(d$s, d$E2F, method="spearman"), 3), HRadj=round(exp(cf["risk", "coef"]), 3), padj=round(cf["risk", "Pr(>|z|)"], 3))
  cc[[n]] <- data.frame(cohort=n, spearman_vs_Q03340=round(cor(a[ids], b[ids], method="spearman"), 3), same_median_group=sum((a[ids] > median(a[ids])) == (b[ids] > median(b[ids]))), n=length(ids)) }
E <- do.call(rbind, res); Cc <- do.call(rbind, cc); K <- fread(file.path(S, "metaNorm_20261010_Npp53Rb1_bcr/km_scan.csv"))[cand == "Q03123"]
fwrite(E, file.path(O, "q03123_e2f.csv")); fwrite(Cc, file.path(O, "q03123_concordance_vs_Q03340.csv")); fwrite(K, file.path(O, "q03123_km_all_splits.csv"))
print(E); print(Cc); print(K[split != "optimal", .(cohort, split, HR=round(HR, 2), p=signif(p, 3), n)]); cat("Jaccard (human):", paste(names(J), J, "shared", sh), sep="\n")
cat("E2F/G2M members:", length(intersect(G$Q03123, HS[["HALLMARK_E2F_TARGETS"]])), length(intersect(G$Q03123, HS[["HALLMARK_G2M_CHECKPOINT"]])), "\n")

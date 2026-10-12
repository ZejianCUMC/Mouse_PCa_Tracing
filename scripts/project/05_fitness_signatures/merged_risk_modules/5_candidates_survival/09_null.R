#!/usr/bin/env Rscript
# opt5 step 09 — chance calibration (usage: 09_null.R <arm>). 1,000 random gene sets per size (15, 30, 60, 120) drawn (set.seed(42))
# from <arm>/pool_human.txt (human orthologs of every gene in the arm's primary FSU / FSD lists) restricted to genes measured in all
# three selection cohorts; scored with km_scan.R (sel, no maxstat) and, for the finalist calibration (step 11), the validation cohorts.
# Pass rates for both risk readings (high score adverse = FSU / BI rule; low score adverse = FSD rule).
# Output <arm>/null_sets.json, km_null.csv, km_null_val.csv, null_pass_rates.csv. Env: rMAST.
suppressPackageStartupMessages({library(jsonlite); library(data.table)})
S <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))); A <- file.path(S, commandArgs(TRUE)[1])
C <- readRDS(file.path(S, "cohorts_corrected.rds")); SEL <- c("SU2C-mCRPC-union", "WCDT-merged", "WCDT-JCO")
pool <- intersect(readLines(file.path(A, "pool_human.txt")), Reduce(intersect, lapply(SEL, function(n) rownames(C[[n]]$mat))))
set.seed(42); sets <- list(); for (sz in c(15, 30, 60, 120)) for (i in 1:1000) sets[[sprintf("rand_s%d_%04d", sz, i)]] <- sample(pool, sz)
cat(basename(A), "pool:", length(pool), "| sets:", length(sets), "\n"); jf <- file.path(A, "null_sets.json"); write_json(sets, jf)
R <- file.path(R.home("bin"), "Rscript")
system2(R, c(file.path(S, "km_scan.R"), jf, file.path(A, "km_null.csv"), "0", "sel")); system2(R, c(file.path(S, "km_scan.R"), jf, file.path(A, "km_null_val.csv"), "0", "val"))
K <- fread(file.path(A, "km_null.csv")); K[, cohort := gsub("-", "_", cohort)]; W <- dcast(K, cand + split ~ cohort, value.var=c("HR", "p"))
hi <- function(c, a) W[[paste0("HR_", c)]] > 1 & W[[paste0("p_", c)]] < a; lo <- function(c, a) W[[paste0("HR_", c)]] < 1 & W[[paste0("p_", c)]] < a
W[, `:=`(strict=hi("SU2C_mCRPC_union", .05) & hi("WCDT_merged", .05), relaxed=hi("SU2C_mCRPC_union", .05) & hi("WCDT_merged", .1),
         strict_low=lo("SU2C_mCRPC_union", .05) & lo("WCDT_merged", .05), relaxed_low=lo("SU2C_mCRPC_union", .05) & lo("WCDT_merged", .1))]
W[, size := as.integer(sub("_.*", "", sub("^rand_s", "", cand)))]
RT <- W[, .(strict=mean(strict, na.rm=TRUE), relaxed=mean(relaxed, na.rm=TRUE), strict_low=mean(strict_low, na.rm=TRUE), relaxed_low=mean(relaxed_low, na.rm=TRUE)), by=.(size, split)]
fwrite(RT, file.path(A, "null_pass_rates.csv"))
AN <- W[, .(s=any(strict, na.rm=TRUE), r=any(relaxed, na.rm=TRUE), sl=any(strict_low, na.rm=TRUE), rl=any(relaxed_low, na.rm=TRUE),
            s_mm=any(strict & split %in% c("median", "mean"), na.rm=TRUE), r_mm=any(relaxed & split %in% c("median", "mean"), na.rm=TRUE),
            sl_mm=any(strict_low & split %in% c("median", "mean"), na.rm=TRUE), rl_mm=any(relaxed_low & split %in% c("median", "mean"), na.rm=TRUE)), by=cand]
cat(sprintf("%s random-set pass rate, any split — high adverse: strict %.4f relaxed %.4f | low adverse: strict %.4f relaxed %.4f\n", basename(A), mean(AN$s), mean(AN$r), mean(AN$sl), mean(AN$rl)))
cat(sprintf("%s random-set pass rate, median/mean — high adverse: strict %.4f relaxed %.4f | low adverse: strict %.4f relaxed %.4f\n", basename(A), mean(AN$s_mm), mean(AN$r_mm), mean(AN$sl_mm), mean(AN$rl_mm)))
fwrite(AN, file.path(A, "null_any_split.csv"))

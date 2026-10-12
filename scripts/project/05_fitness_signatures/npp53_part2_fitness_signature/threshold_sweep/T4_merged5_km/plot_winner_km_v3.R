#!/usr/bin/env Rscript
# Fig 9.5 v3 — split into TWO figures:
#   (A) Fig 9.5: 2x2 MEAN-split KM curves (4 curated cohorts)
#   (B) Fig 9.6: remaining splits (median + q1q4) in a 2 rows x 4 cols grid
#       (row1 = median, row2 = q1q4; cols = TCGA-PRAD PFI | SU2C-mCRPC OS | WCDT-mCRPC OS | DKFZ-primary BCR)
# All settings inherited from v2: TCGA and WCDT capped at 60 mo (ticks 0/15/30/45/60), larger title font.
suppressPackageStartupMessages({
  library(singscore); library(survival); library(survminer)
  library(data.table); library(ggplot2); library(dplyr); library(cowplot)
})
plog <- function(...) cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n", sep="")

NPP53 <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
CACHE <- file.path(NPP53, "analysis/fitness_signatures_survival/cohort_cache_pc.rds")
DIR   <- file.path(NPP53, "analysis/fitness_signatures_survival/threshold_sweep/T4/merged5")
MOD_ID <- "M3_merged"
MOD_LABEL_DISPLAY <- "M3"
MOD_H <- file.path(DIR, "module_3_human.txt")

KEEP <- list(
  list(cohort="TCGA-PRAD",         endpoint="PFI", display="TCGA-PRAD"),
  list(cohort="SU2C-mCRPC-union",  endpoint="OS",  display="SU2C-mCRPC"),
  list(cohort="WCDT-mCRPC-union",  endpoint="OS",  display="WCDT-mCRPC"),
  list(cohort="DKFZ-primary",      endpoint="BCR", display="DKFZ-primary")
)

plog("Loading cache + winner human gene list ...")
COH <- readRDS(CACHE)
MG  <- unique(readLines(MOD_H)); MG <- MG[nchar(MG)>0]
plog("gene set = ", length(MG), " human genes")

score_singscore <- function(expr, genes) {
  p <- intersect(genes, rownames(expr)); if (length(p)<3) return(NULL)
  m <- as.matrix(expr); if (any(duplicated(colnames(m)))) m <- m[,!duplicated(colnames(m))]
  rk <- singscore::rankGenes(m); ss <- singscore::simpleScore(rk, upSet=p, centerScore=FALSE)
  setNames(ss$TotalScore, rownames(ss))
}
build_df <- function(cohort_name, ep) {
  co <- COH[[cohort_name]]
  score <- score_singscore(co$mat, MG); if (is.null(score)) return(NULL)
  if (cohort_name=="TCGA-PRAD") {
    df <- data.frame(sample_id=as.character(names(score)), score=as.numeric(score), stringsAsFactors=FALSE)
    co$cli$sample_id <- as.character(co$cli$sample_id); d <- merge(df, co$cli, by="sample_id")
    tc <- switch(ep, OS=c("OS.time","OS"), PFI=c("PFI.time","PFI"), DFI=c("DFI.time","DFI"), DSS=c("DSS.time","DSS"))
  } else if (grepl("^WCDT", cohort_name)) {
    df <- data.frame(sample_id=as.character(names(score)), score=as.numeric(score), stringsAsFactors=FALSE)
    co$cli$sample_id <- as.character(co$cli$sample_id); d <- merge(df, co$cli, by="sample_id"); tc <- c("OS.time","OS")
  } else {
    patient <- co$map[names(score)]; ok <- !is.na(patient)
    tmp <- data.frame(sample_id=as.character(patient[ok]), score=as.numeric(score)[ok], stringsAsFactors=FALSE) %>%
             group_by(sample_id) %>% summarise(score=mean(score), .groups="drop") %>% as.data.frame()
    co$cli$sample_id <- as.character(co$cli$sample_id); d <- merge(tmp, co$cli, by="sample_id")
    tc <- if (ep=="BCR") c("BCR.time","BCR") else c("OS.time","OS")
  }
  d$time <- d[[tc[1]]]; d$event <- d[[tc[2]]]
  d[complete.cases(d$time,d$event,d$score) & d$time>0, ]
}
apply_split <- function(d, sp) {
  if (sp=="q1q4") {
    q1 <- quantile(d$score,0.25); q4 <- quantile(d$score,0.75)
    d$grp <- factor(ifelse(d$score<=q1,"low", ifelse(d$score>=q4,"high", NA)), levels=c("low","high"))
    d <- d[!is.na(d$grp),]
  } else {
    thr <- if (sp=="mean") mean(d$score,na.rm=TRUE) else median(d$score,na.rm=TRUE)
    d$grp <- factor(ifelse(d$score>thr,"high","low"), levels=c("low","high"))
  }
  d
}
cell_stat <- function(d) {
  if (nrow(d)<10 || length(unique(d$grp))<2) return(list(HR=NA,p=NA))
  sd_ <- try(survdiff(Surv(time,event)~grp, data=d), silent=TRUE)
  cx  <- try(coxph(Surv(time,event)~grp, data=d), silent=TRUE)
  p_lr <- if (inherits(sd_,"try-error")) NA else 1 - pchisq(sd_$chisq, df=length(sd_$n)-1)
  HR <- if (inherits(cx,"try-error")) NA else exp(coef(cx))[1]
  list(HR=HR, p=p_lr)
}
make_plot <- function(cohort_code, cohort_display, ep, split) {
  d <- build_df(cohort_code, ep); if (is.null(d)) return(NULL)
  d <- apply_split(d, split); st <- cell_stat(d)
  if (is.na(st$HR)) return(NULL)
  d$time_mo <- d$time / 30.44
  fit <- survfit(Surv(time_mo, event) ~ grp, data=d)
  sig_mark <- if (!is.na(st$p) && st$p<0.05) " *" else ""
  ttl <- sprintf("%s · %s · %s%s\n(HR=%.2f · p=%.2g · n=%d)",
                 cohort_display, ep, split, sig_mark, st$HR, st$p, nrow(d))
  cap60 <- (cohort_code == "TCGA-PRAD") || grepl("^WCDT", cohort_code) || (cohort_code == "DKFZ-primary")
  xlim_val <- if (cap60) c(0, 60) else NULL
  break_by <- if (cap60) 15 else NULL
  ggsurvplot(fit, data=d,
             pval=TRUE, pval.method=FALSE, pval.size=5, pval.coord=c(3, 0.10),
             conf.int=TRUE, risk.table=TRUE, risk.table.height=0.30,
             risk.table.title="At risk", risk.table.y.text=FALSE,
             palette=c("#4a7fbf","#c44e52"),
             legend.title=MOD_LABEL_DISPLAY, legend.labs=c("low","high"),
             xlab="Time (months)", ylab="Survival probability",
             surv.median.line="none", xlim=xlim_val, break.time.by=break_by,
             ggtheme=theme_bw(base_size=10) +
               theme(plot.title=element_text(hjust=0, size=10, face="bold", color="#0f2747"),
                     legend.position="top", panel.grid.minor=element_blank()),
             title=ttl,
             font.legend=c(9,"plain"), font.tickslab=c(9), font.x=c(10), font.y=c(10))
}
compose <- function(p) {
  km_plt <- p$plot; km_tbl <- p$table + theme(plot.title=element_text(size=9))
  plot_grid(km_plt, km_tbl, ncol=1, rel_heights=c(0.72, 0.28))
}
render <- function(row, sp) {
  p <- try(make_plot(row$cohort, row$display, row$endpoint, sp), silent=TRUE)
  if (inherits(p,"try-error") || is.null(p)) ggplot() + theme_void() else compose(p)
}

# ============================================================================
# (A) Fig 9.5 — MEAN-split only, 2x2 grid
# ============================================================================
plots_mean <- lapply(KEEP, function(r) render(r, "mean"))
title_mean <- sprintf(
  "Fig 9.5 · Module M3 mean-split KM (proliferation-housekeeping · %d human genes) · singscore · 4 cohorts",
  length(MG))
grid_mean <- plot_grid(plotlist=plots_mean, ncol=2, nrow=2)
titled_mean <- plot_grid(
  ggdraw() + draw_label(title_mean, x=0.02, hjust=0, size=20, fontface='bold', colour="#0f2747"),
  grid_mean, ncol=1, rel_heights=c(0.04, 0.96))
fp_mean_pdf <- file.path(DIR, "KM_MERGED_WINNER_mean_2x2.pdf")
fp_mean_png <- file.path(DIR, "KM_MERGED_WINNER_mean_2x2.png")
ggsave(fp_mean_pdf, titled_mean, width=5.5*2, height=5.5*2 + 0.4, limitsize=FALSE)
ggsave(fp_mean_png, titled_mean, width=5.5*2, height=5.5*2 + 0.4, dpi=110, limitsize=FALSE)
plog("[saved] ", fp_mean_png, " (mean 2x2)")

# ============================================================================
# (B) Fig 9.6 — median + q1q4 splits, 2 rows x 4 cols (rows = split, cols = cohort)
# ============================================================================
plots_rest <- list()
for (sp in c("median","q1q4")) {
  for (r in KEEP) plots_rest[[length(plots_rest)+1]] <- render(r, sp)
}
title_rest <- sprintf(
  "Fig 9.6 · Module M3 median-split (top row) and Q1-vs-Q4 split (bottom row) KM (%d human genes) · cols = TCGA-PRAD PFI | SU2C-mCRPC OS | WCDT-mCRPC OS | DKFZ-primary BCR",
  length(MG))
grid_rest <- plot_grid(plotlist=plots_rest, ncol=4, nrow=2)
titled_rest <- plot_grid(
  ggdraw() + draw_label(title_rest, x=0.02, hjust=0, size=18, fontface='bold', colour="#0f2747"),
  grid_rest, ncol=1, rel_heights=c(0.03, 0.97))
fp_rest_pdf <- file.path(DIR, "KM_MERGED_WINNER_median_q1q4_2x4.pdf")
fp_rest_png <- file.path(DIR, "KM_MERGED_WINNER_median_q1q4_2x4.png")
ggsave(fp_rest_pdf, titled_rest, width=5.5*4, height=5.5*2 + 0.35, limitsize=FALSE)
ggsave(fp_rest_png, titled_rest, width=5.5*4, height=5.5*2 + 0.35, dpi=100, limitsize=FALSE)
plog("[saved] ", fp_rest_png, " (median+q1q4 2x4)")
plog("Done.")

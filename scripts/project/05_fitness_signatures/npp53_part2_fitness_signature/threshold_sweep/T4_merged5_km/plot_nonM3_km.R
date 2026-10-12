#!/usr/bin/env Rscript
# Non-M3 modules KM grid: 4 rows (M1/M2/M4/M5) × 4 cols (cohorts), mean-split only.
# Compact size for HTML report insertion after Fig 9.5.
suppressPackageStartupMessages({
  library(singscore); library(survival); library(survminer)
  library(data.table); library(ggplot2); library(dplyr); library(cowplot)
})
plog <- function(...) cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n", sep="")

NPP53 <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
CACHE <- file.path(NPP53, "analysis/fitness_signatures_survival/cohort_cache_pc.rds")
DIR   <- file.path(NPP53, "analysis/fitness_signatures_survival/threshold_sweep/T4/merged5")

MODS <- c(1, 2, 4, 5)
MOD_LABELS <- c("M1 Stromal/ECM", "M2 Glycolysis/Hypoxia",
                "M4 Squamous/Stress", "M5 Lipid/Chaperone")
names(MOD_LABELS) <- as.character(MODS)

KEEP <- list(
  list(cohort="TCGA-PRAD",        endpoint="PFI", display="TCGA-PRAD · PFI"),
  list(cohort="SU2C-mCRPC-union", endpoint="OS",  display="SU2C-mCRPC · OS"),
  list(cohort="WCDT-mCRPC-union", endpoint="OS",  display="WCDT-mCRPC · OS"),
  list(cohort="DKFZ-primary",     endpoint="BCR", display="DKFZ-primary · BCR")
)

plog("Loading cohort cache ...")
COH <- readRDS(CACHE)

# Load gene lists
GENE_LISTS <- list()
for (mid in MODS) {
  fp <- file.path(DIR, sprintf("module_%d_human.txt", mid))
  g <- unique(readLines(fp)); g <- g[nchar(g)>0]
  GENE_LISTS[[as.character(mid)]] <- g
  plog("  M", mid, ": ", length(g), " human genes")
}

score_singscore <- function(expr, genes) {
  p <- intersect(genes, rownames(expr)); if (length(p)<3) return(NULL)
  m <- as.matrix(expr); if (any(duplicated(colnames(m)))) m <- m[,!duplicated(colnames(m))]
  rk <- singscore::rankGenes(m); ss <- singscore::simpleScore(rk, upSet=p, centerScore=FALSE)
  setNames(ss$TotalScore, rownames(ss))
}
build_df <- function(cohort_name, ep, genes) {
  co <- COH[[cohort_name]]
  score <- score_singscore(co$mat, genes); if (is.null(score)) return(NULL)
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
make_plot <- function(mid, cohort_code, cohort_display, ep) {
  genes <- GENE_LISTS[[as.character(mid)]]
  d <- build_df(cohort_code, ep, genes); if (is.null(d)) return(ggplot() + theme_void())
  thr <- mean(d$score, na.rm=TRUE)
  d$grp <- factor(ifelse(d$score > thr, "high", "low"), levels=c("low","high"))
  if (nrow(d)<10 || length(unique(d$grp))<2) return(ggplot() + theme_void())
  sd_ <- try(survdiff(Surv(time,event)~grp, data=d), silent=TRUE)
  cx  <- try(coxph(Surv(time,event)~grp, data=d), silent=TRUE)
  p_lr <- if (inherits(sd_,"try-error")) NA else 1 - pchisq(sd_$chisq, df=length(sd_$n)-1)
  HR <- if (inherits(cx,"try-error")) NA else exp(coef(cx))[1]
  d$time_mo <- d$time / 30.44
  fit <- survfit(Surv(time_mo, event) ~ grp, data=d)
  sig <- if (!is.na(p_lr) && p_lr<0.05) " *" else ""
  mod_tag <- paste0("M", mid)
  ttl <- sprintf("%s · %s%s\nHR=%.2f  p=%.2g  n=%d", mod_tag, cohort_display, sig, HR, p_lr, nrow(d))
  cap60 <- (cohort_code == "TCGA-PRAD") || grepl("^WCDT", cohort_code) || (cohort_code == "DKFZ-primary")
  xlim_val <- if (cap60) c(0, 60) else NULL
  break_by <- if (cap60) 15 else NULL
  p <- ggsurvplot(fit, data=d,
             pval=TRUE, pval.size=4.5, pval.coord=c(1, 0.08), conf.int=TRUE,
             risk.table=FALSE,
             palette=c("#4a7fbf","#c44e52"),
             legend.title=mod_tag, legend.labs=c("low","high"),
             xlab="Time (months)", ylab="Surv. prob.",
             surv.median.line="none", xlim=xlim_val, break.time.by=break_by,
             ggtheme=theme_bw(base_size=8) +
               theme(plot.title=element_text(hjust=0, size=8, face="bold", color="#0f2747"),
                     legend.position="top",
                     legend.key.size=unit(0.3,"cm"),
                     legend.text=element_text(size=6),
                     legend.title=element_text(size=7, face="bold"),
                     legend.margin=margin(0,0,0,0),
                     legend.box.margin=margin(-5,-5,-5,-5),
                     panel.grid.minor=element_blank(),
                     axis.title=element_text(size=7),
                     axis.text=element_text(size=6),
                     plot.margin=margin(2,4,2,2)),
             title=ttl,
             font.legend=c(7,"plain"), font.tickslab=c(7), font.x=c(8), font.y=c(8))
  p$plot
}

# Build 4×4 grid: rows=modules, cols=cohorts
plog("Generating 4x4 KM grid ...")
plots <- list()
idx <- 0
for (mid in MODS) {
  for (k in KEEP) {
    idx <- idx + 1
    plog("  panel ", idx, "/16: M", mid, " × ", k$display)
    plots[[idx]] <- make_plot(mid, k$cohort, k$display, k$endpoint)
  }
}

grid <- plot_grid(plotlist=plots, ncol=4, nrow=4, align="hv")
# Row labels on left margin
row_labels <- plot_grid(
  plotlist=lapply(MODS, function(mid) {
    ggdraw() + draw_label(paste0("M", mid), angle=90, size=11, fontface="bold", colour="#0f2747")
  }), ncol=1, nrow=4)
# Col labels on top
col_labels <- plot_grid(
  plotlist=lapply(KEEP, function(k) {
    ggdraw() + draw_label(k$display, size=9, fontface="bold", colour="#333333")
  }), ncol=4, nrow=1)

titled <- plot_grid(
  ggdraw() + draw_label(
    "Non-M3 modules · mean-split KM (singscore) · 4 cohorts",
    x=0.02, hjust=0, size=12, fontface='bold', colour="#0f2747"),
  col_labels,
  plot_grid(row_labels, grid, ncol=2, rel_widths=c(0.03, 0.97)),
  ncol=1, rel_heights=c(0.035, 0.035, 0.93))

fp_pdf <- file.path(DIR, "KM_nonM3_mean_4x4.pdf")
fp_png <- file.path(DIR, "KM_nonM3_mean_4x4.png")
ggsave(fp_pdf, titled, width=14, height=12, limitsize=FALSE)
ggsave(fp_png, titled, width=14, height=12, dpi=150, limitsize=FALSE)
plog("[saved] ", fp_png)
plog("Done.")

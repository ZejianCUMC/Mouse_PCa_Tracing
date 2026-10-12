#!/usr/bin/env Rscript
# KM survival (singscore only) per-threshold across 10 human PCa cohorts.
# Reuses the same helper structure as km_All_FSU_noNC_min20_8modules.R but limits
# to singscore for speed.
#
# Loops T1..T4 sequentially; writes km_all_modules_singscore.csv into each subfolder.
suppressPackageStartupMessages({
  library(singscore); library(survival); library(data.table); library(dplyr)
})
plog <- function(...) cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n", sep="")

NPP53 <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
CACHE <- file.path(NPP53, "analysis/fitness_signatures_survival/cohort_cache_pc.rds")
SWEEP <- file.path(NPP53, "analysis/fitness_signatures_survival/threshold_sweep")

plog("Loading cohort cache ...")
COH <- readRDS(CACHE)
plog("  cohorts: ", paste(names(COH), collapse=", "))

score_singscore <- function(expr, genes) {
  p <- intersect(genes, rownames(expr)); if (length(p)<3) return(NULL)
  m <- as.matrix(expr); if (any(duplicated(colnames(m)))) m <- m[,!duplicated(colnames(m))]
  rk <- singscore::rankGenes(m)
  ss <- singscore::simpleScore(rk, upSet=p, centerScore=FALSE)
  setNames(ss$TotalScore, rownames(ss))
}

resolve_clin <- function(cohort_name, score) {
  co <- COH[[cohort_name]]
  eps <- if (cohort_name=="TCGA-PRAD") c("OS","PFI","DFI","DSS")
         else if (cohort_name %in% c("DKFZ-primary","MSKCC-primary")) "BCR"
         else "OS"
  out <- list()
  for (ep in eps) {
    if (cohort_name=="TCGA-PRAD") {
      df <- data.frame(sample_id=as.character(names(score)), score=as.numeric(score), stringsAsFactors=FALSE)
      co$cli$sample_id <- as.character(co$cli$sample_id)
      d <- merge(df, co$cli, by="sample_id")
      tc <- switch(ep, OS=c("OS.time","OS"), PFI=c("PFI.time","PFI"),
                   DFI=c("DFI.time","DFI"), DSS=c("DSS.time","DSS"))
    } else if (grepl("^WCDT", cohort_name)) {
      df <- data.frame(sample_id=as.character(names(score)), score=as.numeric(score), stringsAsFactors=FALSE)
      co$cli$sample_id <- as.character(co$cli$sample_id)
      d <- merge(df, co$cli, by="sample_id"); tc <- c("OS.time","OS")
    } else {
      patient <- co$map[names(score)]; ok <- !is.na(patient)
      tmp <- data.frame(sample_id=as.character(patient[ok]),
                        score=as.numeric(score)[ok], stringsAsFactors=FALSE) %>%
               group_by(sample_id) %>% summarise(score=mean(score), .groups="drop") %>% as.data.frame()
      co$cli$sample_id <- as.character(co$cli$sample_id)
      d <- merge(tmp, co$cli, by="sample_id")
      tc <- if (ep=="BCR") c("BCR.time","BCR") else c("OS.time","OS")
    }
    d$time <- d[[tc[1]]]; d$event <- d[[tc[2]]]
    d <- d[complete.cases(d$time,d$event,d$score) & d$time>0, ]
    if (nrow(d)>=10) out[[ep]] <- d
  }
  out
}

apply_split <- function(d, split_method) {
  if (split_method=="q1q4") {
    q1 <- quantile(d$score,0.25); q4 <- quantile(d$score,0.75)
    d$grp <- factor(ifelse(d$score<=q1,"low", ifelse(d$score>=q4,"high", NA)), levels=c("low","high"))
    d <- d[!is.na(d$grp),]
  } else {
    thr <- if (split_method=="mean") mean(d$score,na.rm=TRUE) else median(d$score,na.rm=TRUE)
    d$grp <- factor(ifelse(d$score>thr,"high","low"), levels=c("low","high"))
  }
  d
}
km_row <- function(d, cohort, mid, split, endpoint) {
  d <- apply_split(d, split)
  if (nrow(d)<10 || length(unique(d$grp))<2) return(NULL)
  sd_  <- try(survdiff(Surv(time,event)~grp, data=d), silent=TRUE)
  cx   <- try(coxph(Surv(time,event)~grp, data=d), silent=TRUE)
  p_lr <- if (inherits(sd_,"try-error")) NA_real_ else 1 - pchisq(sd_$chisq, df=length(sd_$n)-1)
  if (inherits(cx,"try-error") || is.null(coef(cx))) { HR<-NA_real_; p_cox<-NA_real_ }
  else { HR<-exp(coef(cx))[1]; p_cox<-summary(cx)$coefficients[1,5] }
  data.frame(cohort=cohort, module=mid, scoring_method="singscore", split=split, endpoint=endpoint,
             n_total=nrow(d), n_high=sum(d$grp=="high"), n_low=sum(d$grp=="low"),
             HR=HR, p_cox=p_cox, p=p_lr, stringsAsFactors=FALSE)
}

for (tid in c("T1","T2","T3","T4")) {
  plog("############## ", tid, " ##############")
  TDIR <- file.path(SWEEP, tid)
  if (!dir.exists(TDIR)) { plog("  skip — ", TDIR, " missing"); next }
  mod_files <- list.files(TDIR, pattern="^module_[0-9]+_human\\.txt$", full.names=TRUE)
  if (length(mod_files)==0) {
    plog("  no module_*_human.txt files found — skipping"); next
  }
  MODS <- list()
  for (fp in mod_files) {
    mid <- sub("^module_([0-9]+)_human\\.txt$", "M\\1", basename(fp))
    g <- unique(readLines(fp)); g <- g[nchar(g)>0]
    MODS[[mid]] <- g
    plog("  ", mid, ": ", length(g), " human genes")
  }

  all_rows <- list()
  COHORTS <- names(COH); SPLITS <- c("median","mean","q1q4")
  for (cohort in COHORTS) {
    plog(">>> ", cohort, " (", ncol(COH[[cohort]]$mat), " samples)")
    for (mid in names(MODS)) {
      score <- try(score_singscore(COH[[cohort]]$mat, MODS[[mid]]), silent=TRUE)
      if (inherits(score,"try-error") || is.null(score)) next
      clin_by_ep <- resolve_clin(cohort, score)
      for (ep in names(clin_by_ep)) for (sp in SPLITS) {
        row <- km_row(clin_by_ep[[ep]], cohort, mid, sp, ep)
        if (!is.null(row)) all_rows[[length(all_rows)+1]] <- row
      }
    }
  }

  if (length(all_rows)==0) {
    plog("  no KM rows produced for ", tid); next
  }
  res <- do.call(rbind, all_rows)
  res$p_bh_within <- ave(res$p, res$cohort, res$scoring_method, res$split, res$endpoint,
                          FUN=function(x) p.adjust(x, method="BH"))
  res <- res[order(res$p, na.last=TRUE), ]
  fwrite(res, file.path(TDIR, "km_all_modules_singscore.csv"))
  fwd <- res[res$HR>1 & is.finite(res$p) & res$p<0.05, ]
  plog("  ", tid, " → ", nrow(res), " KM rows | forward hits (HR>1, p<.05) = ", nrow(fwd))
}
plog("KM sweep done.")

#!/usr/bin/env Rscript
# =============================================================================
# opt5 step 00 — corrected human cohort matrices for singscore (audit, Zejian 2026-10-09: "carefully review all scripts and input
# format, then re-run the whole process with corrected input and scripts").
# singscore ranks genes WITHIN each sample, so inputs must be length-normalised expression (TPM / FPKM / microarray), never CPM or
# per-gene z-scores. All matrices: current HGNC protein-coding symbols x samples, log2 scale.
#   SU2C-mCRPC-union : polyA FPKM (preferred) + capture-only FPKM -> TPM -> log2(TPM + 1); OS from Abida 2019 PNAS Table S1 Sheet2;
#                      sample -> patient map (patient-mean score). (polyA vs capture: same-sample score r 0.91-0.96, no group shift.)
#   WCDT-JCO         : JCO 2018 patients, one biopsy each (baseline preferred), OS updated from EU-Uro 2019 (else JCO);
#                      original RAW counts -> TPM with effective gene lengths derived from the original Counts / TPM files -> log2(TPM+1)
#   WCDT-merged      : all patients of either cohort (JCO 2018 + EU-only), one biopsy each; sva::ComBat_seq on the RAW counts
#                      (Zejian: ComBat-seq needs raw counts), batch = cohort of origin -> TPM (same lengths) -> log2(TPM + 1)
#   TCGA-PRAD        : Xena STAR TPM, which is already log2(TPM + 1) (the old cache logged it twice), primary tumours, PFI (TCGA-CDR)
#   DKFZ-primary     : RPKM -> log2(RPKM + 1), BCR (Gerhauser 2018 cohort; cBioPortal prostate_dkfz_2018)
#   MSKCC-primary    : Agilent microarray (log2 intensities), Entrez -> symbol (org.Hs.eg.db), BCR (file's OS_* columns hold recurrence)
#   Michigan dropped : only per-gene z-scores are available (a within-sample rank of z-scores is not an expression rank).
# Symbols: an outdated / alias symbol is replaced by the current HGNC symbol when the mapping is unique (HGNC complete set 2026-09-30),
# BEFORE the protein-coding filter (HGNC locus_type "gene with protein product") — the old cache filtered first and lost renamed genes
# (e.g. ATP5A1 -> ATP5F1A). Duplicated symbols: row max (as the old cache). WCDT Ensembl features -> symbol: TPM summed.
# Output cohorts_corrected.rds: list(name = list(mat, cli = data.frame(id, time, event), map = c(sample = id), endpoint)),
# cohort_build_report.txt, wcdt_batch_pca.png. Env: rMAST (data.table, readxl, sva, org.Hs.eg.db). Deterministic.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(readxl); library(sva); library(org.Hs.eg.db); library(AnnotationDbi)})
S <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE))))
PA <- "/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas"
WCDT <- "/media/jrchristin/BackUp8T2/project/8-Jinqiu/NPSmad4/WCDT_metastatic_prostate_tumor"
HGNC <- "/media/jrchristin/Backup10TB1/agent/1-project/25-CGC/data/references/hgnc/hgnc_complete_set_20260930.txt"
REP <- c(); say <- function(...) { s <- paste0(...); REP <<- c(REP, s); cat(s, "\n") }

# ---- HGNC symbol harmonisation -------------------------------------------------------------------------------------------------
hg <- fread(HGNC, sep="\t", quote="\"", colClasses="character", select=c("symbol", "locus_type", "prev_symbol", "alias_symbol"))
CUR <- unique(hg$symbol); PC <- unique(hg[locus_type == "gene with protein product", symbol])
al <- rbindlist(lapply(c("prev_symbol", "alias_symbol"), function(col) hg[get(col) != "", .(old=unlist(strsplit(gsub("\"", "", get(col)), "\\|"))), by=symbol]))
al <- unique(al[!old %in% CUR & old != ""]); al <- al[, .(new=unique(symbol)), by=old][, if (.N == 1) .SD, by=old]
A2C <- setNames(al$new, al$old)
harm <- function(sym) { sym <- as.character(sym); i <- !(sym %in% CUR) & sym %in% names(A2C); sym[i] <- A2C[sym[i]]; sym }
finish <- function(m, sym, label, agg=c("max", "sum")) {          # harmonise -> protein-coding -> dedup ; m numeric matrix (rows = sym)
  agg <- match.arg(agg); s0 <- as.character(sym); s1 <- harm(s0); keep <- !is.na(s1) & s1 != "" & s1 %in% PC
  dt <- data.table(SYMBOL=s1[keep], m[keep, , drop=FALSE]); out <- if (agg == "max") dt[, lapply(.SD, max, na.rm=TRUE), by=SYMBOL] else dt[, lapply(.SD, sum, na.rm=TRUE), by=SYMBOL]
  r <- as.matrix(out[, -1, with=FALSE]); rownames(r) <- out$SYMBOL
  say(sprintf("  %-18s %6d input rows -> %6d protein-coding symbols (%d renamed via HGNC prev/alias; %d rows not protein-coding/unknown)",
              label, length(s0), nrow(r), sum(s1 != s0 & keep, na.rm=TRUE), sum(!keep)))
  r }
fpkm_to_tpm <- function(m) sweep(m, 2, colSums(m, na.rm=TRUE), "/") * 1e6
C <- list()

# ---- SU2C-mCRPC (Abida 2019) ---------------------------------------------------------------------------------------------------
say("SU2C-mCRPC-union")
SU <- file.path(PA, "prad_su2c_2019")
rd <- function(f) { x <- fread(file.path(SU, f)); x <- x[!is.na(Hugo_Symbol) & Hugo_Symbol != ""]; m <- as.matrix(x[, -1, with=FALSE]); list(m=m, s=x$Hugo_Symbol) }
pa <- rd("data_mrna_seq_fpkm_polya.txt"); ca <- rd("data_mrna_seq_fpkm_capture.txt")
mpa <- finish(pa$m, pa$s, "polyA FPKM"); mca <- finish(ca$m, ca$s, "capture FPKM")
g <- intersect(rownames(mpa), rownames(mca)); conly <- setdiff(colnames(mca), colnames(mpa))
su <- log2(fpkm_to_tpm(cbind(mpa[g, ], mca[g, conly, drop=FALSE])) + 1)      # TPM computed within the shared protein-coding universe
samp <- fread(file.path(SU, "data_clinical_sample.txt"), skip=4L); smap <- setNames(as.character(samp$PATIENT_ID), samp$SAMPLE_ID)
sx <- as.data.frame(read_excel(file.path(SU, "Abida2019_PNAS_TableS1.xlsx"), sheet="Sheet2")); colnames(sx)[c(1, 3, 4)] <- c("PATIENT_ID", "Death", "OS_months")
scli <- data.frame(id=as.character(sx$PATIENT_ID), time=as.numeric(sx$OS_months) * 30.44, event=as.integer(sx$Death)); scli <- scli[complete.cases(scli) & scli$time > 0, ]
C[["SU2C-mCRPC-union"]] <- list(mat=su, cli=scli, map=smap[colnames(su)], endpoint="OS")
say(sprintf("  %d samples (%d polyA + %d capture-only), %d genes; OS patients %d, with RNA %d", ncol(su), ncol(mpa), length(conly), nrow(su), nrow(scli), sum(scli$id %in% smap[colnames(su)])))

# ---- WCDT (raw counts; one biopsy per patient; OS JCO updated from EU 2019) ---------------------------------------------------------
say("WCDT")
cn <- fread(file.path(WCDT, "Counts_Genes.txt")); tp <- fread(file.path(WCDT, "Counts_TPM_Genes.txt"))
smp <- intersect(setdiff(names(cn), c("FEATURE_ID", "seqid", "gene_type", "gene_name")), names(tp))
stopifnot(identical(cn$FEATURE_ID, tp$FEATURE_ID))
Mc <- as.matrix(cn[, ..smp]); Mt <- as.matrix(tp[, ..smp]); rownames(Mc) <- rownames(Mt) <- cn$FEATURE_ID
stopifnot(all(Mc == round(Mc)))                                               # raw integer counts (ComBat-seq input)
# effective length per feature: counts / TPM = L_f * K_s  ->  two-way median polish on the log ratio (entries with counts >= 5)
lr <- log(Mc / Mt); lr[!(Mc >= 5 & Mt > 0)] <- NA; a <- rep(0, nrow(lr)); b <- rep(0, ncol(lr))
for (it in 1:4) { b <- apply(lr - a, 2, median, na.rm=TRUE); a <- apply(sweep(lr, 2, b), 1, median, na.rm=TRUE) }
res <- sweep(lr - a, 2, b); say(sprintf("  effective lengths: %d features estimated; residual SD of log(count/TPM) %.4f (0 = constant lengths)", sum(is.finite(a)), sd(res, na.rm=TRUE)))
L <- exp(a); names(L) <- rownames(Mc)
sym <- harm(cn$gene_name); fe <- !is.na(sym) & sym %in% PC & is.finite(L)
tpm_of <- function(counts) { x <- counts[fe, , drop=FALSE] / L[fe]; x <- sweep(x, 2, colSums(x), "/") * 1e6
  dt <- data.table(SYMBOL=sym[fe], x)[, lapply(.SD, sum), by=SYMBOL]; r <- as.matrix(dt[, -1, with=FALSE]); rownames(r) <- dt$SYMBOL; log2(r + 1) }
chk <- tpm_of(Mc); orig <- { x <- Mt[fe, ]; x <- sweep(x, 2, colSums(x), "/") * 1e6; dt <- data.table(SYMBOL=sym[fe], x)[, lapply(.SD, sum), by=SYMBOL]; r <- as.matrix(dt[, -1, with=FALSE]); rownames(r) <- dt$SYMBOL; log2(r + 1) }
say(sprintf("  TPM rebuilt from raw counts vs original TPM file: median per-sample Spearman %.5f, max |diff| %.3f (log2)",
            median(sapply(seq_len(ncol(chk)), function(i) cor(chk[, i], orig[rownames(chk), i], method="spearman"))), max(abs(chk - orig[rownames(chk), ]))))
pid <- function(s) sub("-(BL|BL2|2-BL|PRO|PRO2|PRO3|RP)$", "", s); rk <- function(s) match(sub("^DTB-[0-9]+-", "", s), c("BL", "BL2", "2-BL", "PRO", "PRO2", "PRO3", "RP"))
j <- fread(file.path(WCDT, "clinical_metadata/2021_11_21_JCO_2018.txt")); setnames(j, "OS_event_(dead)", "ev")
j <- j[, .(sample_id, patient=pid(sample_id), tJ=suppressWarnings(as.numeric(OS_duration)), dJ=as.integer(grepl("Death", ev, ignore.case=TRUE)))]
e <- fread(file.path(WCDT, "clinical_metadata/2021_05_13_european_urology_2019.txt"))[, .(sample_id, patient=patient_id, tE=suppressWarnings(as.numeric(OS.mCRPC)), dE=as.integer(event.Dead))]
bi <- unique(rbind(j[, .(sample_id, patient)], e[, .(sample_id, patient)]))[sample_id %in% smp]; bi <- bi[order(patient, rk(sample_id))][!duplicated(patient)]
jp <- j[!is.na(tJ) & tJ > 0][!duplicated(patient)]; ep <- e[!is.na(tE) & tE > 0][!duplicated(patient)]
cl <- merge(merge(bi, jp[, .(patient, tJ, dJ)], by="patient", all.x=TRUE), ep[, .(patient, tE, dE)], by="patient", all.x=TRUE)
cl[, `:=`(time=fifelse(!is.na(tE), tE, tJ), event=fifelse(!is.na(tE), dE, dJ), in_JCO=patient %in% j$patient, in_EU=patient %in% e$patient)]
cl <- cl[!is.na(time) & time > 0]; cl[, batch := fifelse(in_JCO, "JCO2018", "EU2019only")]
A <- cl[in_JCO == TRUE]; mA <- tpm_of(Mc[, A$sample_id])
C[["WCDT-JCO"]] <- list(mat=mA, cli=data.frame(id=A$patient, time=A$time, event=A$event), map=setNames(A$patient, A$sample_id), endpoint="OS")
say(sprintf("  WCDT-JCO: %d patients, %d deaths (OS from EU 2019: %d, JCO: %d), %d genes", nrow(A), sum(A$event), sum(!is.na(A$tE)), sum(is.na(A$tE)), nrow(mA)))
Mb <- Mc[fe, cl$sample_id]; Mb <- Mb[rowSums(Mb > 0) >= 10, ]
Mcs <- ComBat_seq(Mb, batch=cl$batch); fe_b <- rownames(Mc) %in% rownames(Mcs)
tpm_b <- { x <- Mcs / L[rownames(Mcs)]; x <- sweep(x, 2, colSums(x), "/") * 1e6; dt <- data.table(SYMBOL=sym[match(rownames(Mcs), rownames(Mc))], x)[, lapply(.SD, sum), by=SYMBOL]
           r <- as.matrix(dt[, -1, with=FALSE]); rownames(r) <- dt$SYMBOL; log2(r + 1) }
C[["WCDT-merged"]] <- list(mat=tpm_b, cli=data.frame(id=cl$patient, time=cl$time, event=cl$event), map=setNames(cl$patient, cl$sample_id), endpoint="OS")
say(sprintf("  WCDT-merged (ComBat-seq on raw counts -> TPM): %d patients (%d JCO incl. %d also in EU; %d EU-only), %d deaths, %d genes",
            nrow(cl), sum(cl$in_JCO), sum(cl$in_JCO & cl$in_EU), sum(!cl$in_JCO), sum(cl$event), nrow(tpm_b)))
pcs <- function(x) { v <- apply(x, 1, var); x <- x[order(-v)[1:2000], ]; p <- prcomp(t(x), scale.=TRUE); p$x[, 1:2] }
before <- { x <- Mb / L[rownames(Mb)]; x <- sweep(x, 2, colSums(x), "/") * 1e6; log2(x + 1) }; after <- { x <- Mcs / L[rownames(Mcs)]; x <- sweep(x, 2, colSums(x), "/") * 1e6; log2(x + 1) }
r2 <- function(p) sapply(1:2, function(k) round(summary(lm(p[, k] ~ cl$batch))$r.squared, 3)); p0 <- pcs(before); p1 <- pcs(after)
say(sprintf("  batch R2 on PC1/PC2 before ComBat-seq %s, after %s", paste(r2(p0), collapse="/"), paste(r2(p1), collapse="/")))
png(file.path(S, "wcdt_batch_pca.png"), width=1300, height=600, res=120); par(mfrow=c(1, 2)); colb <- ifelse(cl$batch == "JCO2018", "#4a7fbf", "#c44e52")
plot(p0, col=colb, pch=19, main=sprintf("log2 TPM before ComBat-seq (batch R2 PC1 %.3f)", r2(p0)[1]), cex.main=.85); legend("topright", c("JCO 2018 patients", "EU 2019-only patients"), col=c("#4a7fbf", "#c44e52"), pch=19, cex=.8)
plot(p1, col=colb, pch=19, main=sprintf("after ComBat-seq on raw counts (batch R2 PC1 %.3f)", r2(p1)[1]), cex.main=.85); invisible(dev.off())

# ---- TCGA-PRAD (validation) ---------------------------------------------------------------------------------------------------
say("TCGA-PRAD")
tc <- fread(cmd=paste("zcat", file.path(PA, "TCGA/TCGA-PRAD.star_tpm.tsv.gz"))); pm <- fread(file.path(PA, "TCGA/gencode.v36.annotation.gtf.gene.probemap"))
tsym <- setNames(pm$gene, pm$id)[tc$Ensembl_ID]; tm <- as.matrix(tc[, -1, with=FALSE])            # already log2(TPM + 1)
tm <- finish(tm, tsym, "TCGA log2TPM"); tm <- tm[, grepl("-01[AB]$", colnames(tm))]; pt <- substr(colnames(tm), 1, 12)
cdr <- suppressWarnings(as.data.frame(read_excel(file.path(PA, "TCGA/TCGA-CDR-SupplementalTableS1.xlsx"), sheet=1))); cdr <- cdr[cdr$type == "PRAD", ]
tcli <- data.frame(id=as.character(cdr$bcr_patient_barcode), time=as.numeric(cdr$PFI.time), event=as.integer(cdr$PFI)); tcli <- tcli[complete.cases(tcli) & tcli$time > 0, ]
C[["TCGA-PRAD"]] <- list(mat=tm, cli=tcli, map=setNames(pt, colnames(tm)), endpoint="PFI")
say(sprintf("  %d tumour samples (%d patients), PFI patients %d", ncol(tm), length(unique(pt)), nrow(tcli)))

# ---- DKFZ-primary (validation) ------------------------------------------------------------------------------------------------
say("DKFZ-primary")
DK <- file.path(PA, "prostate_dkfz_2018"); dx <- fread(file.path(DK, "data_mrna_seq_rpkm.txt")); dx <- dx[!is.na(Hugo_Symbol) & Hugo_Symbol != ""]
dm <- log2(finish(as.matrix(dx[, -(1:2), with=FALSE]), dx$Hugo_Symbol, "DKFZ RPKM") + 1)
dsamp <- fread(file.path(DK, "data_clinical_sample.txt"), skip=4L); dpat <- fread(file.path(DK, "data_clinical_patient.txt"), skip=4L)
dcli <- data.frame(id=as.character(dpat$PATIENT_ID), time=as.numeric(dpat$TIME_FROM_SURGERY_TO_BCR_LASTFU) * 30.44, event=suppressWarnings(as.integer(dpat$BCR_STATUS)))
dcli <- dcli[complete.cases(dcli) & dcli$time > 0, ]
C[["DKFZ-primary"]] <- list(mat=dm, cli=dcli, map=setNames(as.character(dsamp$PATIENT_ID), dsamp$SAMPLE_ID)[colnames(dm)], endpoint="BCR")
say(sprintf("  %d samples, BCR patients with RNA %d", ncol(dm), sum(dcli$id %in% C[["DKFZ-primary"]]$map)))

# ---- MSKCC-primary (validation) -----------------------------------------------------------------------------------------------
say("MSKCC-primary")
MS <- file.path(PA, "prad_mskcc"); mx <- fread(file.path(MS, "data_mrna_agilent_microarray.txt"))
msym <- unname(mapIds(org.Hs.eg.db, keys=as.character(mx$Entrez_Gene_Id), column="SYMBOL", keytype="ENTREZID", multiVals="first"))
mm <- finish(as.matrix(mx[, -1, with=FALSE]), msym, "MSKCC array")
mpat <- fread(file.path(MS, "data_clinical_patient.txt"), sep=","); msamp <- fread(file.path(MS, "data_clinical_sample.txt"), skip="SAMPLE_ID")
mcli <- data.frame(id=as.character(mpat$PATIENT_ID), time=suppressWarnings(as.numeric(mpat$OS_MONTHS)) * 30.44,
                   event=ifelse(grepl("^1:Recurred", mpat$OS_STATUS), 1L, ifelse(grepl("^0:DiseaseFree", mpat$OS_STATUS), 0L, NA)))
mcli <- mcli[complete.cases(mcli) & mcli$time > 0, ]
C[["MSKCC-primary"]] <- list(mat=mm, cli=mcli, map=setNames(as.character(msamp$PATIENT_ID), msamp$SAMPLE_ID)[colnames(mm)], endpoint="BCR")
say(sprintf("  %d samples, BCR patients with RNA %d (events %d)", ncol(mm), sum(mcli$id %in% C[["MSKCC-primary"]]$map), sum(mcli$event[mcli$id %in% C[["MSKCC-primary"]]$map])))

for (n in names(C)) stopifnot(!anyNA(C[[n]]$mat), all(is.finite(C[[n]]$mat)), !anyDuplicated(rownames(C[[n]]$mat)))
saveRDS(C, file.path(S, "cohorts_corrected.rds")); writeLines(REP, file.path(S, "cohort_build_report.txt")); say("saved cohorts_corrected.rds")

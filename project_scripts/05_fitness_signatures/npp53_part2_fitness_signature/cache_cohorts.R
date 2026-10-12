#!/usr/bin/env Rscript
# One-shot: load 5 public prostate cohorts + cache as .rds using data.table::fread + rowsum-style dedup.
# Standard cohort names:
#   TCGA-PRAD        · Cancer Genome Atlas prostate adenocarcinoma (Liu 2018 CDR)
#   SU2C-mCRPC       · SU2C-PCF East Coast Dream Team mCRPC (Abida 2019 PNAS)
#                      MERGED library: prefer polyA, fallback to capture for capture-only samples (union 328)
#                      FPKM → TPM per sample (TPM = FPKM_i / sum(FPKM_i) * 1e6)
#   Michigan-mCRPC   · Grasso 2012 Nature mCRPC (Michigan / Chinnaiyan)
#   DKFZ-primary     · Fraser 2018 primary PCa (BCR)
#   MSKCC-primary    · Taylor 2010 primary PCa (BCR)
suppressPackageStartupMessages({
  library(data.table); library(readxl); library(org.Hs.eg.db); library(AnnotationDbi)
})
# SU2C OS source: Abida 2019 PNAS Table S1 Sheet2 (128 patients with survival)
# — richer than cBioPortal cli OS_STATUS/OS_MONTHS (which is also 128 but less clean)
SU2C_OS_SUPP <- "Abida2019_PNAS_TableS1.xlsx"
plog <- function(...) { cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n", sep=""); flush.console() }

PA    <- "/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas"
NPP53 <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
OUT   <- file.path(NPP53, "analysis/fitness_signatures_survival/cohort_cache.rds")

TCGA  <- file.path(PA, "TCGA"); SU2C <- file.path(PA, "prad_su2c_2019")
MICH  <- file.path(PA, "prad_mich"); DKFZ <- file.path(PA, "prostate_dkfz_2018")
MSKCC <- file.path(PA, "prad_mskcc")
WCDT  <- "/media/jrchristin/BackUp8T2/project/8-Jinqiu/NPSmad4/WCDT_metastatic_prostate_tumor"

# per-symbol max on a data.table — 10-50× faster than dplyr summarise_all(max)
dedup_max <- function(dt, sym_col) {
  setnames(dt, sym_col, "SYMBOL")
  dt <- dt[!is.na(SYMBOL) & SYMBOL != ""]
  if ("Entrez_Gene_Id" %in% names(dt)) dt[, Entrez_Gene_Id := NULL]
  agg <- dt[, lapply(.SD, max, na.rm=TRUE), by=SYMBOL]
  syms <- agg$SYMBOL; agg[, SYMBOL := NULL]
  m <- as.matrix(agg); rownames(m) <- syms
  m
}
# FPKM → TPM: per-sample normalization to 1M
fpkm_to_tpm <- function(m) {
  col_sums <- colSums(m, na.rm=TRUE)
  sweep(m, 2, col_sums, FUN="/") * 1e6
}

# ---- 1. TCGA-PRAD ----
plog("[1/5] TCGA-PRAD")
t0 <- Sys.time()
cdr <- suppressWarnings(as.data.frame(read_excel(file.path(TCGA, "TCGA-CDR-SupplementalTableS1.xlsx"), sheet=1)))
prad <- cdr[cdr$type == "PRAD", ]
prad$sample_id <- as.character(prad$bcr_patient_barcode)
tpm <- fread(cmd = paste("zcat", file.path(TCGA, "TCGA-PRAD.star_tpm.tsv.gz")), showProgress=FALSE)
probemap <- fread(file.path(TCGA, "gencode.v36.annotation.gtf.gene.probemap"))
tpm[, SYMBOL := setNames(probemap$gene, probemap$id)[Ensembl_ID]]
tpm[, Ensembl_ID := NULL]
m <- dedup_max(tpm, "SYMBOL")
tumor <- grepl("-01[AB]$", colnames(m)); m <- m[, tumor]
colnames(m) <- substr(colnames(m), 1, 12)
if (any(duplicated(colnames(m)))) {
  m <- vapply(unique(colnames(m)), function(p) {
                cols <- which(colnames(m) == p)
                if (length(cols)==1) m[,cols] else apply(m[,cols], 1, max)
              }, numeric(nrow(m)))
}
tpm_tumor <- log2(m + 1)
plog("  TCGA-PRAD: ", nrow(tpm_tumor), " × ", ncol(tpm_tumor), "  (", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "s)")

# ---- 2. SU2C-mCRPC — THREE cohorts: polyA + capture + UNION (polyA preferred fallback capture) ----
plog("[2/5] SU2C-mCRPC (polyA + capture + union, FPKM→TPM, OS from Abida 2019 Sheet2)")
t0 <- Sys.time()
polya <- fread(file.path(SU2C, "data_mrna_seq_fpkm_polya.txt"), showProgress=FALSE)
capt  <- fread(file.path(SU2C, "data_mrna_seq_fpkm_capture.txt"), showProgress=FALSE)
su2c_polya   <- log2(fpkm_to_tpm(dedup_max(copy(polya), "Hugo_Symbol")) + 1)
su2c_capture <- log2(fpkm_to_tpm(dedup_max(copy(capt),  "Hugo_Symbol")) + 1)
# Build UNION: polyA preferred, fallback to capture for capture-only samples
common_genes <- intersect(rownames(su2c_polya), rownames(su2c_capture))
capture_only <- setdiff(colnames(su2c_capture), colnames(su2c_polya))
su2c_union <- cbind(su2c_polya[common_genes, ],
                    su2c_capture[common_genes, capture_only, drop=FALSE])
plog("  SU2C-polyA:    ", nrow(su2c_polya),   " × ", ncol(su2c_polya),   " samples")
plog("  SU2C-capture:  ", nrow(su2c_capture), " × ", ncol(su2c_capture), " samples")
plog("  SU2C-union:    ", nrow(su2c_union),   " × ", ncol(su2c_union), " samples (polyA=", ncol(su2c_polya),
     " + capture-only=", length(capture_only), ")")

# Clinical: use Abida 2019 PNAS Table S1 Sheet2 (128 patients with survival)
su2c_samp <- fread(file.path(SU2C, "data_clinical_sample.txt"), skip=4L, showProgress=FALSE)
su2c_map  <- setNames(as.character(su2c_samp$PATIENT_ID), su2c_samp$SAMPLE_ID)
supp <- as.data.frame(read_excel(file.path(SU2C, SU2C_OS_SUPP), sheet="Sheet2"))
colnames(supp)[c(1,3,4)] <- c("PATIENT_ID","Death","OS_ARSI_months")
su2c_cli <- data.table(
  sample_id  = as.character(supp$PATIENT_ID),
  PATIENT_ID = as.character(supp$PATIENT_ID),
  OS         = as.integer(supp$Death),
  OS.time    = as.numeric(supp$OS_ARSI_months) * 30.44
)
su2c_cli <- su2c_cli[!is.na(OS) & !is.na(OS.time) & OS.time > 0]
plog("  clinical (Abida Sheet2): ", nrow(su2c_cli), " patients with OS")
plog("  (", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "s)")

# ---- 3. Michigan-mCRPC (Grasso 2012, RPKM z-scores) ----
plog("[3/5] Michigan-mCRPC (Grasso 2012)")
t0 <- Sys.time()
mich_expr <- fread(file.path(MICH, "data_mrna_seq_rpkm_zscores_ref_diploid_samples.txt"), showProgress=FALSE)
mich_mat  <- dedup_max(copy(mich_expr), "Hugo_Symbol")
mich_cli  <- fread(file.path(MICH, "data_clinical_patient.txt"), skip=4L, showProgress=FALSE)
mich_samp <- fread(file.path(MICH, "data_clinical_sample.txt"), skip=4L, showProgress=FALSE)
mich_map  <- setNames(as.character(mich_samp$PATIENT_ID), mich_samp$SAMPLE_ID)
mich_cli[, OS := as.integer(grepl("DECEASED|1:", OS_STATUS))]
mich_cli[grepl("LIVING|0:", OS_STATUS), OS := 0L]
mich_cli[, OS.time := as.numeric(OS_MONTHS) * 30.44]
mich_cli[, sample_id := as.character(PATIENT_ID)]
plog("  Michigan-mCRPC: ", nrow(mich_mat), " × ", ncol(mich_mat), "  (", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "s)")

# ---- 4. DKFZ Fraser (BCR) ----
plog("[4/5] DKFZ-primary (Fraser 2018)")
t0 <- Sys.time()
dkfz_expr <- fread(file.path(DKFZ, "data_mrna_seq_rpkm.txt"), showProgress=FALSE)
dkfz_mat  <- log2(dedup_max(copy(dkfz_expr), "Hugo_Symbol") + 1)
dkfz_cli  <- fread(file.path(DKFZ, "data_clinical_patient.txt"), skip=4L, showProgress=FALSE)
dkfz_samp <- fread(file.path(DKFZ, "data_clinical_sample.txt"), skip=4L, showProgress=FALSE)
dkfz_map  <- setNames(as.character(dkfz_samp$PATIENT_ID), dkfz_samp$SAMPLE_ID)
dkfz_cli[, BCR := suppressWarnings(as.numeric(BCR_STATUS))]
dkfz_cli[, BCR.time := as.numeric(TIME_FROM_SURGERY_TO_BCR_LASTFU) * 30.44]
dkfz_cli[, sample_id := as.character(PATIENT_ID)]
plog("  DKFZ-primary: ", nrow(dkfz_mat), " × ", ncol(dkfz_mat), "  (", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "s)")

# ---- 5. MSKCC Taylor (BCR) ----
plog("[5/5] MSKCC-primary (Taylor 2010)")
t0 <- Sys.time()
suppressPackageStartupMessages(library(org.Hs.eg.db))
mskcc_expr <- fread(file.path(MSKCC, "data_mrna_agilent_microarray.txt"), showProgress=FALSE)
entrez <- as.character(mskcc_expr$Entrez_Gene_Id)
sym <- unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys=entrez, column="SYMBOL", keytype="ENTREZID", multiVals="first"))
mskcc_expr[, SYMBOL := sym]; mskcc_expr[, Entrez_Gene_Id := NULL]
mskcc_mat  <- dedup_max(copy(mskcc_expr), "SYMBOL")
mskcc_cli  <- fread(file.path(MSKCC, "data_clinical_patient.txt"), sep=",", showProgress=FALSE)
mskcc_samp <- fread(file.path(MSKCC, "data_clinical_sample.txt"), skip="SAMPLE_ID", showProgress=FALSE)
mskcc_map  <- setNames(as.character(mskcc_samp$PATIENT_ID), mskcc_samp$SAMPLE_ID)
mskcc_cli[, BCR := as.integer(grepl("1:Recurred|^1$", OS_STATUS))]
mskcc_cli[grepl("DiseaseFree|^0:|^0$", OS_STATUS), BCR := 0L]
mskcc_cli[, BCR.time := suppressWarnings(as.numeric(OS_MONTHS)) * 30.44]
mskcc_cli[, sample_id := as.character(PATIENT_ID)]
plog("  MSKCC-primary: ", nrow(mskcc_mat), " × ", ncol(mskcc_mat), "  (", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "s)")

# ---- 6. WCDT-mCRPC (West Coast Dream Team) — Aggarwal JCO 2018 OS + Quigley Cell 2018 EU-Uro 2019 OS ----
plog("[6] WCDT-mCRPC (JCO 2018 + EU-Uro 2019 clinical, TPM expression)")
t0 <- Sys.time()
wcdt_expr <- fread(file.path(WCDT, "Counts_TPM_Genes.txt"), showProgress=FALSE)
# File layout: col 1 = FEATURE_ID (Ensembl); cols 2..(N-1) = samples (TPM); LAST col = gene symbol
sym_col <- names(wcdt_expr)[ncol(wcdt_expr)]
wcdt_expr[, FEATURE_ID := NULL]  # drop Ensembl
setnames(wcdt_expr, sym_col, "SYMBOL")
wcdt_mat  <- log2(dedup_max(copy(wcdt_expr), "SYMBOL") + 1)   # TPM → log2(TPM+1)
plog("  WCDT TPM matrix: ", nrow(wcdt_mat), " × ", ncol(wcdt_mat), " samples")

# JCO 2018 clinical (n=118, OS_duration + OS_event_(dead))
jco <- fread(file.path(WCDT, "clinical_metadata/2021_11_21_JCO_2018.txt"))
setnames(jco, "OS_event_(dead)", "OS_event")
jco[, OS      := as.integer(grepl("Death", OS_event, ignore.case=TRUE))]
jco[, OS.time := suppressWarnings(as.numeric(OS_duration))]     # days
jco[, sample_id := as.character(sample_id)]
jco_cli <- jco[!is.na(OS) & !is.na(OS.time) & OS.time > 0,
               .(sample_id, PATIENT_ID = as.character(sample_id), OS, OS.time)]
plog("  JCO 2018 clinical: ", nrow(jco_cli), " samples with valid OS")

# EU-Uro 2019 (n=100, OS.mCRPC + event.Dead)
eur <- fread(file.path(WCDT, "clinical_metadata/2021_05_13_european_urology_2019.txt"))
eur[, OS      := as.integer(event.Dead)]
eur[, OS.time := suppressWarnings(as.numeric(OS.mCRPC))]        # days
eur[, sample_id := as.character(sample_id)]
eur_cli <- eur[!is.na(OS) & !is.na(OS.time) & OS.time > 0,
               .(sample_id, PATIENT_ID = as.character(sample_id), OS, OS.time)]
plog("  EU-Uro 2019 clinical: ", nrow(eur_cli), " samples with valid OS")

# Union (JCO ∪ EU-Uro) — prefer EU-Uro OS.mCRPC when both present (published later)
union_ids <- union(jco_cli$sample_id, eur_cli$sample_id)
wcdt_union_cli <- rbindlist(list(
  eur_cli,
  jco_cli[!sample_id %in% eur_cli$sample_id]
), use.names=TRUE)
plog("  WCDT union clinical: ", nrow(wcdt_union_cli), " samples")

# WCDT map = identity (SAMPLE_ID = "DTB-xxx-BL", also serves as patient key at biopsy level)
wcdt_map <- setNames(colnames(wcdt_mat), colnames(wcdt_mat))
plog("  (", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "s)")

# ---- assemble + save ----
cohorts <- list(
  "TCGA-PRAD"          = list(mat=tpm_tumor,    cli=as.data.frame(prad),      map=NULL,      eps=c("OS","PFI","DFI","DSS")),
  "SU2C-mCRPC-polyA"   = list(mat=su2c_polya,   cli=as.data.frame(su2c_cli),  map=su2c_map,  eps="OS"),
  "SU2C-mCRPC-capture" = list(mat=su2c_capture, cli=as.data.frame(su2c_cli),  map=su2c_map,  eps="OS"),
  "SU2C-mCRPC-union"   = list(mat=su2c_union,   cli=as.data.frame(su2c_cli),  map=su2c_map,  eps="OS"),
  "Michigan-mCRPC"     = list(mat=mich_mat,     cli=as.data.frame(mich_cli),  map=mich_map,  eps="OS"),
  "DKFZ-primary"       = list(mat=dkfz_mat,     cli=as.data.frame(dkfz_cli),  map=dkfz_map,  eps="BCR"),
  "MSKCC-primary"      = list(mat=mskcc_mat,    cli=as.data.frame(mskcc_cli), map=mskcc_map, eps="BCR"),
  "WCDT-mCRPC-JCO"     = list(mat=wcdt_mat,     cli=as.data.frame(jco_cli),   map=wcdt_map,  eps="OS"),
  "WCDT-mCRPC-EU"      = list(mat=wcdt_mat,     cli=as.data.frame(eur_cli),   map=wcdt_map,  eps="OS"),
  "WCDT-mCRPC-union"   = list(mat=wcdt_mat,     cli=as.data.frame(wcdt_union_cli), map=wcdt_map, eps="OS")
)
saveRDS(cohorts, OUT)
plog("[saved] ", OUT, " (", format(file.size(OUT)/1e6, digits=3), " MB)")

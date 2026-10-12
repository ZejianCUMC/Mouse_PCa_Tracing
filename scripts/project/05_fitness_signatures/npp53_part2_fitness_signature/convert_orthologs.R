#!/usr/bin/env Rscript
# Convert Yang-2022 mouse FSU/FSD to human via MGI HOM_MouseHumanSequence.rpt (canonical 1:1 orthology).
# Writes: FitnessPosSignature_human.txt, FitnessNegSignature_human.txt + a mapping report.
suppressPackageStartupMessages({ library(data.table) })
plog <- function(...) { cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n", sep="") }

HOM_FP <- "/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/22Science_Tang_CRPC_ATAC/resources/orthology/HOM_MouseHumanSequence.rpt"
YANG   <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/FitnessSignature"

plog("Loading MGI HOM table ...")
hom <- fread(HOM_FP, sep="\t", quote="")
mouse <- hom[`Common Organism Name` == "mouse, laboratory", .(DBClassKey=`DB Class Key`, mouse_symbol=Symbol)]
human <- hom[`Common Organism Name` == "human",             .(DBClassKey=`DB Class Key`, human_symbol=Symbol)]
mh    <- merge(mouse, human, by="DBClassKey", allow.cartesian=TRUE)
plog("  ", nrow(mh), " mouse↔human ortholog pairs (", uniqueN(mh$mouse_symbol), " unique mouse, ",
     uniqueN(mh$human_symbol), " unique human)")

# Build mouse_symbol → set(human_symbol) map — take FIRST match per mouse gene (1:1)
mh_first <- mh[order(mouse_symbol, human_symbol), .(human_symbol = human_symbol[1]), by=mouse_symbol]
setkey(mh_first, mouse_symbol)

convert <- function(mouse_genes, name) {
  fp_out    <- file.path(YANG, sprintf("yang2022_%s_human.txt", name))
  fp_report <- file.path(YANG, sprintf("yang2022_%s_ortholog_mapping.tsv", name))
  x <- data.table(mouse_symbol = mouse_genes)
  y <- mh_first[x, on = "mouse_symbol"]
  n_total  <- length(mouse_genes)
  n_mapped <- sum(!is.na(y$human_symbol))
  writeLines(unique(y$human_symbol[!is.na(y$human_symbol)]), fp_out)
  fwrite(y, fp_report, sep="\t")
  # Compare with toupper fallback
  n_toupper <- sum(toupper(mouse_genes) %in% mh$human_symbol)
  plog(sprintf("  %s: %d mouse genes → %d human orthologs via MGI (%.1f%% mapped)",
               name, n_total, n_mapped, 100*n_mapped/n_total))
  plog(sprintf("     (toupper() naïve would give: %d — MGI recovers +%d genes)",
               n_toupper, n_mapped - n_toupper))
}

FSU_m <- readLines(file.path(YANG, "FitnessPosSignature_mouse.txt")); FSU_m <- FSU_m[nchar(FSU_m)>0]
FSD_m <- readLines(file.path(YANG, "FitnessNegSignature_mouse.txt")); FSD_m <- FSD_m[nchar(FSD_m)>0]
plog("Input: FSU=", length(FSU_m), " mouse genes, FSD=", length(FSD_m))

convert(FSU_m, "FSU")
convert(FSD_m, "FSD")
plog("Done. Outputs in ", YANG)

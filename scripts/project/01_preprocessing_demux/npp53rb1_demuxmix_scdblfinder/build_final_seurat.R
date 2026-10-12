suppressPackageStartupMessages(library(Seurat))
BASE <- "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing"
FO <- file.path(BASE,"final_objs"); libs <- c("MJZ025","MJZ026","MJZ027","MJZ028")
objs <- list()
for(lib in libs){
  m <- Read10X(file.path(FO,paste0(lib,"_mtx")))
  meta <- read.csv(file.path(FO,paste0(lib,"_final_meta.csv"))); rownames(meta) <- meta$cell
  so <- CreateSeuratObject(counts=m, project=lib)
  so <- AddMetaData(so, meta[colnames(so),])
  so <- PercentageFeatureSet(so, pattern="^mt-", col.name="percent.mt")
  saveRDS(so, file.path(FO,paste0(lib,"_final_singlet.RDS")))
  objs[[lib]] <- so
  cat(lib, ncol(so), "cells |", paste(capture.output(print(table(so$sample_id))),collapse=" "), "\n")
}
merged <- merge(objs[[1]], y=objs[-1], add.cell.ids=libs)
rb1 <- subset(merged, subset = route_project=="Npp53Rb1_Tracing")
saveRDS(rb1, file.path(FO,"Npp53Rb1_merged_singlet.RDS"))
orth <- subset(merged, subset = sample_id=="JZ282")
saveRDS(orth, file.path(FO,"JZ282_orthotopic_forNpp53.RDS"))
write.csv(merged@meta.data, file.path(FO,"Npp53Rb1_all_final_meta.csv"))
cat("\n=== FINAL ===\n")
cat("Npp53Rb1 merged (Np53Rb1 model, excl JZ282):", ncol(rb1), "cells\n")
cat("JZ282 orthotopic (-> Npp53 project):", ncol(orth), "cells\n")
cat("per-sample final counts:\n"); print(table(merged$sample_id, merged$treatment))

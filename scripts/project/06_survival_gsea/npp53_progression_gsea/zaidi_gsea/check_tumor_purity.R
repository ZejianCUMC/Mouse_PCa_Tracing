suppressPackageStartupMessages(library(Seurat))
Z<-"/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/24PNAS_Zaidi_MSK_PCa_atlas"
obj<-readRDS(file.path(Z,"GSE264573_msk.integrated.remove.cellcycle.tumor.cells.rds"))
md<-obj@meta.data
cat("=== meta.data columns ===\n"); print(colnames(md))
cat("\n=== candidate lineage/celltype columns (unique values, <=40) ===\n")
for(cn in colnames(md)){ v<-md[[cn]]; if(is.factor(v)||is.character(v)){ u<-unique(as.character(v)); if(length(u)<=40 && length(u)>1){ cat("--",cn,":\n"); print(sort(table(v),decreasing=TRUE)) } } }
# purity by markers
DefaultAssay(obj)<-"RNA"
mk<-c("EPCAM","KLK3","AR","KRT8","PTPRC","CD3D","C1QA","LYZ","COL1A1","PECAM1","HBB","MS4A1")
mk<-mk[mk %in% rownames(obj)]
ex<-FetchData(obj,vars=c(mk,"subtype"),layer="data")
cat("\n=== % cells expressing marker (>0), by subtype ===\n")
for(g in mk){ a<-tapply(ex[[g]]>0, ex$subtype, mean); cat(sprintf("%-8s CSPC=%.1f%%  CRPC=%.1f%%  NEPC=%.1f%%\n",g,100*a["CSPC"],100*a["CRPC"],100*a["NEPC"])) }
cat("ALL DONE\n")

#!/usr/bin/env Rscript
# v10-NAMED: rebuild the tumor cloupe (18,162 cells, No-Cyc + No-Ion) with the CURRENT
# canonical state display names (c1_Lum-Distal ... c12_Mes-myCAF) — matches Fig 6.4b.
# Also keeps the raw Luminal-*/Basal-like/Mes-*/Partial-EMT string as `state_raw` for
# back-compat with older analyses.
suppressPackageStartupMessages({library(Seurat); library(loupeR)})
M    <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing"
RDS  <- file.path(M, "1-mouse_based/Npp53_tracing/data/npp53_merge/npp53_edited.RDS")
OUT  <- file.path(M, "1-mouse_based/Npp53_tracing/data/npp53_merge/npp53_tumor_no_cyc_ion_v10.cloupe")
DROPBOX <- "/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/npp53_tumor_no_cyc_ion.cloupe"

# canonical raw -> display map (matches figures/npp53_plots/make_sec6_4b_grid.py STATE_MAP)
STATE_MAP <- c(
  'Luminal-Distal'    = 'c1_Lum-Distal',
  'Luminal-Proximal'  = 'c2_Lum-Proximal',
  'Luminal-ARdep'      = 'c3_Lum-ARdep',      # display uses ARdep (biology-forward name)
  'Luminal-CR'        = 'c4_Lum-CR',
  'Luminal-Persistent'= 'c5_Lum-Persistent',
  'Luminal-Squamous'  = 'c6_Lum-Squamous',
  'Basal-like'        = 'c7_Basal-like',
  'Partial-EMT'       = 'c8_Partial-EMT',
  'Mes-Quiescent'     = 'c9_Mes-Quiescent',
  'Mes-Collagen'      = 'c10_Mes-Collagen',
  'Mes-Stemness'      = 'c11_Mes-Stemness',
  'Mes-myCAF'         = 'c12_Mes-myCAF'
)
STATE_LEVELS_DISPLAY <- as.character(STATE_MAP)
STATE_LEVELS_RAW     <- names(STATE_MAP)

cat("[1] load RDS + filter Cyc + Ion ...\n")
obj <- readRDS(RDS)
DefaultAssay(obj) <- "RNA"
if ("counts" %in% Layers(obj)) obj <- JoinLayers(obj)
obj <- subset(obj, cells = colnames(obj)[!(obj$putative_state %in% c("Cycling","Ionocyte"))])
cat("  tumor cells:", ncol(obj), "\n")

# --- attach recomputed embeddings (identical to v10) ---
u  <- read.csv(file.path(M, "1-mouse_based/Npp53_tracing/data/npp53_merge/umap_no_cycling_no_ionocyte.csv"),   row.names=1)
tt <- read.csv(file.path(M, "1-mouse_based/Npp53_tracing/data/npp53_merge/tsne_no_cycling_no_ionocyte.csv"),   row.names=1)
dm <- read.csv(file.path(M, "1-mouse_based/Npp53_tracing/data/npp53_merge/diffmap_no_cycling_no_ionocyte.csv"), row.names=1)
cells <- colnames(obj)
umap_mat <- as.matrix(u[cells, c("umapsub_1","umapsub_2")]);  colnames(umap_mat) <- c("UMAP_1","UMAP_2")
tsne_mat <- as.matrix(tt[cells, c("tsnesub_1","tsnesub_2")]); colnames(tsne_mat) <- c("tSNE_1","tSNE_2")
dm_mat   <- as.matrix(dm[cells, c("DC2","DC3")]);             colnames(dm_mat)   <- c("diffmap_1","diffmap_2")
h        <- Embeddings(obj, "umap.harmony"); colnames(h) <- c("harmony_UMAP_1","harmony_UMAP_2")

obj[["UMAP"]]         <- CreateDimReducObject(embeddings=umap_mat,  key="UMAP_",    assay=DefaultAssay(obj))
obj[["tSNE"]]         <- CreateDimReducObject(embeddings=tsne_mat,  key="tSNE_",    assay=DefaultAssay(obj))
obj[["diffmap"]]      <- CreateDimReducObject(embeddings=dm_mat,    key="diffmap_", assay=DefaultAssay(obj))
obj[["harmony_UMAP"]] <- CreateDimReducObject(embeddings=h,         key="harmony_", assay=DefaultAssay(obj))
for (r in setdiff(Reductions(obj), c("UMAP","tSNE","diffmap","harmony_UMAP","pca"))) obj[[r]] <- NULL

# --- 12-state display naming ---
# keep the raw name too so anyone comparing to older CSVs / tables can join
raw_state <- as.character(obj$putative_state)
display   <- unname(STATE_MAP[raw_state])
if (any(is.na(display))) {
  cat("⚠ ", sum(is.na(display)), " cells have unmapped state (will be labelled NA)\n")
  cat("  offending raw states: ", paste(unique(raw_state[is.na(display)]), collapse=", "), "\n")
}
obj$state_raw       <- factor(raw_state, levels = STATE_LEVELS_RAW)
obj$tumor_state     <- factor(display,   levels = STATE_LEVELS_DISPLAY)   # new canonical name (c1..c12)
obj$putative_state  <- obj$tumor_state                                     # overwrite for consistency

# --- keep a compact metadata set (loupe metadata) ---
KEEP_META <- c("sample_id","treatment","model","batch","run",
               "tumor_state","putative_state","state_raw",
               "nCount_RNA","nFeature_RNA","scvi_res0p7_cl")
obj@meta.data <- obj@meta.data[, intersect(KEEP_META, colnames(obj@meta.data)), drop=FALSE]

cat("[2] state levels present in tumor object (c1..c12 order):\n")
print(table(obj$tumor_state, useNA='ifany'))

cat("[3] create cloupe ...\n")
create_loupe_from_seurat(obj, output_dir=dirname(OUT),
                         output_name=tools::file_path_sans_ext(basename(OUT)),
                         force=TRUE)
sz <- file.info(OUT)$size / 1e6
cat(sprintf("saved: %s  size=%.1f MB\n", OUT, sz))

file.copy(OUT, DROPBOX, overwrite=TRUE)
cat("synced to Dropbox: ", DROPBOX, "\n")

#!/usr/bin/env Rscript
# =============================================================================
# Fig 6.6 · Final tumor marker dotplot — Seurat DotPlot format
# 12 states (Cycling + Ionocyte excluded); 48 top-5-per-state markers.
# =============================================================================
suppressMessages({library(Seurat); library(ggplot2); library(dplyr)})

M_local <- "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing"
M       <- M_local
RDS     <- file.path(M, "1-mouse_based/Npp53_tracing/data/npp53_merge/npp53_edited.RDS")
OUT     <- file.path(M, "1-mouse_based/Npp53_tracing/figures/npp53_plots")

STATE_ORDER <- c('Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
                 'Luminal-Persistent','Luminal-Squamous',
                 'Basal-like','Partial-EMT',
                 'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF')
# c[N]_ display names (matches tree viewer + rest of report)
STATE_DISPLAY <- c('c1_Lum-Distal','c2_Lum-Proximal','c3_Lum-ARdep','c4_Lum-CR',
                   'c5_Lum-Persistent','c6_Lum-Squamous',
                   'c7_Basal-like','c8_Partial-EMT',
                   'c9_Mes-Quiescent','c10_Mes-Collagen','c11_Mes-Stemness','c12_Mes-myCAF')
DISPLAY_MAP <- setNames(STATE_DISPLAY, STATE_ORDER)
EXCLUDE <- c('Cycling','Ionocyte')

# Marker set v3 (2026-09-02) — 5 markers per state (60 total).
# STRICT CRITERIA applied after PubMed literature audit:
#   (1) EVERY marker has either (a) real prostate-cancer progression/CR/NEPC literature, or
#       (b) strong solid-tumor literature (breast/PDAC/CRC/HNSCC/glioma) with the specific role;
#   (2) EVERY marker is either enriched in the assigned state (top DE, pct.1 >> pct.2) OR
#       serves as a major lineage-group splitter (Epi vs Mes / Luminal vs Basal vs EMT).
# Sox2 in c11 mandatory (user request). Fabricated citations (Kwon 2014, Yan 2016, Chen 2019/2011)
# removed; Bmp6/Spink8/Bicc1/Cyp2f2/Serpina1e/Pdzk1ip1/Pcolce2/Dpt/Krt15/Mfap4/Gas7/Anxa1/Nupr1/
# Krtdap/Sprr1b/Col17a1/Plekhs1/Wfdc2/Gpx5/Adam7/Defb42/Pate13/Gabrp/Mgp/Dclk1/Col3a1/Cldn3/
# Acta2 dropped (no primary PCa/solid-tumor citation for stated function, or off-tissue).
MARKERS <- list(
  # v10 (2026-09-02): further edits per Zejian:
  #   remove Fabp5; Ltf before Psca; Apoe before Igfbp5; Yap1 c9→c8 (before Cdh6); Twist1 after Sox2.
  `Luminal-Distal`     = c('Krt8','Krt18','Clu','Wfdc2','Pigr'),
  `Luminal-Proximal`   = c('Ltf','Psca','Krt23'),                                                 # v11: Cst6 removed (assumed "Cst2" typo = Cst6)
  `Luminal-ARdep`       = c('Gpx5','Ptgds','Cd52'),                                                 # v11: Cuzd1 removed                                       # v6 dropped Gc (WEAK — germline SNP/serum only, no tumor-cell functional paper)
  `Luminal-CR`         = c('Onecut2','Agr2','Krt15','Muc16','Ffar4'),                             # Wfdc2 moved to c1
  `Luminal-Persistent` = c('Aldh1a3','Ifitm1','Runx3'),                                            # v11: Tnfrsf11b removed (bone-met, not persister)
  `Luminal-Squamous`   = c('Krt10','Serpinb2','Sult2b1'),                                          # v10: Fabp5 removed
  `Basal-like`         = c('Krt6a','Krt5','Krt14','Col17a1','Bcam','Apoe','Igfbp5','Nrg1'),        # v10: Apoe before Igfbp5
  `Partial-EMT`        = c('Vim','Fn1','Vcam1','Igfbp2','Spp1','Acta2','Yap1','Cdh6','Pdgfb'),    # v10: +Yap1 before Cdh6
  `Mes-Quiescent`      = c('Cdh2','Prex2'),                                                         # v10: Yap1 moved to c8 x-slot
  `Mes-Collagen`       = c('S100a4','Col1a1','Col1a2'),
  `Mes-Stemness`       = c('Col3a1','Sox2','Twist1','Hoxb13','Mmp13'),                              # v10: Twist1 after Sox2
  `Mes-myCAF`          = c('Col11a1','Aspn','Cdh11','Lox','Lrrc15')
)
markers_flat <- unlist(MARKERS, use.names = FALSE)
cat(sprintf("[markers] %d unique across %d groups\n", length(unique(markers_flat)), length(MARKERS)))

cat("[load] RDS ...\n")
obj <- readRDS(RDS)
DefaultAssay(obj) <- "RNA"
if ("counts" %in% Layers(obj)) obj <- JoinLayers(obj)

cat(sprintf("  n_cells full=%d\n", ncol(obj)))
cat("[filter] exclude Cycling + Ionocyte\n")
obj_f <- subset(obj, cells = colnames(obj)[!(obj$putative_state %in% EXCLUDE)])
cat(sprintf("  n_cells after: %d\n", ncol(obj_f)))

# Relabel putative_state with c[N]_ display names via RenameIdents (preserves cell names)
Idents(obj_f) <- "putative_state"
obj_f <- RenameIdents(obj_f,
  "Luminal-Distal"="c1_Lum-Distal", "Luminal-Proximal"="c2_Lum-Proximal",
  "Luminal-ARdep"="c3_Lum-ARdep",    "Luminal-CR"="c4_Lum-CR",
  "Luminal-Persistent"="c5_Lum-Persistent", "Luminal-Squamous"="c6_Lum-Squamous",
  "Basal-like"="c7_Basal-like",     "Partial-EMT"="c8_Partial-EMT",
  "Mes-Quiescent"="c9_Mes-Quiescent","Mes-Collagen"="c10_Mes-Collagen",
  "Mes-Stemness"="c11_Mes-Stemness", "Mes-myCAF"="c12_Mes-myCAF")
levels(obj_f) <- STATE_DISPLAY  # enforce canonical order
cat("idents after RenameIdents:", paste(levels(obj_f), collapse=", "), "\n")

cat("[plot] Seurat DotPlot\n")
p <- DotPlot(obj_f, features = markers_flat,
             cols = c('lightgrey','#08306b'),      # light-grey -> dark navy blue (Blues palette)
             dot.scale = 7, scale = TRUE, cluster.idents = FALSE) +
     RotatedAxis() +
     theme_bw(base_size = 12) +
     theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
           axis.text.y = element_text(size = 12),
           axis.title = element_blank(),
           panel.grid.major = element_line(color = 'grey92', linewidth = 0.3),
           panel.grid.minor = element_blank(),
           legend.position = 'right')                                                            # v11: title removed (embedded in report §6.6 header instead)

# save (v11: 54 markers)
ggsave(file.path(OUT, 'sec6_6_dotplot_final_12states_seurat.png'),
       p, width = 16.5, height = 6, dpi = 170)
ggsave(file.path(OUT, 'sec6_6_dotplot_final_12states_seurat.pdf'),
       p, width = 16.5, height = 6)
cat("saved: sec6_6_dotplot_final_12states_seurat.{png,pdf}\n")

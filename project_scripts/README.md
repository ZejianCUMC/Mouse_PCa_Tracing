# project_scripts: the scripts behind the Tracing reports (snapshot 2026-10-12)

This folder backs up the **current analysis scripts** of the mouse prostate-tumour lineage-tracing project
(Npp53, Npp53Rb1, and the merged Npp53 + Npp53Rb1 cohort). These scripts built the canonical objects and the
figures in the lab's internal reports. It sits next to the teaching workflow (`00_…06_` at the repository root),
which is separate.

**Code only.** No sequencing data, metadata tables, allele tables, character matrices, trees, cell-state
labels, `.h5ad/.RDS/.loom/.cloupe` files, figures, notebooks or report text are stored here. The repository
`.gitignore` whitelist enforces this. Report and HTML builders and three analysis scripts whose text holds
result numbers are kept on the lab drive only (see the end of this file).

How to read it:
- Top-level folders are analysis **stages**. Each subfolder names the model (`npp53`, `npp53rb1`,
  `merged` = Npp53 + Npp53Rb1) and the original task batch.
- Every file is a byte-identical copy of the working script. Original absolute paths (lab workstation and
  the c2b2 HPC) are left as they are, so a script runs only where those paths exist. The conda env is named
  in each script header or `.sbatch` file.
- Seeds are 42 throughout, with 19 and 888 for sensitivity runs.
- Report names: **Npp53 Part 1** (preprocessing and tumour states), **Npp53 Part 2** (evolution pattern
  profiling), **Npp53Rb1 Part 1** (preprocessing and tumour states), **merged report 1** (merged tumour
  states) and **merged report 2** (fitness risk modules, survival, cell differentiation maps).

---

## 01_preprocessing_demux

- `npp53_cellranger_cellbender_htodemux/`: Cell Ranger counts (gene expression + hashtag), CellBender
  ambient-RNA correction, HTODemux (`htodemux_drivers/`, singlets kept), DoubletFinder (`df_run.sh`) and the
  high-quality cell summary (`hq_summary.R`) for runs MJZ019 and MJZ021–024. Makes the per-sample singlet
  objects (Npp53 Part 1 §1–2).
- `npp53rb1_demuxmix_scdblfinder/`: CellBender hashtag input (`redemux_prep.R`, `extract_cbHTO.py`), demuxmix
  demultiplexing (`redemux_demuxmix.R`) with a convergence audit (`demuxmix_diag.R`) and the MJZ028 refit,
  scDblFinder (`scdbl_run.R`), then per-library singlet objects (`build_final_mtx.py`, `build_final_seurat.R`)
  for MJZ025–028. The `*_run.sh`, `submit_*.sh` and `build_*_env.sbatch` files are the HPC wrappers.

## 02_integration_cell_states

**Npp53**
- `npp53_merge_scvi_states/`: the Seurat/scVI chain behind the Npp53 tumour object.
  - `merge_npp53.R` → `scvi_npp53.py` → `rebuild_npp53.R` / `rebuild_scvi_npp53.py`: merge the tumours, then Harmony, scVI and fastMNN with UMAPs.
  - `annotate_npp53.R`, `annotate_detailed_npp53.R`, `relabel_trace_npp53.R`, `scvi_celltype_npp53.py`, `subcluster_mes_npp53.R`, `refine_tumor_npp53.R`: cell types, reporter-positive tumour cells, host-fibroblast removal.
  - `count_edited_cells.R`: cells with edited lineage barcodes.
  - `phaseA_edited_substates.R` → `substate_embeddings.R` → `relabel_scvi_native.R` → `set_scvi_primary.R` → `scvi_dense_umap.R` → `scvi_dotplots_markers_hpc.R` → `v3_updates.R`: build and update the edited-cell tumour object and its states (Part 1 Fig 4.1).
  - `scvi_resolution_bench.R`, `scvi_bootstrap_step01.R`, `scvi_bootstrap_stability.R`, `percluster_jaccard_res04_08.R`: resolution choice (benchmark and bootstrap stability).
  - `plots_*`, `umap_clusters_npp53.R`, `rerender_*.R`, `build_cloupe_npp53.R`: figures and the Loupe file.
- `npp53_part1_figures/`: UMAP / tSNE / Diffmap grids (`make_sec6_4*_grid*.py`, Part 1 §6.4) and the final
  marker dotplot (`sec6_6_dotplot_final*.{py,R}`, Part 1 §6.8).
- `npp53_major_cluster_stability/`: stability of the four major groups by five methods (hierarchical
  clustering, silhouette/gap, chooseR-style bootstrap, Jaccard, canonical markers) and Part 1 Fig 6.7.
- `npp53_ionocyte_per_arm/`: Ionocyte validation, cladograms and pruning from the trees; per-arm UMAP /
  Diffmap panels (Part 1 §6.6).
- `npp53_loupe_export/`: Loupe file of tumour cells without Cycling / Ionocyte, with the current state names.

**Npp53Rb1**
- `npp53rb1_integration/`
  - `rebuild_legacy_prepare.py` (with `feature_name_mapping.py`, `sample_retention.py`): merge and the per-sample cell filter.
  - `rebuild_legacy_scvi.py`, `rebuild_legacy_evaluate.py`: Harmony / fastMNN / scVI benchmark.
  - `review_min400.py`, `harmony_resolution_sweep.py`, `publish_harmony_major_types.R`, `plot_harmony_treatment.R`, `cluster_tumor_classification.py`: major cell types and the tumour-cell rule (Npp53Rb1 Part 1 §2).
  - `tumor_umap_*`: tumour-only re-integration, batch-correction comparison and the selected scVI latent.
- `npp53rb1_cell_states/`
  - `tumor_states_filtered.py`: cohort preparation and first embeddings (Fig 3.4).
  - `tumor_states_choir.R` with `prepare_choir_label_fix.py` and `choir_negative_level_fix.R`: CHOIR with the label-parser fix. `audit_choir_labels.R` and `compare_choir_trees.R` verify it.
  - `rank_corrected_choir.py`: AMI ranking of resolutions (Table 3.1).
  - `tumor_states_scshc.R`: run-aware sc-SHC.
  - `adopt_corrected_state_results.py` → `plot_top5_evidence.py`: Table 3.2 and Fig 3.5.
  - `profile_filtered_candidates.py`, `tumor_states_expression_review.py`, `review_filtered_scshc.py`, `cluster_specific_enrichment.py`: candidate profiling (Table 3.3S evidence).
  - `name_inferred_states.py`: state names and order.
  - `plot_filtered_candidate_review.py` → `plot_filtered_dotplot_seurat.R`: Fig 3.6 and 3.7.
  - `rebuild_without_c17_c19.py`: final cohort with a new graph and UMAP (Fig 3.8).
  - `final_state_projections.py`: tSNE and Diffmap.
  - `plot_final_state_splits.py`: Fig 3.9–3.11.
  - `export_final_loupe_inputs.py` → `build_final_state_cloupe.R`: final Loupe file.
  - Screening figures: `plot_tumor_states_reviewed.py` (Fig 3.2), `review_tumor_islands.py` and `plot_tumor_island_evidence.py` (Fig 3.3).

**Merged Npp53 + Npp53Rb1**
- `merged_integration/`: freeze cohorts and shared genes, choose HVGs, train the joint scVI (HPC) and run
  integration diagnostics (`prepare.py`, `hvg.py`, `train_scvi.py`, `evaluate.py`, `run_local.py`, `stage.sh`,
  `submit_initial.sh`). Resolution tests with `choir.R`, `scshc.R`, `after_choir.py` and `after_statistics.py`
  (merged report 1 Table 3.1).
- `merged_states/`
  - `restore_old_c18_umap.py`, `recluster_restored_c18.py` (Leiden, seeds 42 / 19 / 888), `merge_reviewed_states.py`, `prepare_reviewed10.py`, `profile_final_cleanup.py`, `name_reclustered.py` (with `plot_reviewed_umap.py`), `final_cleanup_dotplot.R`: the ten reviewed states.
  - `score_reviewed10_pathways.py`, `finalize_seed42_pathways.py`, `render_reviewed10_pathways.py`: AUCell pathway statistics (Fig 3.8).
  - `export_state_markers_excel.py`: marker workbook.
  - `loupe_*`, `prepare_reviewed10_loupe.py`, `build_merged_cloupe.R`: Loupe file.
- `merged_states_relabel_20261008/`: `relabel_and_plot.py` applies the approved state names and order (the
  canonical merged labels, numbered UMAP, dotplot). `render_fig36bc_fig38.py` draws Fig 3.6b/c and 3.8. The
  `*_loupe_relabel.py` and `export_markers_relabel.py` scripts rebuild the relabelled Loupe file and marker
  workbook.

## 03_lineage_trees

- `cassiopeia_preprocessing/`
  - `original_MJZ19_25_recipe/`: lineage (TSL) FASTQ → allele table; the recipe used for MJZ019–028 (copied from the HPC).
  - `standard/`: the current standard (`min_reads_per_umi = 1`, faster faithful intBC error correction).
  - `MJZ019_rerun_minreads1_20261004/`: the MJZ019 re-run and its QC against the original.
- `allele_tables/`: per-mouse allele tables on the final cohorts, plus the MJZ026–028 allele-table QC.
- `intbc_error_scan_20261004/`: finds near-duplicate (error) intBCs and writes tree inputs without them.
- `npp53_trees_v8/`
  - `rebuild_tree_one_sample.py`: trees with five solvers; Greedy is canonical.
  - `run_jungle_one_sample.py`: LBI fitness.
  - `run_phylovelo_one_sample.py`: PhyloVelo pseudotime, with per-sample parameters.
  - `rebuild_sc_fitness_plasticity_v10.py`: the per-cell LBI / plasticity / pseudotime table.
- `npp53_tree_patches/`: rebuild the tree, LBI and PhyloVelo for single samples, then patch the Npp53 viewer:
  - JZ136 without error intBCs;
  - JZ201–205 from the MJZ019 re-run;
  - JZ201 without one outlier cell.

  The Part 2 cohort analyses use the v8 trees.
- `npp53rb1_trees/`
  - `cassiopeia_eight_20260923`: Greedy trees for the first eight mice.
  - `tree_viewer_20260923`: JZ135 subsamples, alternative solvers, LBI, PhyloVelo and the first viewer builder.
  - `cassiopeia_mjz026_028_20260925` and `tree_viewer_mjz026_028_20260925`: the MJZ026–028 mice.
  - `intbcfix_20261004`: JZ137, JZ138 and the JZ135 subsamples rebuilt without error intBCs.
  - `JZ272_state_min20_20261006` and `state_min20_rebuild_20261006`: rebuilds keeping only states with at least 20 cells.

## 04_barcode_qc

- `npp53_editing_qc/`: Cas9 editing saturation and rate, allele QC, Manhattan plots and the combined
  Fig 6.1 (Npp53 Part 1 §6.1); clonal dominance; raw sequencing QC.
- `barcode_heatmaps/`: tree + intBC allele heatmap per mouse. Version 4 is current and imports v3, which
  imports v2. Versions 1 and 2 import helpers from `23-ancestorOT/report/gen_barcode_joint.py` (another lab
  project, not copied).
- `barcode_error_fixing_20261004/`: the tables and "before" panels of the barcode error-fixing report.

## 05_fitness_signatures

- `npp53_part2_fitness_signature/` (Npp53 Part 2 §9)
  - `derive_FitnessSignature.py`: per-tumour gene ~ LBI models.
  - `cache_cohorts.R`, `filter_cache_protein_coding.R`, `convert_orthologs.R`: human cohorts and mouse → human orthologs.
  - `plot_yang2022_*.py`: Fig 9.1.
  - `threshold_sweep/`: the meta-vote sweep and Hotspot; `T4_merged5_generate.py` (Fig 9.2–9.3); `T4_merged5_finalize.py` and `T4_merged5_xlsx.py` (Table 9.1 and the gene workbook); `T4_merged5_km/` (Kaplan–Meier, Fig 9.5–9.6).
  - `npp53_5FSU_module_activity_per_sample.py`: per-tumour module activity (§9.6).
- `merged_risk_modules/` (merged report 2 §4), in run order:
  1. `1_cohorts/`: human cohorts and orthologs.
  2. `2_fitness_model/`: Greedy-tree subclones and the per-tumour expression ~ Blom(LBI) model on log1p CP10k.
  3. `3_meta_lists/`: random-effects meta-analysis (DerSimonian–Laird, Hartung–Knapp, leave-one-tumour-out) and the FSU / FSD lists.
  4. `4_hotspot/`: Hotspot gene modules; the round driver is `run_meta_round.sh`.
  5. `5_candidates_survival/`: candidate gene sets, singscore splits, the random-set null, E2F adjustment and the castration effect.
  6. `6_selection/`: the fixed selection rule, chance level, and the pair statistics.
  7. `7_pair_figures/`: gene tables, heatmaps, volcanos and median-split KM curves.
  8. `8_report2_figures/01_fig43_fig42b.py`: Fig 4.3 and the Fig 4.2b legend.
  9. `9_final_modules/`: merge sub-modules by shared GO / Hallmark terms, then plot module scores per state and on the UMAP (Fig 4.5, Fig 4.6, the Table 4.3 modules).

## 06_survival_gsea

- `npp53_progression_gsea/`: FSU / FSD gene-set enrichment in human CRPC vs primary single-cell data
  (Npp53 Part 2).
- `merged_report2_survival/`
  - `01_km_all_cohorts.R`: five-cohort Kaplan–Meier with median and mean splits (§4.8 Fig 4.7–4.8, Tables 4.4–4.5, Appendix A).
  - `appC_pooled/`: pooled-tumour arm; candidate selection, cohorts, GO and rank-1 figures (Appendix C).
  - `appD_cycling_out/`: the same without cycling cells (Appendix D).
- `merged_report2_gsea_progression/`: progression GSEA with the current FSU / FSD lists (Appendix E).

## 07_evolution_features (Npp53 Part 2 §10)

- `npp53_PATH/`: PATH heritability and PATHpro γ on unit, edit-distance and pseudotime trees; tree
  annotation, per-arm violins, heritability vs plasticity and γ plots (§10.4).
- `npp53_cohort_heterogeneity/`: within- vs between-arm similarity (Fig 10.1); per-arm coverage and state
  universality (Fig 10.2).
- `npp53_tarca/`: TarCA founding-progenitor analysis.
- `npp53_evopath_trajectory/`: genes, pathways, regulons and FSU modules along pseudotime (Fig 10.3.10–10.3.14).

## 08_fate_maps_carta

- `npp53_carta_fine12/`: 12-state labels, pruned Greedy and NJ trees, CARTA (Gurobi) for k = 1–8,
  cross-solver stability and the chord-elbow k*.
- `npp53_carta_review/`: DAG extraction and checks, k* curves (Fig 10.3.1), per-sample DAGs and Greedy vs NJ
  (Fig 10.3.2–10.3.4), c7 / Fate-2 GSEA and dotplots (Fig 10.3.7, 10.3.9).
- `npp53_datadriven_fatemap/`: combined pseudotime, per-arm map (Fig 10.3.5) and cohort map (Fig 10.3.6).
  The pooling step that wrote the per-arm consensus file was never saved. Use
  `merged_report2_s5/step0_baseline/reconstruct_consensus_json.py` to reconstruct it.
- `merged_report2_s5/`: cell differentiation maps for merged report 2 §5 (Npp53, Npp53Rb1 and both; cycling
  cells kept or removed). Step 0 rebuilds the Part 2 Fig 10.3.6 baseline. Steps 1–5 are inputs, CARTA arrays,
  the tree position of cycling cells, pseudotime and the maps (Fig 5.x, Tables 5.1–5.2).

## 09_regulons_pathways

- `npp53_pyscenic/`: pySCENIC chain (loom, GRN, cisTarget, AUCell) and the top regulons per state (Npp53 Part 1 §6.9).
- `npp53_pathway_gmt/`: the mouse pathway GMT used by the AUCell and GSEA scripts.

## 10_viewers_videos

- `npp53_viewer/`: refreshes the CARTA panel of the Npp53 interactive tree viewer.
- `npp53rb1_viewer/`: builds the Npp53Rb1 interactive tree viewer.
- `npp53_videos/`: tumour-evolution animations.

---

## Not in this folder

- **Kept on the lab drive only:** report and HTML builders and patchers (they contain report text), and three
  analysis scripts whose text holds result numbers (the scripts that draw merged report 2 Fig 4.1, 4.2 and 4.4,
  and one Part 2 GSEA dotplot). The lab copy is `2-Tracing/github_script_backup_20261012/local_only/`.
- **Not found on disk:**
  - the step that wrote the Npp53 `putative_state` labels and exported the canonical Npp53 `.h5ad`;
  - Npp53 Part 2 Fig 9.1.2 and Table 9.2;
  - Npp53 Part 1 §6.10 (pyVIPER / ARACNe) and §7 (AR AUCell);
  - the TreeVAE runner;
  - the original Npp53 viewer builder.
- **Left out on purpose:** superseded rounds (for example the raw-count fitness rounds and earlier demultiplexing
  and integration rounds), tests, one-off HPC clean-up scripts, and the Jupyter notebooks of the pilot cohort.

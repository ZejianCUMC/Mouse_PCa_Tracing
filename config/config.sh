#!/bin/bash
# =============================================================================
# config.sh — ONE place for every path / sample setting used by the tutorial.
# Every step script does:   source config/config.sh
# Always submit jobs from the repository root (so relative paths resolve).
# =============================================================================

# ── 1. Your own output folder (edit this) ────────────────────────────────────
export WORK=/groups/ms3625_gp/${USER}/tracing_tutorial        # all outputs go here

# ── 2. Shared lab resources on c2b2 (read-only, do not edit) ─────────────────
export LAB=/groups/ms3625_gp/zw2994
export ENVS=${LAB}/mamba/envs                                  # conda envs
export CONDA_SH=${LAB}/mamba/etc/profile.d/conda.sh
export CELLRANGER=${LAB}/tools/cellranger-10.0.0/cellranger    # Cell Ranger 10.0.0
export GEX_REF=${LAB}/0-database/GRCm39_tracing                # GRCm39 + GFP/mCherry/Cas9/Cre transgenes
export TSL_REF=${LAB}/1-project/12-Jingbo/Nov13/refOct1.fa     # lineage-recorder target-site reference
export WHITELIST_10XV3=${LAB}/1-project/12-Jingbo/Nov13/3M-3pgex-may-2023.txt   # 10x 3' v3 barcode whitelist

# ── 3. Demo (subsampled) raw data — 2 Npp53 mice from 10x run MJZ019 ─────────
# Made by 00_make_demo_data/ (Zejian ran this once; students just read it).
export DEMO=${LAB}/1-project/12-Jingbo/2-Tracing/teaching_demo_MJZ019
export GEX_FASTQ_DIR=${DEMO}/fastq/GEX     # MJZ019sub_S1_L001_R{1,2}_001.fastq.gz   (10x 3' GEX)
export HTO_FASTQ_DIR=${DEMO}/fastq/HTO     # MJZ019Fsub_S1_L001_R{1,2}_001.fastq.gz  (TotalSeq-B hashtags)
export TSL_FASTQ_DIR=${DEMO}/fastq/TSL     # MJZ019TSLsub_R{1,2}.fastq.gz            (lineage amplicon, AVITI)
export GEX_SAMPLE=MJZ019sub                # fastq prefix of the GEX library
export HTO_SAMPLE=MJZ019Fsub               # fastq prefix of the hashtag library
export TSL_LIB=MJZ019TSLsub                # fastq prefix of the lineage library
export RUN=MJZ019sub                       # run / Cell Ranger id used in all output names

# Samples (mice) to carry downstream; the others in the pool are dropped after demux
export KEEP_SAMPLES="JZ202,JZ204"

# ── 4. Output sub-folders (derived; no need to edit) ─────────────────────────
export OUT_CR=${WORK}/01_cellranger
export CR_OUTS=${OUT_CR}/${RUN}/outs          # Cell Ranger outs (filtered_feature_bc_matrix/ lives here)
export OUT_CB=${WORK}/02_cellbender
export OUT_QC=${WORK}/03_demux_qc
export OUT_ST=${WORK}/04_cell_states
export OUT_CAS=${WORK}/05_cassiopeia
export OUT_TREE=${WORK}/06_greedy_tree

# ── 5. Reproducibility ───────────────────────────────────────────────────────
export SEED=42                             # lab rule: 42 primary; 19 / 888 for sensitivity

# ── 6. OPTION B: start at step 03 with the provided tutorial data folder ─────
# Zejian hands over "Rotation_tutorial_data/" (on a hard drive). Per run it holds:
#   Rotation_tutorial_data/<RUN>/01_scRNA/    Cell Ranger outs (filtered_feature_bc_matrix/, web_summary.html)
#                                             + CellBender <RUN>_cellbender_filtered.h5
#   Rotation_tutorial_data/<RUN>/02_barcode/  <RUN>_allele_table.csv  (Cassiopeia, one table per run)
# To use it: set TUTORIAL_DATA and RUN below, uncomment, then run steps 03 -> 04 -> 06
# (on HPC with the sbatch files, or on your own machine with run_local_from_step03.sh).
#
# export TUTORIAL_DATA=/path/to/Rotation_tutorial_data
# export RUN=MJZ019                                   # MJZ019 (5 Npp53 mice) or MJZ008 (3 pilot mice)
# export TSL_LIB=${RUN}
# export KEEP_SAMPLES="JZ201,JZ202,JZ203,JZ204,JZ205" # MJZ008: "JZ136,JZ137,JZ138"
# export WORK=${HOME}/tracing_tutorial/${RUN}         # your outputs (03, 04, 06)
# export CR_OUTS=${TUTORIAL_DATA}/${RUN}/01_scRNA     # Cell Ranger outs          -> step 03
# export OUT_CB=${TUTORIAL_DATA}/${RUN}/01_scRNA      # CellBender filtered .h5   -> step 03
# export OUT_CAS=${TUTORIAL_DATA}/${RUN}/02_barcode   # allele table              -> step 06
# export OUT_QC=${WORK}/03_demux_qc
# export OUT_ST=${WORK}/04_cell_states
# export OUT_TREE=${WORK}/06_greedy_tree

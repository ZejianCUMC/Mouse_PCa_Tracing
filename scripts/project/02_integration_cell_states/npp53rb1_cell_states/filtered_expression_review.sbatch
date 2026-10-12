#!/bin/bash
#SBATCH --job-name=rb1-filtered-expression
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/expression_%j.log
set -euo pipefail
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing
date --iso-8601=seconds
hostname
sha256sum scripts/tumor_states_expression_review_filtered.py \
 analysis/min400_20260921/tumor_umap_20260921/tumor_full_counts.h5ad \
 analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/metadata_PRIVATE.csv \
 analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/resolution_assignments_PRIVATE.csv \
 analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/Npp53_reference_top25.csv
/groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python -u scripts/tumor_states_expression_review_filtered.py \
 --source analysis/min400_20260921/tumor_umap_20260921 \
 --output analysis/tumor_states_20260922/exclude_cl5_cl8_seed42 \
 --reference analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/Npp53_reference_top25.csv \
 --seed 42 --skip-islands
date --iso-8601=seconds

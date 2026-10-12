#!/bin/bash
set -euo pipefail
export RB1_RUN=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/analysis/min400_20260921
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8
scripts="$RB1_RUN/tumor_umap_20260921/scripts"
case "$1" in
 prepare) /groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python "$scripts/tumor_umap_prepare.py";;
 harmony) /groups/ms3625_gp/zw2994/mamba/envs/ST/bin/Rscript "$scripts/tumor_umap_harmony.R";;
 scvi)
  module load cuda
  modes=(run sample run_sample)
  /groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python "$scripts/tumor_umap_scvi.py" "${modes[$SLURM_ARRAY_TASK_ID]}";;
 evaluate) /groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python "$scripts/tumor_umap_evaluate.py";;
 *) exit 2;;
esac

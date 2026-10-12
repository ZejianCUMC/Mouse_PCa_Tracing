#!/bin/bash
set -euo pipefail
export RB1_BASE=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing
export RB1_RUN="$RB1_BASE/analysis/min400_20260921"
export RB1_MIN_VALID_CELLS=400
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8
scripts="$RB1_RUN/scripts"
case "$1" in
 prepare_integrate)
  /groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python "$scripts/rebuild_legacy_prepare.py"
  /groups/ms3625_gp/zw2994/mamba/envs/ST/bin/Rscript "$scripts/rebuild_legacy_integrate.R"
  ;;
 scvi)
  module load cuda
  /groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python "$scripts/rebuild_legacy_scvi.py"
  ;;
 evaluate)
  /groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python "$scripts/rebuild_legacy_evaluate.py"
  # Manual review is mandatory: old cluster-number labels must NOT be reused.
  ;;
 *) exit 2;;
esac

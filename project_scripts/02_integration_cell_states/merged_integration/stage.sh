#!/bin/bash
set -euo pipefail
project=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp_tumor_merged
out="$project/analysis/merged_20260923_seed42"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export NUMBA_CACHE_DIR="$out/cache/numba" MPLCONFIGDIR="$out/cache/mpl"
mkdir -p "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR"
python=/groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python
case "$1" in
 hvg) "$python" -u "$project/scripts/hvg.py" "$out";;
 train) modes=(run run_sample sample); "$python" -u "$project/scripts/train_scvi.py" "$out" "${modes[$SLURM_ARRAY_TASK_ID]}";;
 evaluate) "$python" -u "$project/scripts/evaluate.py" "$out";;
 hvg_endogenous) "$python" -u "$project/scripts/hvg.py" "$out" endogenous;;
 train_endogenous) "$python" -u "$project/scripts/train_scvi.py" "$out" "$2" 42 endogenous;;
 evaluate_endogenous) "$python" -u "$project/scripts/evaluate.py" "$out" endogenous "$2";;
 states) "$python" -u "$project/scripts/prepare_states.py" "$out" "$2";;
 seed_sensitivity) seeds=(19 888); "$python" -u "$project/scripts/train_scvi.py" "$out" run_sample "${seeds[$SLURM_ARRAY_TASK_ID]}" endogenous;;
 seed_graphs) seeds=(19 888); "$python" -u "$project/scripts/seed_graphs.py" "$out" "${seeds[$SLURM_ARRAY_TASK_ID]}";;
 *) exit 2;;
esac

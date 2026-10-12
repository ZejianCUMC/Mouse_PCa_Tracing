#!/bin/bash
#SBATCH --job-name=Rb1_Harmony_res
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
export RB1_RUN=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/analysis/min400_20260921
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8
/groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python scripts/harmony_resolution_sweep.py

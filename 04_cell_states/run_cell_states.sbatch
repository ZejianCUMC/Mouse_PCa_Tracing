#!/bin/bash
# =============================================================================
# 04 · Tumour cells + tumour cell states  (python, env scvi)
#   04a: major types, reporter (GFP & mCherry) -> tumour cells
#   04b: tumour-only scVI, Leiden sweep + bootstrap stability, markers -> cell_states.csv
# Submit from repo root:  sbatch 04_cell_states/run_cell_states.sbatch
# Re-run 04b only with a chosen resolution, e.g.:
#   sbatch --export=ALL,RES=0.7 04_cell_states/run_cell_states.sbatch
# Demo size (~1-2k cells) trains fine on CPU; full cohorts go to the gpu partition.
# =============================================================================
#SBATCH --job-name=cell_states
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
source config/config.sh
export PATH=${ENVS}/scvi/bin:$PATH
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4 MKL_NUM_THREADS=4

if [ -z "${RES:-}" ]; then
    python -u 04_cell_states/04a_tumor_cells.py
    python -u 04_cell_states/04b_tumor_states.py
else
    python -u 04_cell_states/04b_tumor_states.py --res "${RES}"
fi

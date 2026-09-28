#!/bin/bash
# =============================================================================
# 02 · CellBender remove-background: remove ambient RNA (and ambient hashtag)
# In : ${OUT_CR}/${RUN}/outs/raw_feature_bc_matrix.h5  (RAW, not filtered!)
# Out: ${OUT_CB}/${RUN}_cellbender.h5, ..._filtered.h5 (downstream count matrix), .pdf (QC)
# Submit from repo root:  sbatch 02_cellbender/run_cellbender.sbatch   (GPU job)
# =============================================================================
#SBATCH --job-name=cellbender
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-gpu=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -eo pipefail                      # not -u: conda activation scripts use unset vars
cd "${SLURM_SUBMIT_DIR}"
source config/config.sh
source ${CONDA_SH}
export PATH=${ENVS}/cellbender/bin:$PATH      # never `conda activate` inside batch jobs
module load cuda
mkdir -p ${OUT_CB} && cd ${OUT_CB}

cellbender remove-background --cuda \
    --input  ${OUT_CR}/${RUN}/outs/raw_feature_bc_matrix.h5 \
    --output ${OUT_CB}/${RUN}_cellbender.h5
# Look at ${RUN}_cellbender.pdf: the learning curve should converge and the
# cell-probability plot should separate cells from empty droplets.

#!/bin/bash
# =============================================================================
# 03 · Demultiplex hashtags + QC + DoubletFinder  (R / Seurat, env ST)
# Submit from repo root:  sbatch 03_demux_qc/run_demux_qc.sbatch
# =============================================================================
#SBATCH --job-name=demux_qc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
source config/config.sh
bash config/hashtags.sh                      # (re)writes ${WORK}/metadata/hash_map.csv
export PATH=${ENVS}/ST/bin:$PATH
ulimit -s unlimited                          # DoubletFinder paramSweep can overflow the stack
Rscript --max-ppsize=500000 03_demux_qc/demux_qc_doublets.R

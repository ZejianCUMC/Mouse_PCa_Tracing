#!/bin/bash
# =============================================================================
# 05 · Cassiopeia preprocessing of the lineage (TSL) library -> allele table
# Independent of steps 02-04: can be submitted right after the demo data exist.
# Submit from repo root:  sbatch 05_cassiopeia/run_cassiopeia.sbatch
# Full libraries: 16 CPU and --mem=384G (error_correct_umis thrashes at 96G on big libs).
# =============================================================================
#SBATCH --job-name=cassiopeia_pp
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
source config/config.sh
export PATH=${ENVS}/cassiopeia_env/bin:$PATH
python -u 05_cassiopeia/run_cassiopeia_pp.py ${TSL_LIB} \
    ${TSL_FASTQ_DIR}/${TSL_LIB}_R1.fastq.gz ${TSL_FASTQ_DIR}/${TSL_LIB}_R2.fastq.gz \
    ${OUT_CAS} ${SLURM_CPUS_PER_TASK} ${TSL_REF} ${WHITELIST_10XV3}

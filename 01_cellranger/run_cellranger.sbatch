#!/bin/bash
# =============================================================================
# 01 · Cell Ranger count: GEX + hashtag (Antibody Capture) -> count matrices
# In : ${GEX_FASTQ_DIR}, ${HTO_FASTQ_DIR}   Out: ${OUT_CR}/${RUN}/outs/
# Key outputs: raw_feature_bc_matrix.h5 (-> CellBender), filtered_feature_bc_matrix/
#              (-> demux), web_summary.html (check this first!)
# Submit from repo root:  sbatch 01_cellranger/run_cellranger.sbatch
# Full-size real runs used 12 CPU / 96 GB / --create-bam true; the demo needs far less.
# =============================================================================
#SBATCH --job-name=cellranger
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source config/config.sh
bash config/hashtags.sh                      # writes ${WORK}/metadata/feature_ref.csv
mkdir -p ${OUT_CR}

# library CSV: which fastq folder is which library type
cat > ${WORK}/metadata/library.csv <<CSV
fastqs,sample,library_type
${GEX_FASTQ_DIR},${GEX_SAMPLE},Gene Expression
${HTO_FASTQ_DIR},${HTO_SAMPLE},Antibody Capture
CSV

cd ${OUT_CR}
${CELLRANGER} count --id=${RUN} \
    --libraries=${WORK}/metadata/library.csv \
    --feature-ref=${WORK}/metadata/feature_ref.csv \
    --transcriptome=${GEX_REF} \
    --localcores=8 --localmem=60 \
    --create-bam false
# NB: --feature-ref must list EVERY hashtag used in the pool (here 5; newer runs use all 6, B301-B306).

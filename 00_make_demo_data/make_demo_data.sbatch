#!/bin/bash
# =============================================================================
# 00 · Build the tutorial demo data (run ONCE by Zejian; students skip this)
# Keeps reads of ~800 cells x 2 mice (JZ202 Cas+Enz, JZ204 Castration) + 300
# doublet/negative barcodes + 15k empty droplets from 10x run MJZ019, in all
# three libraries (GEX, HTO, TSL). Output -> ${DEMO}
# Submit from repo root:  sbatch 00_make_demo_data/make_demo_data.sbatch
# =============================================================================
#SBATCH --job-name=make_demo_data
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source config/config.sh

SRC=${LAB}/1-project/12-Jingbo/2-Tracing
GEX=${SRC}/0-raw/260624_MICHAEL_JINGBO_1_MOUSE_10X/MJZ019/fastq/SC_3PGEX_10X/260813_LH00364_0685_A252GYTLT4
HTO=${SRC}/0-raw/260624_MICHAEL_JINGBO_1_MOUSE_10X/MJZ019F/fastq/FB_TSB_10X/260813_LH00364_0685_A252GYTLT4
TSL=${SRC}/0-raw/MJZ019_21_22_bases2fastq/Samples/DefaultProject/MJZ019
mkdir -p ${DEMO}

cat > ${DEMO}/sources.json <<EOF
{
  "hash_meta": "${SRC}/seurat_objs/MJZ019_cell_hash_meta.csv",
  "cr_outs":   "${SRC}/1-counts/MJZ019/outs",
  "fastqs": {
    "GEX": ["${GEX}/MJZ019_S2_L001_R1_001.fastq.gz",  "${GEX}/MJZ019_S2_L001_R2_001.fastq.gz",
            "${GEX_SAMPLE}_S1_L001_R1_001.fastq.gz",   "${GEX_SAMPLE}_S1_L001_R2_001.fastq.gz"],
    "HTO": ["${HTO}/MJZ019F_S3_L001_R1_001.fastq.gz", "${HTO}/MJZ019F_S3_L001_R2_001.fastq.gz",
            "${HTO_SAMPLE}_S1_L001_R1_001.fastq.gz",   "${HTO_SAMPLE}_S1_L001_R2_001.fastq.gz"],
    "TSL": ["${TSL}/MJZ019_R1.fastq.gz",               "${TSL}/MJZ019_R2.fastq.gz",
            "${TSL_LIB}_R1.fastq.gz",                  "${TSL_LIB}_R2.fastq.gz"]
  }
}
EOF

export OMP_NUM_THREADS=1
${ENVS}/scvi/bin/python -u 00_make_demo_data/subsample_fastq_by_barcode.py \
    --sources ${DEMO}/sources.json --outdir ${DEMO}
ls -la ${DEMO}/fastq/*

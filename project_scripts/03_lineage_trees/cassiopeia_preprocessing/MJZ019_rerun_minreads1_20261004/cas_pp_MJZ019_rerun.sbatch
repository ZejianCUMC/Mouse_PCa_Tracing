#!/bin/bash
#SBATCH --job-name=cas_pp_MJZ019_rerun
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=96:00:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_MJZ019_rerun_minreads1_20261004/logs/cas_pp_%j.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_MJZ019_rerun_minreads1_20261004/logs/cas_pp_%j.err
# MJZ019 TSL re-run (Zejian 2026-10-04): original recipe, min_reads_per_umi fixed to 1. New folder; the existing
# cassiopeia_MJZ19_22/MJZ019 results are not touched. Uninterrupted run (no checkpoint resume).
set -eo pipefail
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/cassiopeia_env/bin:$PATH
BASE=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_MJZ019_rerun_minreads1_20261004
FQ=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/0-raw/MJZ019_21_22_bases2fastq/Samples/DefaultProject
LIB=MJZ019
test ! -e "$BASE/$LIB/${LIB}_allele_table.csv"
echo "host=$(hostname) lib=$LIB cpus=$SLURM_CPUS_PER_TASK start=$(date)"
python "$BASE/scripts/run_cassiopeia_pp_minreads1.py" "$LIB" "$FQ/$LIB/${LIB}_R1.fastq.gz" "$FQ/$LIB/${LIB}_R2.fastq.gz" "$BASE/$LIB" "$SLURM_CPUS_PER_TASK"
echo "end=$(date)"
if [ -s "$BASE/$LIB/${LIB}_allele_table.csv" ]; then echo "== SUCCESS: $(ls -la --block-size=M $BASE/$LIB/${LIB}_allele_table.csv)"; else echo "== NO ALLELE TABLE"; exit 3; fi

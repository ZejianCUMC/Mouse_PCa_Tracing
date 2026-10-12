#!/bin/bash
#SBATCH --job-name=cas_pp_MJZ19_22
#SBATCH --array=0-2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=36:00:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_MJZ19_22/logs/cas_pp_%A_%a.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_MJZ19_22/logs/cas_pp_%A_%a.err
set -eo pipefail
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/cassiopeia_env/bin:$PATH
BASE=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_MJZ19_22
FQ=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/0-raw/MJZ019_21_22_bases2fastq/Samples/DefaultProject
LIBS=(MJZ019 MJZ021 MJZ022)
LIB=${LIBS[$SLURM_ARRAY_TASK_ID]}
echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID lib=$LIB cpus=$SLURM_CPUS_PER_TASK start=$(date)"
python "$BASE/scripts/run_cassiopeia_pp.py" "$LIB" "$FQ/$LIB/${LIB}_R1.fastq.gz" "$FQ/$LIB/${LIB}_R2.fastq.gz" "$BASE/$LIB" "$SLURM_CPUS_PER_TASK"
echo "end=$(date)"

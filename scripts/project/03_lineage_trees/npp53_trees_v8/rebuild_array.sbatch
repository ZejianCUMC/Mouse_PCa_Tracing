#!/bin/bash
#SBATCH --job-name=cas_rebuild_v8
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/logs/rebuild_%x_%A_%a.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/logs/rebuild_%x_%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-24%8
set -euo pipefail
mkdir -p /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/logs

SAMPLES=(JZ136 JZ201 JZ202 JZ203 JZ204 JZ205 JZ206 JZ207 JZ208 JZ209 JZ210 JZ211 JZ212 JZ213 JZ214 JZ224 JZ225 JZ226 JZ227 JZ228 JZ229 JZ230 JZ231 JZ232 JZ234)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
SOLVER=${SOLVER:-greedy}
OUTDIR=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_trees_v8_rebuilt/${SOLVER}

echo "== host=$(hostname) sample=${SAMPLE} solver=${SOLVER} start=$(date) =="
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/cassiopeia_env/bin:$PATH

python /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/rebuild_tree_one_sample.py \
    ${SAMPLE} ${SOLVER} ${OUTDIR}

echo "== end=$(date) =="

#!/bin/bash
#SBATCH --job-name=phylovelo_retry
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/logs/phylo_retry_%x_%A_%a.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/logs/phylo_retry_%x_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-1%2
set -euo pipefail
SAMPLES=(JZ136 JZ203)   # the 2 failed with 3D diffmap
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "== host=$(hostname) sample=${SAMPLE} tool=phylovelo-2D start=$(date) =="
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/PhyloVelo/bin:$PATH
python /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/run_phylovelo_one_sample.py ${SAMPLE}
echo "== end=$(date) =="

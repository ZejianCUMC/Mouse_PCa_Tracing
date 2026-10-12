#!/bin/bash
#SBATCH --job-name=phylovelo_JZ213
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/logs/JZ213_rerun_%A.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/logs/JZ213_rerun_%A.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
set -euo pipefail
echo "== host=$(hostname) sample=JZ213 tool=phylovelo-NNfix start=$(date) =="
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/PhyloVelo/bin:$PATH
# Remove stale flag files first so a successful run is unambiguous
rm -f /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_trees_v8_rebuilt/phylovelo_out/JZ213_*_RUNTIME_FAIL.flag
python /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_pathway_analysis/hpc_tree_rebuild/run_phylovelo_one_sample.py JZ213
echo "== end=$(date) =="

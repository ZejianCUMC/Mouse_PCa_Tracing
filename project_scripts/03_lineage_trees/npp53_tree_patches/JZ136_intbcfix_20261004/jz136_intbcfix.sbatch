#!/bin/bash
# JZ136 error-intBC rebuild (Zejian 2026-10-04). Same scripts/params as v8 (rebuild_array.sbatch +
# downstream_array.sbatch); only paths differ. Task = one step:
#   STEP=tree  SOLVER=greedy|nj|maxcut|spectral|percolation   (cassiopeia_env)
#   STEP=jungle                                             (cassiopeia_env)
#   STEP=phylovelo                                          (PhyloVelo env)
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
set -euo pipefail
RUN=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/lineage_tracing/JZ136_intbcfix_20261004
M=/groups/ms3625_gp/zw2994/mamba/envs
echo "== host=$(hostname) step=${STEP} solver=${SOLVER:-} start=$(date) =="
case "$STEP" in
  tree)      export PATH=$M/cassiopeia_env/bin:$PATH; python $RUN/scripts/rebuild_tree_one_sample.py JZ136 ${SOLVER} $RUN/${SOLVER} ;;
  jungle)    export PATH=$M/cassiopeia_env/bin:$PATH; python $RUN/scripts/run_jungle_one_sample.py JZ136 ;;
  phylovelo) export PATH=$M/PhyloVelo/bin:$PATH;      python $RUN/scripts/run_phylovelo_one_sample.py JZ136 ;;
  *) echo "unknown STEP=$STEP"; exit 1 ;;
esac
echo "== end=$(date) =="

#!/bin/bash
# JZ201-205 from the MJZ019 re-run (Zejian 2026-10-05). Same scripts/params as v8 (rebuild_array.sbatch + downstream_array.sbatch);
# only paths differ. One task = one sample x one step:  STEP=tree SOLVER=greedy|nj|maxcut|spectral|percolation | STEP=jungle | STEP=phylovelo
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=6:00:00
set -euo pipefail
RUN=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/lineage_tracing/JZ201_excl_outlier_20261006
M=/groups/ms3625_gp/zw2994/mamba/envs
echo "== host=$(hostname) sample=${SAMPLE} step=${STEP} solver=${SOLVER:-} start=$(date) =="
case "$STEP" in
  tree)      export PATH=$M/cassiopeia_env/bin:$PATH; python $RUN/scripts/rebuild_tree_one_sample.py ${SAMPLE} ${SOLVER} $RUN/${SOLVER} ;;
  jungle)    export PATH=$M/cassiopeia_env/bin:$PATH; python $RUN/scripts/run_jungle_one_sample.py ${SAMPLE} ;;
  phylovelo) export PATH=$M/PhyloVelo/bin:$PATH;      python $RUN/scripts/run_phylovelo_one_sample.py ${SAMPLE} ;;
  *) echo "unknown STEP=$STEP"; exit 1 ;;
esac
echo "== end=$(date) =="

#!/bin/bash
# Per-sample state >= 20 cells rebuild (Zejian 2026-10-06; same rule as the adopted JZ272 test). Same scripts/params as the original
# Rb1 MJZ026-028 runs; only change = casdir path in viewer/scripts/run_tree_viewer_analysis.py.  Env: SAMPLE, STEP.
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
set -euo pipefail
R=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/lineage_tracing/state_min20_rebuild_20261006/${SAMPLE}
E=/groups/ms3625_gp/zw2994/mamba/envs
export PYTHONHASHSEED=42 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MPLCONFIGDIR=$R/mpl_${STEP}
echo "== sample=$SAMPLE step=$STEP host=$(hostname) start=$(date -Is) =="
case "$STEP" in
  Greedy)    $E/cassiopeia_env/bin/python $R/cassiopeia/scripts/run_cassiopeia_eight.py --run-dir $R/cassiopeia --sample $SAMPLE ;;
  PhyloVelo) $E/PhyloVelo/bin/python      $R/viewer/scripts/run_tree_viewer_analysis.py --root $R/viewer --sample $SAMPLE --method PhyloVelo ;;
  LBI|NJ|Spectral|MaxCut|Percolation) $E/cassiopeia_env/bin/python $R/viewer/scripts/run_tree_viewer_analysis.py --root $R/viewer --sample $SAMPLE --method $STEP ;;
  *) echo "unknown STEP=$STEP"; exit 1 ;;
esac
echo "== end=$(date -Is) =="

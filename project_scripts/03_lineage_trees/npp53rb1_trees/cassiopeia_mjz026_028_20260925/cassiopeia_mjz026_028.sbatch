#!/bin/bash
#SBATCH --job-name=rb1_cas14_greedy
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=0-13%4
set -euo pipefail
RUN_DIR=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/lineage_tracing/cassiopeia_mjz026_028_20260925_seed42
CAS_PY=/groups/ms3625_gp/zw2994/mamba/envs/cassiopeia_env/bin/python
SAMPLES=(JZ265 JZ266 JZ267 JZ268 JZ269 JZ270 JZ271 JZ272 JZ273 JZ274 JZ276 JZ277 JZ279 JZ280)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
export PYTHONHASHSEED=42 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR="$RUN_DIR/logs/mpl_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$MPLCONFIGDIR"
echo "START sample=$SAMPLE job=$SLURM_JOB_ID host=$(hostname) at=$(date -Is)"
"$CAS_PY" "$RUN_DIR/scripts/run_cassiopeia_eight.py" --run-dir "$RUN_DIR" --sample "$SAMPLE"
echo "END sample=$SAMPLE at=$(date -Is)"

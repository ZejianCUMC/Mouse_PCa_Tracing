#!/bin/bash
# opt5 step 05 on c2b2 HPC — one SLURM array task per list; hotspot_tasks.txt line = "<arm> <list_id> <k>" (both arms in one queue).
# Same 05_hotspot_run.py as locally, HPC hotspot env (versions identical to local). Sized from the opt4 run (job 846378: <= 45 s,
# <= 1 GB per task at ~20k cells): 5 CPUs (Hotspot jobs=5), 4 GB, 30 min. Submit: sbatch --array=1-<n>%80 05h_hotspot_array.sbatch
#SBATCH --job-name=opt5_hotspot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=5
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/hs_%A_%a.out
set -euo pipefail
cd /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp_tumor_merged/opt5_corrected_20261009
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/hotspot/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=42
TASKS=${TASKS:-hotspot_tasks.txt}   # task file (env; default = main run)
read -r ARM LID K TOPN <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKS")"; TOPN=${TOPN:-500}; export TOPN; SUF=""; [ "$TOPN" != "500" ] && SUF="_t$TOPN"
echo "task ${SLURM_ARRAY_TASK_ID}: ${ARM} ${LID} k=${K} top=${TOPN} on $(hostname) at $(date)"
if [ -f "${ARM}/hotspot_out/${LID}_k${K}${SUF}.json" ]; then echo "already done"; exit 0; fi
python 05_hotspot_run.py "$ARM" "$LID" "$K"
echo "finished $(date)"

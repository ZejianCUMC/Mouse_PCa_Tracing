#!/bin/bash
# =============================================================================
# 06 · One VanillaGreedy lineage tree per mouse (array job, one task per mouse)
# Needs step 04 (cell_states.csv) and step 05 (allele table).
# Submit from repo root:  sbatch 06_greedy_tree/run_greedy_tree.sbatch
# Other solvers (robustness):  sbatch --export=ALL,SOLVER=nj 06_greedy_tree/run_greedy_tree.sbatch
# =============================================================================
#SBATCH --job-name=greedy_tree
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
source config/config.sh
export PATH=${ENVS}/cassiopeia_env/bin:$PATH
IFS=',' read -ra SAMPLES <<< "${KEEP_SAMPLES}"      # array size must equal the number of mice
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
python -u 06_greedy_tree/build_greedy_tree.py --sample ${SAMPLE} --solver ${SOLVER:-greedy}

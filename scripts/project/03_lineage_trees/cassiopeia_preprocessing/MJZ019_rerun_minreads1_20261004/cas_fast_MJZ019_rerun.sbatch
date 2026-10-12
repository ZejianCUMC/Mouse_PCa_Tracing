#!/bin/bash
#SBATCH --job-name=cas_fast_MJZ019_rerun
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=06:00:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_MJZ019_rerun_minreads1_20261004/logs/cas_fast_%j.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_MJZ019_rerun_minreads1_20261004/logs/cas_fast_%j.err
# Fast branch (Zejian 2026-10-04): starts from the main job's error_correct_umis checkpoint (read-only) and finishes
# steps 9-10 with the validated faithful intBC error correction. Writes ONLY to fast_faithful/MJZ019/.
set -eo pipefail
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/cassiopeia_env/bin:$PATH
BASE=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cassiopeia_MJZ019_rerun_minreads1_20261004
SRC=$BASE/MJZ019/MJZ019_error_correct_umis.csv
[ -s "$SRC" ] || { echo "no checkpoint"; exit 2; }
echo "host=$(hostname) start=$(date) src_size=$(stat -c %s $SRC)"
python "$BASE/scripts/fast_finish_from_ecumis.py" MJZ019 "$SRC" "$BASE/fast_faithful/MJZ019" "$SLURM_CPUS_PER_TASK"
echo "end=$(date)"
[ -s "$BASE/fast_faithful/MJZ019/MJZ019_allele_table.csv" ] && echo "== SUCCESS" || { echo "== NO ALLELE TABLE"; exit 3; }

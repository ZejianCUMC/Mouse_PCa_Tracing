#!/bin/bash
#SBATCH --job-name=carta_fatemap_c10
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/carta_out/fate_map_cycling_20261011/logs/%x_%A_%a.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/carta_out/fate_map_cycling_20261011/logs/%x_%A_%a.err
#SBATCH --array=0-55%12
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=cpu

# CARTA for the Fig 10.3.6-style cohort maps (fate_map_cycling_20261011, local step2_carta/)
# 47 tumors (Npp53 25 + Npp53Rb1 22), merged 10-state labels, current Greedy trees (report 2 Appendix D manifest).
# Tasks (inputs/carta_tasks.tsv): 47 "inc" (c10 kept) + 9 "exc" (c10 removed; only tumors with >=20 c10 cells,
# elsewhere the exc inputs are byte-identical to inc and reuse the inc result).
# CARTA settings as Npp53 carta_greedy_fine12.sbatch; k runs 1..min(8, n_states) because the chord-elbow
# selection never uses k > n_states (plateau), so the extra k values of the original only cost CPU.
set -euo pipefail
BASE=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/carta_out/fate_map_cycling_20261011
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/carta/bin:$PATH
export GRB_LICENSE_FILE=/groups/ms3625_gp/zw2994/gurobi.lic
export OMP_NUM_THREADS=4

LINE=$(awk -v i=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==i' "$BASE/inputs/carta_tasks.tsv")
VER=$(printf '%s\n' "$LINE" | cut -f1)
SAMPLE=$(printf '%s\n' "$LINE" | cut -f2)
[[ -n "$SAMPLE" ]] || { echo "no task for array index $SLURM_ARRAY_TASK_ID"; exit 1; }
TREE=$BASE/inputs/$VER/${SAMPLE}_tree_greedy.nwk
LABEL=$BASE/inputs/$VER/${SAMPLE}_labels.tsv
STATES=$BASE/inputs/states.txt
OUTDIR=$BASE/results/$VER/$SAMPLE
NST=$(awk -F'\t' 'NR>1 {print $2}' "$LABEL" | sort -u | wc -l)
KMAX=$(( NST < 8 ? NST : 8 ))
KLIST=$(seq 1 "$KMAX" | tr '\n' ' ')
mkdir -p "$OUTDIR"
if [[ -e "$OUTDIR/${SAMPLE}_progenitors_${KMAX}.txt" ]]; then echo "already complete: $OUTDIR"; exit 0; fi
printf '%s\t%s\n' "$TREE" "$LABEL" > "$OUTDIR/file_locations.tsv"

echo "[$(date)] START idx=$SLURM_ARRAY_TASK_ID version=$VER sample=$SAMPLE n_states=$NST klist=$KLIST"
cd "$OUTDIR"
# shellcheck disable=SC2086
carta --prefix "$OUTDIR/$SAMPLE" --file_locations "$OUTDIR/file_locations.tsv" --states_file "$STATES" \
      --klist $KLIST --normalize_method no_normalization --time_limit_sec 3000 --enforce_tree -m exact
echo "[$(date)] DONE version=$VER sample=$SAMPLE"
ls -la "$OUTDIR"

#!/bin/bash
#SBATCH --job-name=carta_fine12_greedy
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/carta_out/logs/carta_fine12_greedy_%A_%a.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/carta_out/logs/carta_fine12_greedy_%A_%a.err
#SBATCH --array=0-24%12
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --partition=cpu

# CARTA fine-12 (verbatim putative_state) x 25-sample x GREEDY-tree sweep
# ----------------------------------------------------------------------
# Same knobs as k4rec_greedy (which was validated); only differences:
#   - LABEL_DIR / STATES / OUTDIR point to fine12 inputs / fine12_greedy
#   - Array size 25 (all Npp53 tumor samples >=3 states after per-sample >=20-cell filter)
#   - Wall time 6h (12-state input; runtime scales ~n_states^2 * n_leaves)

set -euo pipefail

source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh || true
conda activate /groups/ms3625_gp/zw2994/mamba/envs/carta
export GRB_LICENSE_FILE=/groups/ms3625_gp/zw2994/gurobi.lic
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/carta/bin:$PATH

BASE=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing
SOLVER=greedy
SAMPLES=(JZ136 JZ201 JZ202 JZ203 JZ204 JZ205 JZ206 JZ207 JZ208 JZ209 JZ210 JZ211 JZ212 JZ213 JZ214 JZ224 JZ225 JZ226 JZ227 JZ228 JZ229 JZ230 JZ231 JZ232 JZ234)
STATES=$BASE/carta_out/inputs_fine12/states_fine12.txt
LABEL_DIR=$BASE/carta_out/inputs_fine12

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
[[ -z "$SAMPLE" ]] && { echo "ERROR: no sample for array idx $SLURM_ARRAY_TASK_ID"; exit 1; }

TREE=$LABEL_DIR/${SAMPLE}_tree_${SOLVER}_fine12.nwk
LABEL=$LABEL_DIR/${SAMPLE}_labels_fine12.tsv
OUTDIR=$BASE/carta_out/fine12_${SOLVER}/${SAMPLE}
mkdir -p "$OUTDIR"

FILEMAP=$OUTDIR/file_locations.tsv
echo -e "${TREE}\t${LABEL}" > "$FILEMAP"

echo "[$(date)] START . idx=${SLURM_ARRAY_TASK_ID} sample=$SAMPLE solver=$SOLVER"
echo "  tree   : $TREE"
echo "  labels : $LABEL"
echo "  states : $STATES"
echo "  outdir : $OUTDIR"
ls -la "$TREE" "$LABEL" "$STATES"
echo "  leaf tokens in tree (approx):"
head -c 4096 "$TREE" | tr ',' '\n' | wc -l
echo "  rows in labels tsv:"; wc -l "$LABEL"
echo "  fine-12 label distribution (sample):"; awk -F'\t' 'NR>1 {print $2}' "$LABEL" | sort | uniq -c

which carta

echo "[$(date)] launching carta ..."
cd "$OUTDIR"
carta \
  --prefix "$OUTDIR/${SAMPLE}" \
  --file_locations "$FILEMAP" \
  --states_file "$STATES" \
  --klist 1 2 3 4 5 6 \
  --normalize_method no_normalization \
  --time_limit_sec 3000 \
  --enforce_tree \
  -m exact 2>&1

echo "[$(date)] DONE . sample=$SAMPLE"
ls -la "$OUTDIR"

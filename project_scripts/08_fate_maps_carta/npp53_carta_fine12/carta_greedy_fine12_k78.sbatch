#!/bin/bash
#SBATCH --job-name=carta_fine12_greedy_k78
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/carta_out/logs/carta_fine12_greedy_k78_%A_%a.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/carta_out/logs/carta_fine12_greedy_k78_%A_%a.err
#SBATCH --array=0-24%12
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --partition=cpu

# CARTA fine-12 GREEDY-tree sweep -- EXTENSION for k=7 and k=8
# ---------------------------------------------------------------
# Adds k=7 and k=8 to the existing k=1..6 sweep (job 323320, all samples
# COMPLETED). Same inputs, same env, same solver settings; only differences:
#   - --klist 7 8
#   - Output goes to the SAME OUTDIR ($OUTDIR/${SAMPLE}_results_{7,8}.txt etc.)
#     CARTA writes prefix-tagged output per k, so k=1..6 files are NOT touched.
# ---------------------------------------------------------------
# n_states per sample (from fine12_greedy_summary.tsv):
#   n=3  : JZ203                                         -> k=7,8 will PLATEAU
#   n=4  : JZ136, JZ227                                  -> k=7,8 will PLATEAU
#   n=5  : JZ210,11,12,13,14,24,25,26,28,29,30,31,32,34  -> k=7,8 will PLATEAU
#   n=6  : JZ201, JZ202, JZ205, JZ209                    -> k=7,8 will PLATEAU
#   n=7  : JZ204, JZ207, JZ208                           -> k=7 informative, k=8 plateaus
#   n=9  : JZ206                                         -> k=7 and k=8 both informative
# Plateau runs still complete (fast for small samples), producing degenerate
# solutions that will be flagged post-hoc in compute_best_k.py.

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

echo "[$(date)] START k78 . idx=${SLURM_ARRAY_TASK_ID} sample=$SAMPLE solver=$SOLVER"
echo "  tree   : $TREE"
echo "  labels : $LABEL"
echo "  states : $STATES"
echo "  outdir : $OUTDIR"
ls -la "$TREE" "$LABEL" "$STATES"
echo "  existing k=1..6 files (should be preserved):"
ls "$OUTDIR" | grep -E "_(results|progenitors|nodeLabels)_[1-6]\.txt$" | wc -l

which carta

echo "[$(date)] launching carta --klist 7 8 ..."
cd "$OUTDIR"
carta \
  --prefix "$OUTDIR/${SAMPLE}" \
  --file_locations "$FILEMAP" \
  --states_file "$STATES" \
  --klist 7 8 \
  --normalize_method no_normalization \
  --time_limit_sec 3000 \
  --enforce_tree \
  -m exact 2>&1

echo "[$(date)] DONE k78 . sample=$SAMPLE"
echo "  k=7,8 outputs written:"
ls "$OUTDIR" | grep -E "_(results|progenitors|nodeLabels)_[78]\.txt$"
echo "  full outdir listing:"
ls -la "$OUTDIR"

#!/bin/bash
# =============================================================================
# run_local_from_step03.sh — run steps 03, 04, 06 on YOUR OWN machine (no SLURM),
# starting from the provided Cell Ranger / CellBender / allele-table folder.
#
# Before running:
#   1. bash envs/create_envs.sh seurat scvi cassiopeia     (once)
#   2. in config/config.sh, fill in and uncomment section 6 (OPTION B)
# Then, from the repository root:
#   bash run_local_from_step03.sh            # all three steps
#   bash run_local_from_step03.sh 04         # one step: 03 | 04 | 04b | 06
#   RES=0.7 bash run_local_from_step03.sh 04b   # re-run states at a chosen resolution
# Needs ~16-32 GB RAM for a full run (MJZ019 ~13k cells). No GPU needed.
# =============================================================================
set -eo pipefail
source config/config.sh
RUNNER="$(command -v mamba || command -v conda) run --no-capture-output -n"
STEP=${1:-all}
mkdir -p ${WORK}

if [[ $STEP == all || $STEP == 03 ]]; then
    echo ">>> 03 demux + QC + DoubletFinder (${RUN})"
    bash config/hashtags.sh
    ulimit -s unlimited 2>/dev/null || true
    ${RUNNER} tracing_seurat Rscript --max-ppsize=500000 03_demux_qc/demux_qc_doublets.R
fi
if [[ $STEP == all || $STEP == 04 ]]; then
    echo ">>> 04a tumour cells"
    ${RUNNER} tracing_scvi python -u 04_cell_states/04a_tumor_cells.py
fi
if [[ $STEP == all || $STEP == 04 || $STEP == 04b ]]; then
    echo ">>> 04b tumour cell states"
    ${RUNNER} tracing_scvi python -u 04_cell_states/04b_tumor_states.py ${RES:+--res ${RES}}
fi
if [[ $STEP == all || $STEP == 06 ]]; then
    IFS=',' read -ra SAMPLES <<< "${KEEP_SAMPLES}"
    for s in "${SAMPLES[@]}"; do
        echo ">>> 06 greedy tree ${s}"
        ${RUNNER} tracing_cassiopeia python -u 06_greedy_tree/build_greedy_tree.py --sample ${s} \
            || echo "    (${s} skipped — see message above)"
    done
fi
echo "done. outputs in ${WORK}"

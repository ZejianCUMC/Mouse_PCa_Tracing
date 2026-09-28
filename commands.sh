#!/bin/bash
# =============================================================================
# commands.sh — the whole tutorial as a command sheet.
# Run from the repository root on the c2b2 login node (only submits jobs; no compute here).
# Step by step (recommended the first time): copy one block at a time and check outputs.
# =============================================================================
set -euo pipefail
source config/config.sh
mkdir -p ${WORK} logs

# 0 · (Zejian only, once) build the subsampled demo FASTQs -> ${DEMO}
# sbatch 00_make_demo_data/make_demo_data.sbatch

# 1 · Cell Ranger (GEX + hashtags)                         ~1-2 h
J1=$(sbatch --parsable 01_cellranger/run_cellranger.sbatch)

# 5 · Cassiopeia on the lineage library (independent of 1-4) ~2-4 h
J5=$(sbatch --parsable 05_cassiopeia/run_cassiopeia.sbatch)

# 2 · CellBender (GPU), needs 1                            ~30 min
J2=$(sbatch --parsable --dependency=afterok:${J1} 02_cellbender/run_cellbender.sbatch)

# 3 · HTODemux + QC + DoubletFinder, needs 1 + 2           ~20 min
J3=$(sbatch --parsable --dependency=afterok:${J2} 03_demux_qc/run_demux_qc.sbatch)

# 4 · tumour cells + cell states, needs 3                  ~20-40 min
J4=$(sbatch --parsable --dependency=afterok:${J3} 04_cell_states/run_cell_states.sbatch)

# 6 · greedy tree per mouse, needs 4 + 5                   ~5 min
J6=$(sbatch --parsable --dependency=afterok:${J4}:${J5} 06_greedy_tree/run_greedy_tree.sbatch)

echo "submitted: cellranger=${J1} cassiopeia=${J5} cellbender=${J2} demux=${J3} states=${J4} tree=${J6}"
echo "watch:     squeue -u ${USER};  tail -f logs/<job>_<id>.out"

# After looking at resolution_metrics.csv / markers, re-run states at your chosen resolution:
#   sbatch --export=ALL,RES=0.7 04_cell_states/run_cell_states.sbatch
#   sbatch 06_greedy_tree/run_greedy_tree.sbatch

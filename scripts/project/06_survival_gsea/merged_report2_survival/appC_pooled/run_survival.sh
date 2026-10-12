#!/bin/bash
# Pooled arm, survival steps as run_meta_round.sh (km_scan sel without maxstat, then 09_null random-set calibration). 2026-10-11.
set -euo pipefail; cd "$(dirname "$0")/.."; export NCORES=12
/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript km_scan.R metaNorm_20261011_all_bcr/candidates_human.json metaNorm_20261011_all_bcr/km_scan.csv 0 sel && echo "KM_SCAN DONE"
/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript 09_null.R metaNorm_20261011_all_bcr && echo "NULL DONE"

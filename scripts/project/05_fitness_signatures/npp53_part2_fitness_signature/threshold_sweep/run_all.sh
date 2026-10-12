#!/usr/bin/env bash
# Threshold sweep driver — sequential execution of derive+hotspot+aucell → KM → xref.
set -euo pipefail
SWEEP=/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/fitness_signatures_survival/threshold_sweep
LOG=$SWEEP/logs

PY=/home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python
R=/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript

echo "[$(date +%F' '%T)] === Step 1-3: derive + hotspot + AUCell (T1..T4 sequentially) ==="
$PY $SWEEP/derive_and_hotspot.py 2>&1 | tee $LOG/step1_derive_hotspot_aucell.log

echo "[$(date +%F' '%T)] === Step 4: singscore KM sweep (T1..T4) ==="
$R $SWEEP/km_singscore.R 2>&1 | tee $LOG/step2_km_singscore.log

echo "[$(date +%F' '%T)] === Step 5: cross-reference + SUMMARY + README ==="
$PY $SWEEP/xref_m8like.py 2>&1 | tee $LOG/step3_xref_m8like.log

echo "[$(date +%F' '%T)] === Done ==="

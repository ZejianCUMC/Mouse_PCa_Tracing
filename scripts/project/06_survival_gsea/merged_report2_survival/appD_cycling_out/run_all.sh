#!/bin/bash
# Pooled Npp53 + Npp53Rb1 arm WITHOUT cycling cells (report 2 Appendix D, Zejian 2026-10-11): cells = cycOUT non-cycling LBI cells (36,485);
# same steps and settings as metaNorm_20261011_all_bcr: subclones (20) -> per-tumour log1p CP10k Blom/CR1 effects (21) -> meta lists, set all (22)
# -> Hotspot inputs (04) -> Hotspot k30 / top 1000 locally (05) -> candidates (06) -> km_scan sel (no maxstat) -> random-set null (09).
set -euo pipefail; cd "$(dirname "$0")/.."; A=metaNorm_20261011_all_cycout_bcr
PY=/home/jrchristin/zejian/miniforge3/envs/moscot/bin/python; HP=/home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python; RS=/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript
export PYTHONHASHSEED=42 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NCORES=12
$PY 20_meta_clades.py $A > $A/clades.log 2>&1 && echo "STEP clades done"
EXPR=cp10k NPROC=10 $PY 21_meta_effects.py $A blom > $A/effects.log 2>&1 && echo "STEP effects done"
META_SETS=all $PY 22_meta_lists.py $A cr > $A/meta_lists.log 2>&1 && echo "STEP lists done"
$PY 04_hotspot_inputs.py $A all > $A/hotspot_inputs.log 2>&1 && echo "STEP hotspot inputs done"
TOPN=1000 $HP 05_hotspot_run.py $A META_all_FSU_top1000_n3000 30 > $A/hotspot_k30_t1000.log 2>&1 && echo "STEP hotspot done"
$PY 06_candidates.py $A > $A/candidates.log 2>&1 && echo "STEP candidates done"
$RS km_scan.R $A/candidates_human.json $A/km_scan.csv 0 sel > $A/km_scan.log 2>&1 && echo "STEP km_scan done"
$RS 09_null.R $A > $A/null.log 2>&1 && echo "ALL DONE"

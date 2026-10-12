#!/bin/bash
# Normalized-expression re-run (Zejian 2026-10-10: "yes, u need to use normalized counts"). Same shared design as the raw-UMI arms
# (fitness Blom rank b / z-scored LBI z x per-tumour SE subclone-robust cr / OLS ols; same lists; Hotspot k15/30/50 x top500/1000 on raw
# counts as danb requires), for Npp53 and Npp53Rb1, with the per-tumour model on log1p CP10k (21_meta_effects.py EXPR=cp10k).
set -euo pipefail; cd "$(dirname "$0")"; PY=/home/jrchristin/zejian/miniforge3/envs/moscot/bin/python
for t in b z; do E=metaNorm_20261010_eff_$t; mkdir -p $E; cp metaR1/cells_n*_PRIVATE.csv metaR1/clades_PRIVATE.csv metaR1/tree_manifest.csv $E/
  [ -f $E/effects.csv.gz ] || EXPR=cp10k NPROC=10 $PY 21_meta_effects.py $E $([ $t = b ] && echo blom || echo z) > $E/effects.log 2>&1; done
for g in Npp53 Npp53Rb1; do for a in bcr zcr bols zols; do A=metaNorm_20261010_${g}_$a; mkdir -p $A; t=${a:0:1}; se=${a:1}
  cp metaR1/cells_n*_PRIVATE.csv metaR1/clades_PRIVATE.csv metaR1/tree_manifest.csv metaNorm_20261010_eff_$t/effects.csv.gz $A/
  META_SETS=$g $PY 22_meta_lists.py $A $se > $A/meta_lists.log 2>&1; done; done
echo "$(date) effects + lists done"
run() { KS="15 30 50" TOPS="500 1000" bash run_meta_round.sh $1 > $1/round.log 2>&1; }
for a in bcr zcr bols zols; do run metaNorm_20261010_Npp53_$a & sleep 20; run metaNorm_20261010_Npp53Rb1_$a & wait; done
echo "ALL NORM ARMS DONE $(date)"

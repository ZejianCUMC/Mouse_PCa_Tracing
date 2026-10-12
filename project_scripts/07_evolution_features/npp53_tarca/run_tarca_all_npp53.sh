#!/bin/bash
# TarCA on 25 Npp53 tumor trees (v8-greedy) using cell-state labels from sc_fitness_plasticity_v10.csv.
# Uses ancestorOT's runner (rMAST env — TarCA v0.1.0 already installed there).
set -e

ROOT=/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing
TREES=$ROOT/analysis/hpc_tree_rebuild/v8_trees/greedy
ANN_DIR=$ROOT/analysis/evolution_features_v2/tarca_npp53/inputs
OUT_BASE=$ROOT/analysis/evolution_features_v2/tarca_npp53/results
SCFP=$ROOT/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv
RUNNER=/media/jrchristin/Backup10TB1/agent/1-project/23-ancestorOT/ancestorOT/scripts/run_tarca.R
RSCRIPT=/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript
PY=/home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python

mkdir -p "$ANN_DIR" "$OUT_BASE"

echo "[prep] building per-sample annotation CSVs from scfp"
$PY - <<'PYEOF'
from pathlib import Path
import pandas as pd
SCFP='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT=Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/tarca_npp53/inputs')
df = pd.read_csv(SCFP)
for samp, g in df.groupby('sample'):
    (g.rename(columns={'cellBC':'TipLabel','putative_state':'TipAnn'})
       [['TipLabel','TipAnn']]
       .to_csv(OUT/f'{samp}_ann.csv', index=False))
print(f'wrote {df["sample"].nunique()} annotation CSVs')
PYEOF

echo "[run] TarCA per sample"
for tree in "$TREES"/*_tree_greedy.nwk; do
  samp=$(basename "$tree" _tree_greedy.nwk)
  ann=$ANN_DIR/${samp}_ann.csv
  out=$OUT_BASE/$samp
  mkdir -p "$out"
  if [ -f "$out/tarca_EffN.csv" ]; then
    echo "  [skip] $samp — already done"
    continue
  fi
  if [ ! -f "$ann" ]; then
    echo "  [skip] $samp — no annotation CSV"
    continue
  fi
  printf "[%s] " "$samp"
  $RSCRIPT "$RUNNER" --tree "$tree" --ann "$ann" --out "$out" 2>&1 | tail -3
  echo ""
done
echo "[done] TarCA finished on $(ls $OUT_BASE/*/tarca_EffN.csv 2>/dev/null | wc -l) / 25 samples"

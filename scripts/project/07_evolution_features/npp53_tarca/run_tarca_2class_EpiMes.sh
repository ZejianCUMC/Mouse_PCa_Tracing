#!/bin/bash
# TarCA on 25 Npp53 tumor trees with the 2-class Epithelial/Mesenchymal mapping
# (from major_cluster_stability analysis — k=2 is the data-supported hierarchy).
set -e
ROOT=/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing
TREES=$ROOT/analysis/hpc_tree_rebuild/v8_trees/greedy
ANN_DIR=$ROOT/analysis/evolution_features_v2/tarca_npp53/inputs_2class
OUT_BASE=$ROOT/analysis/evolution_features_v2/tarca_npp53/results_2class
RUNNER=/media/jrchristin/Backup10TB1/agent/1-project/23-ancestorOT/ancestorOT/scripts/run_tarca.R
RSCRIPT=/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript
PY=/home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python

mkdir -p "$ANN_DIR" "$OUT_BASE"

# 1. Build per-sample annotation CSVs with 2-class labels
$PY - <<'PYEOF'
from pathlib import Path
import pandas as pd
SCFP='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv'
OUT=Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/tarca_npp53/inputs_2class')
FINE_TO_2CLASS = {
    'Luminal-Distal':'Epithelial','Luminal-Proximal':'Epithelial','Luminal-ARdep':'Epithelial',
    'Luminal-CR':'Epithelial','Luminal-Persistent':'Epithelial','Luminal-Squamous':'Epithelial',
    'Basal-like':'Epithelial',
    'Partial-EMT':'Mesenchymal','Mes-Quiescent':'Mesenchymal','Mes-Collagen':'Mesenchymal',
    'Mes-Stemness':'Mesenchymal','Mes-myCAF':'Mesenchymal',
}
df = pd.read_csv(SCFP)
df['macro'] = df['putative_state'].map(FINE_TO_2CLASS)
n_unmapped = df['macro'].isna().sum()
if n_unmapped > 0:
    print(f'WARN: {n_unmapped} cells with unmapped fine state:', df.loc[df['macro'].isna(),'putative_state'].unique())
df = df.dropna(subset=['macro'])
for samp, g in df.groupby('sample'):
    (g.rename(columns={'cellBC':'TipLabel','macro':'TipAnn'})[['TipLabel','TipAnn']]
       .to_csv(OUT/f'{samp}_ann.csv', index=False))
print(f'wrote {df["sample"].nunique()} annotation CSVs with 2-class Epi/Mes labels')
print(df.groupby(['sample','macro']).size().unstack(fill_value=0).head(3))
PYEOF

# 2. Run TarCA per sample
for tree in "$TREES"/*_tree_greedy.nwk; do
  samp=$(basename "$tree" _tree_greedy.nwk)
  ann=$ANN_DIR/${samp}_ann.csv
  out=$OUT_BASE/$samp
  mkdir -p "$out"
  [ -f "$out/tarca_EffN.csv" ] && { echo "[skip] $samp"; continue; }
  [ ! -f "$ann" ] && { echo "[skip] $samp — no ann"; continue; }
  printf "[%s] " "$samp"
  $RSCRIPT "$RUNNER" --tree "$tree" --ann "$ann" --out "$out" 2>&1 | tail -2
done
echo "[done] 2-class TarCA finished on $(ls $OUT_BASE/*/tarca_EffN.csv 2>/dev/null | wc -l)/25 samples"

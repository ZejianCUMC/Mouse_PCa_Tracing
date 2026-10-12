#!/bin/bash
# Per sample: Greedy first; LBI, PhyloVelo, NJ, Spectral, MaxCut, Percolation after Greedy succeeds (they read its character matrix / tree).
set -euo pipefail
G=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/lineage_tracing/state_min20_rebuild_20261006
test ! -e $G/SUBMITTED.txt
for S in JZ267 JZ269 JZ274; do
  test ! -e $G/$S/cassiopeia/results
  cd $G/$S/viewer/logs
  g=$(sbatch --parsable -J ${S}m20_Greedy --export=ALL,SAMPLE=$S,STEP=Greedy -o %x_%j.out -e %x_%j.err $G/state_min20.sbatch); ids="Greedy=$g"
  for m in LBI PhyloVelo NJ Spectral MaxCut Percolation; do
    ids="$ids $m=$(sbatch --parsable -J ${S}m20_$m --dependency=afterok:$g --export=ALL,SAMPLE=$S,STEP=$m -o %x_%j.out -e %x_%j.err $G/state_min20.sbatch)"
  done
  echo "$S $ids" | tee -a $G/SUBMITTED.txt
done

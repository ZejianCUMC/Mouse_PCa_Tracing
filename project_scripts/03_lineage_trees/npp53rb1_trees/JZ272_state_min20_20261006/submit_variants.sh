#!/bin/bash
# Alternative trees for the adopted JZ272 state>=20 tree (Greedy already COMPLETE). Refuses to run twice.
set -euo pipefail
R=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/lineage_tracing/JZ272_state_min20_test_20261006
test -e $R/cassiopeia/results/JZ272/COMPLETE.json && test ! -e $R/SUBMITTED_variants.txt
cd $R/viewer/logs; ids=""
for m in NJ Spectral MaxCut Percolation; do ids="$ids $m=$(sbatch --parsable -J JZ272m20_$m --export=ALL,STEP=$m -o %x_%j.out -e %x_%j.err $R/jz272_min20.sbatch)"; done
echo "$ids" | tee $R/SUBMITTED_variants.txt

#!/bin/bash
# Greedy first; LBI and PhyloVelo after Greedy succeeds. Refuses to run twice.
set -euo pipefail
R=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/lineage_tracing/JZ272_state_min20_test_20261006
test ! -e $R/SUBMITTED.txt && test ! -e $R/cassiopeia/results
cd $R/viewer/logs
G=$(sbatch --parsable -J JZ272m20_Greedy --export=ALL,STEP=Greedy -o %x_%j.out -e %x_%j.err $R/jz272_min20.sbatch)
L=$(sbatch --parsable -J JZ272m20_LBI --dependency=afterok:$G --export=ALL,STEP=LBI -o %x_%j.out -e %x_%j.err $R/jz272_min20.sbatch)
P=$(sbatch --parsable -J JZ272m20_PhyloVelo --dependency=afterok:$G --export=ALL,STEP=PhyloVelo -o %x_%j.out -e %x_%j.err $R/jz272_min20.sbatch)
echo "Greedy=$G LBI=$L PhyloVelo=$P" | tee $R/SUBMITTED.txt

#!/bin/bash
# JZ201 without the isolated start cell CACTGTCCACATGCAG: 5 trees in parallel; jungle + PhyloVelo after Greedy.
set -euo pipefail
RUN=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/lineage_tracing/JZ201_excl_outlier_20261006
test ! -e $RUN/greedy && test ! -e $RUN/SUBMITTED.txt
cd $RUN/logs; S=JZ201
G=$(sbatch --parsable -J ${S}x_greedy --export=ALL,SAMPLE=$S,STEP=tree,SOLVER=greedy -o %x_%j.out -e %x_%j.err $RUN/scripts/jz201_downstream.sbatch)
O=""; for s in nj maxcut spectral percolation; do O="$O $(sbatch --parsable -J ${S}x_$s --export=ALL,SAMPLE=$S,STEP=tree,SOLVER=$s -o %x_%j.out -e %x_%j.err $RUN/scripts/jz201_downstream.sbatch)"; done
J=$(sbatch --parsable -J ${S}x_jungle    --dependency=afterok:$G --export=ALL,SAMPLE=$S,STEP=jungle    -o %x_%j.out -e %x_%j.err $RUN/scripts/jz201_downstream.sbatch)
P=$(sbatch --parsable -J ${S}x_phylovelo --dependency=afterok:$G --export=ALL,SAMPLE=$S,STEP=phylovelo -o %x_%j.out -e %x_%j.err $RUN/scripts/jz201_downstream.sbatch)
echo "$S greedy=$G other_trees=$O jungle=$J phylovelo=$P" | tee -a $RUN/SUBMITTED.txt

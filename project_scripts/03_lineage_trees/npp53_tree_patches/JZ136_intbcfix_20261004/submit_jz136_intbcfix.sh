#!/bin/bash
# 5 trees in parallel; jungle + PhyloVelo after the Greedy tree succeeds.
set -euo pipefail
RUN=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/lineage_tracing/JZ136_intbcfix_20261004
test ! -e $RUN/greedy && test ! -e $RUN/SUBMITTED.txt
cd $RUN/logs
G=$(sbatch --parsable -J jz136fix_greedy --export=ALL,STEP=tree,SOLVER=greedy -o %x_%j.out -e %x_%j.err $RUN/scripts/jz136_intbcfix.sbatch)
for s in nj maxcut spectral percolation; do
  sbatch --parsable -J jz136fix_$s --export=ALL,STEP=tree,SOLVER=$s -o %x_%j.out -e %x_%j.err $RUN/scripts/jz136_intbcfix.sbatch >> $RUN/SUBMITTED.txt
done
J=$(sbatch --parsable -J jz136fix_jungle    --dependency=afterok:$G --export=ALL,STEP=jungle    -o %x_%j.out -e %x_%j.err $RUN/scripts/jz136_intbcfix.sbatch)
P=$(sbatch --parsable -J jz136fix_phylovelo --dependency=afterok:$G --export=ALL,STEP=phylovelo -o %x_%j.out -e %x_%j.err $RUN/scripts/jz136_intbcfix.sbatch)
echo "greedy=$G jungle=$J phylovelo=$P" | tee -a $RUN/SUBMITTED.txt

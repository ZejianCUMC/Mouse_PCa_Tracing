#!/bin/bash
# Per sample: 5 trees in parallel; jungle + PhyloVelo after its Greedy tree succeeds. No automatic retries.
set -euo pipefail
RUN=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/lineage_tracing/MJZ019_rerun_trees_20261005
test ! -e $RUN/greedy && test ! -e $RUN/SUBMITTED.txt
cd $RUN/logs
for S in JZ201 JZ202 JZ203 JZ204 JZ205; do
  G=$(sbatch --parsable -J ${S}r_greedy --export=ALL,SAMPLE=$S,STEP=tree,SOLVER=greedy -o %x_%j.out -e %x_%j.err $RUN/scripts/jz201_205_downstream.sbatch)
  O=""
  for s in nj maxcut spectral percolation; do O="$O $(sbatch --parsable -J ${S}r_$s --export=ALL,SAMPLE=$S,STEP=tree,SOLVER=$s -o %x_%j.out -e %x_%j.err $RUN/scripts/jz201_205_downstream.sbatch)"; done
  J=$(sbatch --parsable -J ${S}r_jungle    --dependency=afterok:$G --export=ALL,SAMPLE=$S,STEP=jungle    -o %x_%j.out -e %x_%j.err $RUN/scripts/jz201_205_downstream.sbatch)
  P=$(sbatch --parsable -J ${S}r_phylovelo --dependency=afterok:$G --export=ALL,SAMPLE=$S,STEP=phylovelo -o %x_%j.out -e %x_%j.err $RUN/scripts/jz201_205_downstream.sbatch)
  echo "$S greedy=$G other_trees=$O jungle=$J phylovelo=$P" | tee -a $RUN/SUBMITTED.txt
done

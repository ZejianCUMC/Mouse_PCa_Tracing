#!/bin/bash
# Chain 3 pySCENIC stages: grn -> ctx -> aucell (afterok dependency).
# Run FROM the HPC work dir:
#   /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion
set -eo pipefail
cd "$(dirname "$0")"
mkdir -p /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion/logs

j1=$(sbatch --parsable 01_grn.sbatch)
j2=$(sbatch --parsable --dependency=afterok:$j1 02_ctx.sbatch)
j3=$(sbatch --parsable --dependency=afterok:$j2 03_aucell.sbatch)

echo "submitted chain:"
echo "  1 grn    : $j1"
echo "  2 ctx    : $j2  (afterok:$j1)"
echo "  3 aucell : $j3  (afterok:$j2)"
echo
echo "watch: squeue -u zw2994 -o '%.10i %.16j %.10T %.10M %R'"
echo "$j1 $j2 $j3" > .job_ids

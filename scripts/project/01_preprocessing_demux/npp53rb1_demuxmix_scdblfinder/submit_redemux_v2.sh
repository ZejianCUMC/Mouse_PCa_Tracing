#!/bin/bash
set -euo pipefail
S=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/scripts
cd /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/logs
for lib in MJZ025 MJZ027 MJZ028; do
  p=$(sbatch --parsable -J redemux_prep_$lib $S/redemux_prep_run.sh $lib)
  echo "$lib prep=$p"
  d=$(sbatch --parsable -J redemux_dmm_$lib --dependency=afterok:$p $S/redemux_demuxmix_run.sh $lib)
  echo "$lib demuxmix=$d (afterok:$p)"
done

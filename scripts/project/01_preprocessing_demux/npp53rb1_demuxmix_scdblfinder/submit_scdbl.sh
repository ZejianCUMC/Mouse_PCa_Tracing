#!/bin/bash
set -euo pipefail
S=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/scripts
cd /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/logs
env=$(sbatch --parsable $S/build_scdbl_env.sbatch)
echo "scdbl_env_build=$env"
for lib in MJZ025 MJZ026 MJZ027 MJZ028; do
  j=$(sbatch --parsable -J scdbl_$lib --dependency=afterok:$env $S/scdbl_run.sh $lib)
  echo "$lib scdbl=$j (afterok:$env)"
done

#!/bin/bash
set -euo pipefail
base=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing
out="$base/analysis/min400_20260921"
previous="$base/analysis/rebuild_with_legacy_20260919"
test ! -e "$out/submission.tsv" || { echo 'Refusing duplicate submission' >&2; exit 1; }
test ! -e "$out/merged_counts.h5ad" || { echo 'Refusing to overwrite prepared run' >&2; exit 1; }
test -d "$previous/inputs"
test -f "$previous/cellbender/MJZ005_filtered.h5"
mkdir -p "$out/scripts" "$base/logs"
test -e "$out/inputs" || ln -s "$previous/inputs" "$out/inputs"
test -e "$out/cellbender" || ln -s "$previous/cellbender" "$out/cellbender"
for name in rebuild_legacy_prepare.py rebuild_legacy_integrate.R rebuild_legacy_scvi.py rebuild_legacy_evaluate.py feature_name_mapping.py sample_retention.py min400_stage.sh; do
 cp "$base/scripts/$name" "$out/scripts/$name"
done
stage="$out/scripts/min400_stage.sh"
merge=$(sbatch --parsable --job-name=Rb1_400_integrate --cpus-per-task=8 --mem=128G --time=08:00:00 --output="$base/logs/%x_%j.out" --error="$base/logs/%x_%j.err" "$stage" prepare_integrate)
printf 'prepare_integrate\t%s\n' "$merge" | tee -a "$out/submission.tsv"
scvi=$(sbatch --parsable --dependency="afterok:$merge" --job-name=Rb1_400_scVI --partition=gpu --ntasks=1 --gpus-per-task=1 --cpus-per-task=8 --mem-per-gpu=96G --time=12:00:00 --output="$base/logs/%x_%j.out" --error="$base/logs/%x_%j.err" "$stage" scvi)
printf 'scvi\t%s\n' "$scvi" | tee -a "$out/submission.tsv"
evaluation=$(sbatch --parsable --dependency="afterok:$scvi" --job-name=Rb1_400_evaluate --cpus-per-task=8 --mem=128G --time=08:00:00 --output="$base/logs/%x_%j.out" --error="$base/logs/%x_%j.err" "$stage" evaluate)
printf 'evaluate\t%s\n' "$evaluation" | tee -a "$out/submission.tsv"

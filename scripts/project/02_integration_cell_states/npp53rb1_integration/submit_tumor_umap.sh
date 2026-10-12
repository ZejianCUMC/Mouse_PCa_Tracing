#!/bin/bash
set -euo pipefail
base=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing
out="$base/analysis/min400_20260921/tumor_umap_20260921"
test ! -e "$out/submission.tsv" || { echo 'Refusing duplicate submission' >&2; exit 1; }
mkdir -p "$out/scripts"
for name in tumor_umap_prepare.py tumor_umap_harmony.R tumor_umap_scvi.py tumor_umap_evaluate.py tumor_umap_stage.sh;do
 cp "$base/scripts/$name" "$out/scripts/$name"
done
stage="$out/scripts/tumor_umap_stage.sh"
prepare=$(sbatch --parsable --job-name=Rb1_tumor_prepare --cpus-per-task=8 --mem=96G --time=01:00:00 --output="$base/logs/%x_%j.out" --error="$base/logs/%x_%j.err" "$stage" prepare)
printf 'prepare\t%s\n' "$prepare" | tee -a "$out/submission.tsv"
harmony=$(sbatch --parsable --dependency="afterok:$prepare" --job-name=Rb1_tumor_Harmony --cpus-per-task=8 --mem=64G --time=01:00:00 --output="$base/logs/%x_%j.out" --error="$base/logs/%x_%j.err" "$stage" harmony)
printf 'harmony\t%s\n' "$harmony" | tee -a "$out/submission.tsv"
scvi=$(sbatch --parsable --dependency="afterok:$prepare" --job-name=Rb1_tumor_scVI --array=0-2 --partition=gpu --ntasks=1 --gpus-per-task=1 --cpus-per-task=8 --mem-per-gpu=64G --time=02:00:00 --output="$base/logs/%x_%A_%a.out" --error="$base/logs/%x_%A_%a.err" "$stage" scvi)
printf 'scvi_array\t%s\n' "$scvi" | tee -a "$out/submission.tsv"
evaluation=$(sbatch --parsable --dependency="afterok:$harmony:$scvi" --job-name=Rb1_tumor_UMAP --cpus-per-task=8 --mem=64G --time=01:00:00 --output="$base/logs/%x_%j.out" --error="$base/logs/%x_%j.err" "$stage" evaluate)
printf 'evaluate\t%s\n' "$evaluation" | tee -a "$out/submission.tsv"

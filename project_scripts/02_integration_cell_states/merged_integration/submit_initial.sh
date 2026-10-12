#!/bin/bash
set -euo pipefail
project=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp_tumor_merged
out="$project/analysis/merged_20260923_seed42"
test ! -e "$out/INITIAL_JOBS.tsv"
test -f "$out/PREPARED.json"
expected=7a472649cd18500a9eef7bfd49bfacdd38c53fc4865a1cdbcc223bdbedfb4d2a
actual=$(sha256sum "$out/shared_counts.h5ad")
test "${actual%% *}" = "$expected"
hvg=$(sbatch --parsable --job-name=NppM_HVG --cpus-per-task=8 --mem=64G --time=01:00:00 --output="$out/logs/%x_%j.out" --error="$out/logs/%x_%j.err" "$project/scripts/stage.sh" hvg)
printf 'stage\tjob_id\nHVG\t%s\n' "$hvg" > "$out/INITIAL_JOBS.tsv"
train=$(sbatch --parsable --dependency="afterok:$hvg" --kill-on-invalid-dep=yes --job-name=NppM_scVI --array=0-2%2 --partition=gpu --ntasks=1 --gpus-per-task=1 --cpus-per-task=8 --mem-per-gpu=64G --time=04:00:00 --output="$out/logs/%x_%A_%a.out" --error="$out/logs/%x_%A_%a.err" "$project/scripts/stage.sh" train)
printf 'scVI\t%s\n' "$train" >> "$out/INITIAL_JOBS.tsv"
cat "$out/INITIAL_JOBS.tsv"

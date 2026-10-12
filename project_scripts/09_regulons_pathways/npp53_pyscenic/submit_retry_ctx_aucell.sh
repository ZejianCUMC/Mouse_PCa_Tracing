#!/bin/bash
# Resubmit after 321730 (ctx) TIMEOUT+OOM. 321731 (aucell) is stuck PENDING
# with DependencyNeverSatisfied — left alone (queue-safe, no compute reserved).
# New pair: ctx (128G/8h) → aucell (afterok:new_ctx).
set -eo pipefail
cd "$(dirname "$0")"

j2=$(sbatch --parsable 02_ctx.sbatch)
j3=$(sbatch --parsable --dependency=afterok:$j2 03_aucell.sbatch)

echo "resubmitted:"
echo "  2 ctx    : $j2  (replaces 321730 TIMEOUT+OOM; mem bumped 64G->128G, wall 4:30->8:00)"
echo "  3 aucell : $j3  (afterok:$j2, replaces 321731 DependencyNeverSatisfied)"
echo "$j2 $j3" > .retry_job_ids

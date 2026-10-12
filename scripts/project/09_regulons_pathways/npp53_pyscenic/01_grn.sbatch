#!/bin/bash
#SBATCH --job-name=npp53_grn
#SBATCH --partition=cpu
#SBATCH --account=ms3625_gp
#SBATCH --cpus-per-task=30
#SBATCH --mem=128G
#SBATCH --time=09:30:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion/logs/01_grn.%j.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion/logs/01_grn.%j.err
set -eo pipefail
WORK=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion
RES=/groups/ms3625_gp/zw2994/Resource/pyscenic
ENV=/groups/ms3625_gp/zw2994/mamba/envs/pyscenic
export PATH=$ENV/bin:$PATH
echo "[$(date '+%F %T')] pyscenic grn (GRNBoost2, 30 workers, seed 777)"
$ENV/bin/pyscenic grn \
    "$WORK/npp53_tumor_no_cyc_ion.loom" \
    "$RES/mm_mgi_tfs.txt" \
    -o "$WORK/adjacencies.tsv" \
    --num_workers 30 --seed 777
echo "[$(date '+%F %T')] done: $(wc -l < $WORK/adjacencies.tsv) adjacency rows"

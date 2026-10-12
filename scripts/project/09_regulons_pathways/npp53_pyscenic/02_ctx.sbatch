#!/bin/bash
#SBATCH --job-name=npp53_ctx
#SBATCH --partition=cpu
#SBATCH --account=ms3625_gp
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion/logs/02_ctx.%j.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion/logs/02_ctx.%j.err
set -eo pipefail
WORK=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion
RES=/groups/ms3625_gp/zw2994/Resource/pyscenic
ENV=/groups/ms3625_gp/zw2994/mamba/envs/pyscenic
export PATH=$ENV/bin:$PATH
# NOTE: 128G (up from 64G) after 321730 OOM'd + timeout at 4:30h. Karthaus 13.6k
# cells needed 96G; 18.2k Npp53 tumor needs more. Wall bumped to 8h for headroom;
# actual computation typically <30 min once mem is sufficient (karthaus ~6 min).
DBS="$RES/mm9-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather \
     $RES/mm9-tss-centered-5kb-7species.mc9nr.genes_vs_motifs.rankings.feather \
     $RES/mm9-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather"
echo "[$(date '+%F %T')] pyscenic ctx (cisTarget, 3 mm9-mc9nr DBs, 20 workers, 128G mem)"
$ENV/bin/pyscenic ctx \
    "$WORK/adjacencies.tsv" $DBS \
    --annotations_fname "$RES/motifs-v9-nr.mgi-m0.001-o0.0.tbl" \
    --expression_mtx_fname "$WORK/npp53_tumor_no_cyc_ion.loom" \
    --output "$WORK/regulons.csv" \
    --num_workers 20
echo "[$(date '+%F %T')] done: $(wc -l < $WORK/regulons.csv) regulon rows"

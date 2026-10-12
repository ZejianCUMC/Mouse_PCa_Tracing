#!/bin/bash
#SBATCH --job-name=npp53_aucell
#SBATCH --partition=cpu
#SBATCH --account=ms3625_gp
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion/logs/03_aucell.%j.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion/logs/03_aucell.%j.err
set -eo pipefail
WORK=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion
ENV=/groups/ms3625_gp/zw2994/mamba/envs/pyscenic_aucell
export PATH=$ENV/bin:$PATH
echo "[$(date '+%F %T')] pyscenic aucell (env: pyscenic_aucell, 8 workers)"
$ENV/bin/pyscenic aucell \
    "$WORK/npp53_tumor_no_cyc_ion.loom" \
    "$WORK/regulons.csv" \
    --output "$WORK/npp53_tumor_pyscenic_output.loom" \
    --num_workers 8
echo "[$(date '+%F %T')] loom done: $(du -h $WORK/npp53_tumor_pyscenic_output.loom | cut -f1)"

echo "[$(date '+%F %T')] extract auc_mtx.csv (cells x regulons)"
$ENV/bin/python - <<'PY'
import loompy as lp, pandas as pd, numpy as np
loom = "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion/npp53_tumor_pyscenic_output.loom"
out  = "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/pysenic/npp53_tumor_no_cyc_ion/auc_mtx.csv"
with lp.connect(loom, mode="r", validate=False) as ds:
    reg_names = list(ds.ca["RegulonsAUC"].dtype.names)
    auc = pd.DataFrame(ds.ca["RegulonsAUC"], index=ds.ca["CellID"])
    auc.columns = reg_names
auc.to_csv(out)
print(f"[OK] wrote {out}  shape={auc.shape}")
med = auc.median(axis=0).sort_values(ascending=False)
print("top-5 regulons by median AUCell:")
for name, v in med.head(5).items():
    print(f"  {name}\t{v:.4f}")
n_act = sum(1 for c in auc.columns if c.endswith("(+)"))
n_rep = sum(1 for c in auc.columns if c.endswith("(-)"))
print(f"n_regulons = {auc.shape[1]}  ({n_act} activating + {n_rep} repressing)")
PY
echo "[$(date '+%F %T')] all done"

"""Preserve the warned run and stage a separate, hash-bound correction run."""
import hashlib,json,shutil
from pathlib import Path
project=Path(__file__).resolve().parents[1]
source=project/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42'
out=source/'choir_label_fix_seed42'
assert (source/'CHOIR_COMPLETE.json').exists()
out.mkdir(exist_ok=True)
assert not any(out.iterdir()),'Refuse overwriting an existing correction attempt'
inputs=['PREPARED.json','hvg_counts_PRIVATE.mtx.gz','features.txt',
        'barcodes_PRIVATE.txt','metadata_PRIVATE.csv','scvi_latent_PRIVATE.csv']
hashes={}
for name in inputs:
    # Workspace volume does not support symbolic links; preserve exact copies.
    shutil.copy2(source/name,out/name)
    hashes[name]=hashlib.sha256((source/name).read_bytes()).hexdigest()
(out/'CORRECTION_PLAN.json').write_text(json.dumps(dict(
    reason='Diagnosed CHOIR negative-level ID parser bug; 564 cells received temporary _NA labels during the warned run',
    correction='One project-local pruneTree regex literal accepts signed L levels; no installed package modification',
    preserve='Original warned CHOIR and preliminary AMI remain available; no publication acceptance',
    computation='Repeat CHOIR build+prune with unchanged seed42/settings/input hashes; do not repeat scVI, projections or Leiden',
    input_sha256=hashes,
    script_sha256={name:hashlib.sha256((project/'scripts'/name).read_bytes()).hexdigest()
                   for name in ['choir_negative_level_fix.R','tumor_states_choir.R']},
    resources=dict(cores=8,estimated_memory_GB=32,estimated_hours=1.5,timeout_hours=6,GPU_hours=0)
),indent=2))
print(out)

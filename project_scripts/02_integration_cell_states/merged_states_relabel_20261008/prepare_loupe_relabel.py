"""Loupe export with the approved 2026-10-08 labels: reuse the verified identical-cohort counts/projections from
reviewed10/loupe_export_oxidative (byte-identical copies) and replace only putative_state / final_cluster. Conda env: moscot."""
from pathlib import Path
import json, hashlib, shutil
import pandas as pd
N = Path(__file__).resolve().parent; src = N.parent / "loupe_export_oxidative"; out = N / "loupe_export"
assert not out.exists() or not any(out.iterdir()); out.mkdir(exist_ok=True)
sha = lambda f: hashlib.file_digest(open(f, "rb"), "sha256").hexdigest()
assert json.loads((src / "FINAL_EXPORT_VERIFIED.json").read_text())["status"] == "PASS"
for name in ["all_counts_PRIVATE.h5", "features.csv", "source_feature_mapping.csv", "COUNTS_COMPLETE.json", "umap_PRIVATE.csv", "tsne_PRIVATE.csv",
             "diffmap_PRIVATE.csv", "diffmap_all_components_PRIVATE.csv", "PROJECTIONS_COMPLETE.json"]:
    shutil.copy2(src / name, out / name); assert sha(src / name) == sha(out / name)
meta = pd.read_csv(src / "metadata_PRIVATE.csv", index_col=0); lab = pd.read_csv(N / "cell_labels_PRIVATE.csv", index_col=0)
assert meta.index.equals(lab.index) and (meta.final_cluster.astype(str) == lab.old_cluster).all() and (meta.putative_state == lab.old_state).all()
before = meta.drop(columns=["putative_state", "final_cluster"]).copy()
meta["putative_state"] = lab.new_state; meta["final_cluster"] = lab.new_cluster
assert meta.drop(columns=["putative_state", "final_cluster"]).equals(before)
meta.to_csv(out / "metadata_PRIVATE.csv")
print("PASS: identical counts/projections copied; labels replaced:", meta.groupby(["final_cluster", "putative_state"]).size().to_dict())

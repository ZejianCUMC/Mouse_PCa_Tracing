#!/usr/bin/env python
# Step 3a: per-tumor leaf -> merged state labels for ALL tree-matched cells (no per-state filter),
# same leaf/barcode matching as step1_inputs/build_inputs.py. Input for the PATH c10 check.
# Env: /home/jrchristin/miniconda3/envs/Cassiopeia/bin/python
from pathlib import Path
import pandas as pd
from ete3 import Tree
T1 = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
RV10 = T1 / "Npp_tumor_merged/analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42/approved_merges/reviewed10"
CELLS = RV10 / "fsu_m3_merged_20261008/cell_inputs_PRIVATE.csv"
MAN = RV10 / "fsu_m3_merged_20261008/opt5_corrected_20261009/metaNorm_20261011_all_cycout_bcr/tree_manifest.csv"
OUT = Path(__file__).resolve().parent / "labels_all"
cells = pd.read_csv(CELLS, usecols=["sample_id", "new_state", "bc16", "model", "treatment"])
rows = []
for r in pd.read_csv(MAN).itertuples():
    leaves = [l.name for l in Tree(str(T1 / r.tree), format=1).get_leaves()]
    lf = pd.DataFrame({"leaf": leaves, "bc16": [str(x).split("_")[-1].split("-")[0] for x in leaves]}).drop_duplicates("bc16")
    m = lf.merge(cells[cells.sample_id == r.tumour], on="bc16")
    m[["leaf", "new_state"]].to_csv(OUT / f"{r.tumour}.tsv", sep="\t", index=False)
    rows.append(dict(tumor=r.tumour, model=m.model.iloc[0], treatment=m.treatment.iloc[0], tree=str(T1 / r.tree),
                     n_leaves=len(leaves), n_matched=len(m), n_c10=int((m.new_state == "c10_Mes-Proliferation").sum())))
pd.DataFrame(rows).to_csv(OUT.parent / "tumors.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).groupby("model")[["n_matched", "n_c10"]].sum())

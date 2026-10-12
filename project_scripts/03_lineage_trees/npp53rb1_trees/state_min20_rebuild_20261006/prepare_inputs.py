# Prepare one sample for the per-sample state >= 20 cells rebuild (generalised from JZ272_state_min20_test_20261006/prepare_inputs.py).
# Usage: python prepare_inputs.py JZ267
import json, hashlib, shutil, numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse
RL = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/lineage_tracing")
import sys
OC, OV, S, MIN = RL / "cassiopeia_mjz026_028_20260925_seed42", RL / "tree_viewer_mjz026_028_20260925_seed42", sys.argv[1], 20
R = RL / f"state_min20_rebuild_20261006/{S}"; assert not R.exists()
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
meta = pd.read_csv(OC / "inputs" / S / "cell_metadata.csv", index_col=0); vc = meta.putative_state.value_counts()
keep_states = vc[vc >= MIN].index; keep = meta.index[meta.putative_state.isin(keep_states)]
removed = meta.loc[~meta.index.isin(keep)]
# tree inputs (cassiopeia run-dir layout)
ci = R / "cassiopeia/inputs" / S; ci.mkdir(parents=True); (R / "cassiopeia/scripts").mkdir(); (R / "cassiopeia/logs").mkdir()
a = pd.read_csv(OC / "inputs" / S / "allele_table.csv.gz", dtype={"cellBC": str, "intBC": str}); a = a[a.cellBC.isin(keep)]
a.to_csv(ci / "allele_table.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}); meta.loc[keep].to_csv(ci / "cell_metadata.csv")
qc = json.loads((OC / "inputs" / S / "INPUT_QC.json").read_text())
qc.update(matched_cells=len(keep), allele_rows=len(a), cell_intbc_pairs=len(a), allele_table_sha256=sha(ci / "allele_table.csv.gz"), metadata_sha256=sha(ci / "cell_metadata.csv"),
          no_new_cell_exclusions=False, state_filter=f"Zejian 2026-10-06 rebuild (JZ272 rule): keep states with >= {MIN} cells in {S}; removed {len(removed)} cells")
(ci / "INPUT_QC.json").write_text(json.dumps(qc, indent=2)); (R / "cassiopeia/PREPARED.json").write_text(json.dumps(dict(samples=[S]), indent=2))
shutil.copy2(OC / "scripts/run_cassiopeia_eight.py", R / "cassiopeia/scripts/")
# PhyloVelo / LBI inputs (tree_viewer root layout)
vi = R / "viewer/inputs" / S; vi.mkdir(parents=True); (R / "viewer/scripts").mkdir(); (R / "viewer/logs").mkdir()
cells = pd.read_csv(OV / "inputs" / S / "cells.csv", index_col=0).index; m = cells.isin(keep)
x = sparse.load_npz(OV / "inputs" / S / "counts.npz").tocsr()[np.where(m)[0]]; sparse.save_npz(vi / "counts.npz", x)
pd.DataFrame(index=pd.Index(cells[m], name="cell")).to_csv(vi / "cells.csv"); shutil.copy2(OV / "inputs" / S / "genes.csv", vi / "genes.csv")
u = pd.read_csv(OV / "inputs" / S / "umap.csv", index_col=0); u.loc[cells[m]].to_csv(vi / "umap.csv")
(R / "viewer/PREPARED.json").write_text(json.dumps(dict(samples=[S], inputs=[dict(sample=S, cells=int(m.sum()), files={f: sha(vi / f) for f in ["counts.npz", "cells.csv", "genes.csv", "umap.csv"]})]), indent=2))
src = (OV / "scripts/run_tree_viewer_analysis.py").read_text()
old = "root.parent/'cassiopeia_mjz026_028_20260925_seed42'"; assert src.count(old) == 1
(R / "viewer/scripts/run_tree_viewer_analysis.py").write_text(src.replace(old, "root.parent/'cassiopeia'   # state>=20 rebuild (only change)"))
removed[["putative_state"]].to_csv(R / "REMOVED_CELLS.csv")
print("kept", len(keep), "removed", len(removed)); print(removed.putative_state.value_counts().to_string())

"""Independent check of the relabelled Loupe export (adapted from scripts/verify_merged_cloupe.py --reviewed10 --oxidative).
Labels are checked against the approved 2026-10-08 mapping instead of the old graph labels; counts/projections against the
previously verified export (all_counts_PRIVATE.h5 is a byte-identical copy). Conda env: moscot."""
from pathlib import Path
import json, hashlib
import numpy as np, pandas as pd, h5py, anndata as ad
from scipy.sparse import csc_matrix
N = Path(__file__).resolve().parent; d = N / "loupe_export"; old = N.parent / "loupe_export_oxidative"; n = 55133
sha = lambda f: hashlib.file_digest(open(f, "rb"), "sha256").hexdigest()
done = json.loads((d / "CLOUPE_COMPLETE.json").read_text()); counts = json.loads((d / "COUNTS_COMPLETE.json").read_text()); pr = json.loads((d / "PROJECTIONS_COMPLETE.json").read_text())
for f in ["all_counts_PRIVATE.h5", "features.csv", "umap_PRIVATE.csv", "tsne_PRIVATE.csv", "diffmap_PRIVATE.csv"]: assert sha(d / f) == sha(old / f), f
assert json.loads((old / "FINAL_EXPORT_VERIFIED.json").read_text())["status"] == "PASS"
meta = pd.read_csv(d / "metadata_PRIVATE.csv", index_col=0); f = pd.read_csv(d / "features.csv"); smap = pd.read_csv(N / "STATE_MAP.csv").set_index("old_cluster")
graph = N.parent / "named_graph_PRIVATE.h5ad"; assert sha(graph) == pr["source_sha256"]
a = ad.read_h5ad(graph); assert a.obs_names.equals(meta.index)
np.testing.assert_array_equal(a.obs.final_cluster.astype(str).map(smap.new_cluster), meta.final_cluster.astype(str))
np.testing.assert_array_equal(a.obs.final_cluster.astype(str).map(smap.new_state), meta.putative_state.astype(str))
assert done["cells"] == counts["cells"] == n and done["features"] == counts["features"] == len(f) and done["states"] == 10
coords = {k: pd.read_csv(d / v, index_col=0) for k, v in [("UMAP", "umap_PRIVATE.csv"), ("tSNE", "tsne_PRIVATE.csv"), ("Diffmap", "diffmap_PRIVATE.csv")]}
np.testing.assert_allclose(coords["UMAP"], a.obsm["X_umap"], rtol=1e-7, atol=1e-7)
S = lambda ds: [v.decode() if isinstance(v, bytes) else str(v) for v in np.asarray(ds).ravel()]
with h5py.File(d / "all_counts_PRIVATE.h5") as src, h5py.File(d / "louper_input_PRIVATE.h5") as ex:
    assert S(ex["matrix/barcodes"]) == meta.index.tolist() and set(ex["projections"]) == {"UMAP", "tSNE", "Diffmap"}
    for k, x in coords.items():
        y = ex[f"projections/{k}/data"][:]; y = y.T if y.shape == (2, n) else y; np.testing.assert_allclose(y, x, rtol=1e-13, atol=1e-13)
    for key in ["id", "name"]: assert S(src[f"matrix/features/{key}"]) == S(ex[f"matrix/features/{key}"])
    assert np.array_equal(src["matrix/shape"][:], ex["matrix/shape"][:]); ptr = src["matrix/indptr"][:]; assert np.array_equal(ptr, ex["matrix/indptr"][:])
    for first in range(0, n, 512):
        last = min(first + 512, n); lo, hi = map(int, [ptr[first], ptr[last]]); xx = []
        for h in [src, ex]:
            z = csc_matrix((h["matrix/data"][lo:hi], h["matrix/indices"][lo:hi], ptr[first:last + 1] - lo), shape=(len(f), last - first)); z.sort_indices(); xx.append(z)
        assert np.array_equal(xx[0].indices, xx[1].indices) and np.array_equal(xx[0].data, xx[1].data)
    cats = {}
    for key, col in [("Putative_cell_state", "putative_state"), ("Cluster_ID", "final_cluster"), ("Model", "model"), ("Treatment", "treatment"), ("Run_ID", "run_id"), ("Sample_ID", "sample_id")]:
        names = np.array(S(ex[f"clusters/{key}/group_names"])); ids = ex[f"clusters/{key}/assignments"][:].ravel().astype(int)
        np.testing.assert_array_equal(names[ids], meta[col].astype(str)); cats[key] = names.tolist()
cl = d / "Npp_tumor_merged_final_states_all_genes.cloupe"; assert cl.stat().st_size == done["bytes"]
rec = dict(status="PASS", cells=n, features=len(f), states=10, state_order_in_loupe=cats["Putative_cell_state"], cluster_ids=cats["Cluster_ID"],
           counts_identical_to_verified_export=True, louper_input_counts_exact=True, projections_exact=True, labels_match_approved_mapping=True,
           loupe_browser_gui_open_not_tested=True, cloupe_sha256=sha(cl), bytes=cl.stat().st_size)
(d / "FINAL_EXPORT_VERIFIED.json").write_text(json.dumps(rec, indent=2)); print(json.dumps(rec, indent=2))

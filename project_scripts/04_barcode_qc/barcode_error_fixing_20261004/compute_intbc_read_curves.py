# =============================================================================
# Per-sample intBC read-count curves + knee (cf. mESC clone 1D5 "10 major intBCs" plot)
# =============================================================================
# Description : For each of our tree datasets, total reads (sum of readCount) per intBC in the tree-input allele
#               table (tree cells only; all intBCs, incl. error/low intBCs, so the tail is visible), ranked high to
#               low. Knee = largest step down in log10 reads between consecutive ranks -> number of "major" intBCs.
# Input       : Npp53 allele_tables/per_sample/<S>_allele_table.csv (v8 tree input, Ionocyte/Cycling removed);
#               Np53Rb1 frozen Cassiopeia inputs (pre-fix inputs for JZ135_n3000/JZ137/JZ138 so error intBCs show)
# Output      : data/intbc_read_curves.csv (sample, rank, intBC, reads, umis, cells, knee_k, ...)
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python
# Validated   : 2026-10-04
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import json
from pathlib import Path
import numpy as np, pandas as pd

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
B = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
N, RL = B / "Npp53_tracing", B / "Npp53Rb1_Tracing/lineage_tracing"
DEFAULTS = {
    "npp53_dir": N / "lineage_tracing/allele_tables/per_sample",
    "npp53_samples": ["JZ136"] + [f"JZ{i}" for i in list(range(201, 215)) + list(range(224, 233)) + [234]],
    "rb1_inputs": {**{s: RL / "cassiopeia_eight_20260923_seed42/inputs" / s for s in ["JZ137", "JZ138", "JZ256", "JZ259", "JZ261", "JZ262", "JZ263"]},
                   "JZ135_n3000": RL / "tree_viewer_20260923_seed42/subsample_cassiopeia/inputs/JZ135_n3000",
                   **{s: RL / "cassiopeia_mjz026_028_20260925_seed42/inputs" / s for s in
                      ["JZ265", "JZ266", "JZ267", "JZ268", "JZ269", "JZ270", "JZ271", "JZ272", "JZ273", "JZ274", "JZ276", "JZ277", "JZ279", "JZ280"]}},
    "exclude_states": {"Ionocyte", "Cycling"},   # as the v8 Npp53 tree inputs
    "knee_min": 3,
    "out": Path(__file__).resolve().parent / "data/intbc_read_curves.csv",
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def knee(reads, kmin):
    """Major-intBC knee = largest drop in log10(reads + 1) between consecutive ranks (rank >= kmin); returns number of
    intBCs before the drop. Chosen over the chord rule: our real intBCs decline gradually in reads, so the meaningful
    break is the step down to error / low-level intBCs (e.g. exactly 25 for JZ137/JZ138)."""
    y = np.log10(np.asarray(reads, float) + 1)
    if len(y) <= kmin: return len(y)
    dy = y[:-1] - y[1:]; dy[:kmin - 1] = -np.inf
    return int(np.argmax(dy) + 1)

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}; tabs = {}
    for s in p["npp53_samples"]:
        a = pd.read_csv(p["npp53_dir"] / f"{s}_allele_table.csv", dtype={"cellBC": str, "intBC": str})
        tabs[s] = ("Npp53", a[~a.putative_state.isin(p["exclude_states"])])
    for s, d in p["rb1_inputs"].items():
        tabs[s] = ("Np53Rb1", pd.read_csv(d / "allele_table.csv.gz", dtype={"cellBC": str, "intBC": str}))
    rows = []
    for s, (proj, a) in tabs.items():
        g = a.groupby("intBC").agg(reads=("readCount", "sum"), umis=("UMI", "sum"), cells=("cellBC", "nunique")).sort_values("reads", ascending=False)
        k = knee(g.reads.values, p["knee_min"])
        gap = float(np.max(np.diff(-np.log10(g.reads.values.astype(float) + 1)))) if len(g) > 1 else 0.0
        for r, (ib, x) in enumerate(g.iterrows(), 1):
            rows.append(dict(sample=s, project=proj, rank=r, intBC=ib, reads=int(x.reads), umis=int(x.umis), cells=int(x.cells),
                             n_cells=a.cellBC.nunique(), n_intBCs=len(g), knee_k=k, knee_gap_log10=round(gap, 3), is_major=r <= k))
    d = pd.DataFrame(rows); p["out"].parent.mkdir(exist_ok=True); d.to_csv(p["out"], index=False)
    return {"output_path": str(p["out"]), "params_used": {k: str(v) for k, v in p.items()},
            "summary": d.drop_duplicates("sample")[["sample", "n_intBCs", "knee_k"]].to_dict("records")}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    r = run_pipeline(); print(pd.DataFrame(r["summary"]).to_string(index=False))

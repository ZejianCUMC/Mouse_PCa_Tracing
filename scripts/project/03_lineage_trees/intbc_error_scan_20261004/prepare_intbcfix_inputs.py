# =============================================================================
# Error-intBC removal: tree inputs for JZ136 / JZ137 / JZ138 / JZ135_n1000-n3000
# =============================================================================
# Description : Zejian 2026-10-04: rebuild these trees after REMOVING error intBCs. Error intBC = an
#               intBC whose detected-cell count is <= 30 % of a more-detected intBC within Levenshtein
#               <= 2 in the same input table (sequencing/PCR variants of the 25 founder intBCs; see
#               intbc_scan.csv). All rows of an error intBC are dropped; nothing else changes (cells,
#               metadata, read/UMI values, row order). Cells left with no rows are reported (and would
#               stop the run - none expected).
# Input       : Np53Rb1 frozen tree inputs (cassiopeia_eight_20260923_seed42/inputs/{JZ137,JZ138},
#               tree_viewer_20260923_seed42/subsample_cassiopeia/inputs/JZ135_n*), Npp53
#               lineage_tracing/allele_tables/per_sample/JZ136_allele_table.csv (v8 tree input)
# Output      : Np53Rb1 lineage_tracing/intbcfix_20261004_seed42/cassiopeia/{inputs/<S>/, PREPARED.json}
#               Npp53   lineage_tracing/JZ136_intbcfix_20261004/input/JZ136_allele_table.csv
#               + ERROR_INTBCS_REMOVED.csv / INTBCFIX_SUMMARY.csv in each
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python
# Validated   : 2026-10-04
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import json, hashlib, shutil
from pathlib import Path
import pandas as pd

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
MB = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
RL = MB / "Npp53Rb1_Tracing/lineage_tracing"
DEFAULTS = {
    "max_lev": 2,            # error intBC: <= 2 edits from a more-detected intBC ...
    "max_ratio": 0.30,       # ... detected in <= 30 % as many cells (all 108 scanned errors <= 0.298)
    "rb1_sources": {"JZ137": RL / "cassiopeia_eight_20260923_seed42/inputs/JZ137",
                    "JZ138": RL / "cassiopeia_eight_20260923_seed42/inputs/JZ138",
                    **{f"JZ135_n{n}": RL / f"tree_viewer_20260923_seed42/subsample_cassiopeia/inputs/JZ135_n{n}" for n in (1000, 2000, 3000)}},
    "rb1_out": RL / "intbcfix_20261004_seed42/cassiopeia",
    "npp53_source": MB / "Npp53_tracing/lineage_tracing/allele_tables/per_sample/JZ136_allele_table.csv",
    "npp53_out": MB / "Npp53_tracing/lineage_tracing/JZ136_intbcfix_20261004/input",
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""): h.update(b)
    return h.hexdigest()

def lev(a, b):
    d = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        prev, d[0] = d[0], i
        for j, y in enumerate(b, 1):
            prev, d[j] = d[j], min(d[j] + 1, d[j - 1] + 1, prev + (x != y))
    return d[-1]

def find_errors(a, p, cell_col="cellBC"):
    det = a.groupby("intBC")[cell_col].nunique().sort_values(ascending=False)
    out = []
    for g, n in det.items():
        cand = [(h, lev(g, h)) for h, m in det.items() if m > n]
        cand = [(h, d) for h, d in cand if d <= p["max_lev"] and n <= p["max_ratio"] * det[h]]
        if cand:
            h, d = min(cand, key=lambda x: (x[1], -det[x[0]]))
            out.append(dict(intBC=g, cells=int(n), parent_intBC=h, distance=d, parent_cells=int(det[h]), ratio=round(n / det[h], 4)))
    return pd.DataFrame(out, columns=["intBC", "cells", "parent_intBC", "distance", "parent_cells", "ratio"]), det

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    summary, errs = [], []
    # ---- Np53Rb1 (frozen Cassiopeia inputs; same format as run_cassiopeia_eight.py expects) ----
    out = p["rb1_out"]; (out / "inputs").mkdir(parents=True, exist_ok=False)
    records = []
    for s, src in p["rb1_sources"].items():
        a = pd.read_csv(src / "allele_table.csv.gz", dtype=str, keep_default_na=False)   # strings in = strings out
        meta = pd.read_csv(src / "cell_metadata.csv", index_col=0)
        e, det = find_errors(a, p)
        b = a[~a.intBC.isin(set(e.intBC))]
        lost = sorted(set(a.cellBC) - set(b.cellBC))
        assert not lost, f"{s}: {len(lost)} cells would lose all rows"
        dest = out / "inputs" / s; dest.mkdir()
        b.to_csv(dest / "allele_table.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
        shutil.copy2(src / "cell_metadata.csv", dest / "cell_metadata.csv")
        q = json.loads((src / "INPUT_QC.json").read_text())
        q.update(allele_rows=int(len(b)), allele_table_sha256=sha(dest / "allele_table.csv.gz"), metadata_sha256=sha(dest / "cell_metadata.csv"),
                 matched_cells=int(b.cellBC.nunique()), intbcfix_source=str(src / "allele_table.csv.gz"),
                 intbcfix_source_sha256=sha(src / "allele_table.csv.gz"), error_intbcs_removed=int(len(e)),
                 error_rows_removed=int(len(a) - len(b)), intbcs_kept=int(b.intBC.nunique()),
                 intbcfix_rule=f"drop intBC with <= {p['max_ratio']:.0%} of the cells of a more-detected intBC within Levenshtein <= {p['max_lev']}")
        (dest / "INPUT_QC.json").write_text(json.dumps(q, indent=2)); records.append(q)
        errs.append(e.assign(project="Np53Rb1", sample=s))
        summary.append(dict(project="Np53Rb1", sample=s, cells=int(len(meta)), rows_in=int(len(a)), rows_out=int(len(b)),
                            intBCs_in=int(det.size), error_intBCs_removed=int(len(e)), intBCs_out=int(b.intBC.nunique())))
        print(summary[-1], flush=True)
    (out / "PREPARED.json").write_text(json.dumps(dict(
        samples=list(p["rb1_sources"]), seed=42, solver="VanillaGreedySolver", allele_rep_thresh=.99, collapse_duplicates=True,
        mutation_priors=None, collapse_mutationless_edges=False, held_runs=[], records=records,
        purpose="Zejian 2026-10-04: rebuild JZ137/JZ138/JZ135 subsample trees after removing error intBCs (all else identical)"), indent=2))
    # ---- Npp53 JZ136 (v8 input format: per-sample csv, default read in the v8 script) ----
    out = p["npp53_out"]; out.mkdir(parents=True, exist_ok=False)
    a = pd.read_csv(p["npp53_source"], dtype={"cellBC": str, "intBC": str}, keep_default_na=False)
    a_v8 = a[~a.putative_state.isin({"Ionocyte", "Cycling"})]          # detection counted on the v8 tree cells
    e, det = find_errors(a_v8, p)
    b = a[~a.intBC.isin(set(e.intBC))]
    assert not (set(a_v8.cellBC) - set(b.cellBC)), "JZ136: cells would lose all rows"
    raw = pd.read_csv(p["npp53_source"], dtype=str, keep_default_na=False)   # write back byte-faithfully (same columns/order/strings)
    raw[~raw.intBC.isin(set(e.intBC))].to_csv(out / "JZ136_allele_table.csv", index=False)
    errs.append(e.assign(project="Npp53", sample="JZ136"))
    summary.append(dict(project="Npp53", sample="JZ136", cells=int(a_v8.cellBC.nunique()), rows_in=int(len(a)), rows_out=int(len(b)),
                        intBCs_in=int(det.size), error_intBCs_removed=int(len(e)), intBCs_out=int(b.intBC.nunique())))
    print(summary[-1], flush=True)
    E, S = pd.concat(errs), pd.DataFrame(summary)
    for o in [p["rb1_out"].parent, p["npp53_out"].parent]:
        E.to_csv(o / "ERROR_INTBCS_REMOVED.csv", index=False); S.to_csv(o / "INTBCFIX_SUMMARY.csv", index=False)
    return {"output_path": [str(p["rb1_out"]), str(p["npp53_out"])], "params_used": {k: str(v) for k, v in p.items()}, "summary": summary}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    run_pipeline()

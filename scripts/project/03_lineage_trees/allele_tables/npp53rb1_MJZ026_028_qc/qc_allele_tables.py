# =============================================================================
# Allele-table QC: MJZ026-028 vs previous TSL libraries
# =============================================================================
# Description : Library-level QC of Cassiopeia allele tables (cells, intBCs, UMIs/cell,
#               intBCs/cell, reads/UMI, lineage groups, per-site editing) plus overlap
#               with the current Np53Rb1 final RNA cohort (per mouse).
# Input       : per-library *_allele_table.csv; final-cohort metadata_PRIVATE.csv
# Output      : qc_library_summary.csv, qc_per_mouse_overlap.csv (next to this script)
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/scvi/bin/python
# Validated   : 2026-09-25
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import json
from pathlib import Path
import numpy as np
import pandas as pd

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
TR = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
HERE = Path(__file__).resolve().parent
DEFAULTS = {
    "new_tables": {L: str(HERE / f"{L}_allele_table.csv") for L in ["MJZ026", "MJZ027", "MJZ028"]},  # this run
    "ref_tables": {L: str(TR / f"Npp53_tracing/lineage_tracing/allele_tables/per_library/{L}_allele_table.csv")
                   for L in ["MJZ019", "MJZ021", "MJZ022", "MJZ023", "MJZ024", "MJZ025"]},            # earlier libraries
    "rna_meta": str(TR / "Npp53Rb1_Tracing/analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/"
                         "exclude_current_c17_c19_seed42/metadata_PRIVATE.csv"),                       # final Np53Rb1 cohort (40,141 cells)
    "sites": ["r1", "r2", "r3"],                # Cas9 cut sites
    "unedited_token": "None",                   # site string containing this = unedited (e.g. '[None]')
    "out_dir": str(HERE),
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def library_qc(name, path, p):
    a = pd.read_csv(path, keep_default_na=False)
    per_cell = a.groupby("cellBC").agg(umis=("UMI", "sum"), intbcs=("intBC", "nunique"), reads=("readCount", "sum"))
    # dominant allele per (cell, intBC) by UMI, then per-site editing (covered sites only)
    dom = a.sort_values("UMI", ascending=False).drop_duplicates(["cellBC", "intBC"])
    ed = {}
    for s in p["sites"]:
        v = dom[s].astype(str)
        cov = v.str.len() > 0
        ed[f"pct_edited_{s}"] = round(100 * (~v[cov].str.contains(p["unedited_token"])).mean(), 1)
    multi = (a.groupby(["cellBC", "intBC"]).size() > 1).sum()
    return {
        "library": name, "rows": len(a), "cells": a.cellBC.nunique(), "intBCs": a.intBC.nunique(),
        "lineage_groups": a.lineageGrp.nunique(),
        "median_UMIs_per_cell": float(per_cell.umis.median()), "median_intBCs_per_cell": float(per_cell.intbcs.median()),
        "median_reads_per_UMI": round(float((per_cell.reads / per_cell.umis).median()), 2),
        "multi_allele_cell_intBC_pairs": int(multi), **ed,
    }, a

def mouse_overlap(name, a, meta):
    m = meta[meta.run_id == name]
    if m.empty:
        return []
    cb = set(a.cellBC.unique())
    rows = []
    for (sid, trt), g in m.groupby(["sample_id", "treatment"]):
        n = len(g); k = int(g.cellBC.isin(cb).sum())
        rows.append({"library": name, "sample_id": sid, "treatment": trt, "final_rna_cells": n,
                     "with_allele_data": k, "pct_matched": round(100 * k / n, 1)})
    return rows

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    out = Path(output_dir or p["out_dir"])
    meta = pd.read_csv(p["rna_meta"], usecols=["sample_id", "run_id", "treatment", "source_barcode"])
    meta["cellBC"] = meta.source_barcode.str.replace(r"-\d+$", "", regex=True)
    lib_rows, mouse_rows = [], []
    for grp, tables in [("new", p["new_tables"]), ("reference", p["ref_tables"])]:
        for name, path in tables.items():
            if not Path(path).exists():
                print(f"[skip] {name}: {path} missing"); continue
            r, a = library_qc(name, path, p); r["set"] = grp
            lib_rows.append(r); mouse_rows += mouse_overlap(name, a, meta)
            print(json.dumps(r))
    lib = pd.DataFrame(lib_rows); mo = pd.DataFrame(mouse_rows)
    lib.to_csv(out / "qc_library_summary.csv", index=False); mo.to_csv(out / "qc_per_mouse_overlap.csv", index=False)
    return {"output_path": str(out), "params_used": p, "summary": {"libraries": len(lib), "mice": len(mo)}}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    res = run_pipeline()
    print(res["summary"])

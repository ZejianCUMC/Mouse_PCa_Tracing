# =============================================================================
# Per-intBC information/missingness: ours vs KP-Tracer (Yang et al., Cell 2022)
# =============================================================================
# Description : For every dataset (our 47 tree datasets; KP-Tracer tumours with >= MIN_CELLS cells) and every
#               intBC: fraction of cells detected, edited and informative, and the missing fraction of its
#               cut sites. Informative = cell carries an indel allele shared by >=2 cells and by <95% of the
#               cells detected at that site. Then, per dataset, the missing-site rate when keeping the best
#               10, best 15, or a per-dataset knee (inflection point of the sorted informative-fraction curve)
#               of the most informative intBCs.
# Input       : our v3 display matrices (-1 missing, 0 uncut, >0 edit, tuple = ambiguous -> observed);
#               KPTracer.alleleTable.FINAL.txt (blank site = missing, "[None]" = uncut)
# Output      : data/intbc_info_both.csv, data/selection_missing_both.csv
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python
# Validated   : 2026-10-04
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import glob, json, os
from pathlib import Path
import numpy as np, pandas as pd

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
B = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
DEFAULTS = {
    "ours": [B / "Npp53_tracing/figures/barcode_heatmaps_v3_20261004/_display_matrices",
             B / "Npp53Rb1_Tracing/figures/barcode_heatmaps_v3_20261004/_display_matrices"],
    "kp_table": Path("/media/jrchristin/Backup10TB1/agent/1-project/23-ancestorOT/data/yang2022/unpack/KPTracer-Data/KPTracer.alleleTable.FINAL.txt"),
    "min_cells": 200,          # KP-Tracer tumours kept (our smallest tree dataset has 243 cells)
    "shared_min": 2,           # an informative allele is carried by >= 2 cells ...
    "ubiquity_max": 0.95,      # ... and by < 95% of cells detected at that site
    "ks": [10, 15],            # fixed selections
    "knee_min": 3,             # smallest knee allowed
    "out": Path(__file__).resolve().parent / "data",
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def site_stats(states, p):
    """states: Series over all cells; NaN = missing, 'U' = uncut, other = allele (str). Returns detected, edited, informative masks."""
    det = states.notna()
    ed = det & (states != "U") & (states != "AMB")
    vc = states[ed].value_counts()
    good = vc[(vc >= p["shared_min"]) & (vc < p["ubiquity_max"] * det.sum())].index
    return det, ed, states.isin(good)

def intbc_rows(dataset, source, group, cells, sites, p):
    """sites: {intBC: [Series per cut site indexed by cells]}"""
    rows = []
    for g, ss in sites.items():
        d = e = inf = pd.Series(False, index=cells); miss = []
        for st in ss:
            dd, ee, ii = site_stats(st, p); d |= dd; e |= ee; inf |= ii; miss.append(1 - dd.mean())
        rows.append(dict(dataset=dataset, source=source, group=group, intBC=g, cells=len(cells), detected=d.mean(),
                         edited=e.mean(), informative=inf.mean(), missing_sites=float(np.mean(miss)), n_sites=len(ss)))
    return rows

def knee(y, kmin):
    """Inflection point of the decreasing sorted informative-fraction curve: the point farthest (either side) from the
    chord joining the first and last point - works for plateau-then-drop (concave) and steep-then-flat (convex) curves."""
    y = np.asarray(y, float); n = len(y)
    if n <= kmin: return n
    x = np.linspace(0, 1, n); yn = (y - y.min()) / (np.ptp(y) or 1)
    chord = yn[0] + (yn[-1] - yn[0]) * x
    return int(max(kmin, np.argmax(np.abs(chord - yn)) + 1))

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}; rows = []
    mjz019 = {"JZ201", "JZ202", "JZ203", "JZ204", "JZ205"}
    for d in p["ours"]:
        proj = "Npp53" if "Npp53_tracing" in str(d) else "Np53Rb1"
        for f in sorted(glob.glob(f"{d}/*_display_matrix.csv")):
            s = os.path.basename(f).split("_display")[0]
            cm = pd.read_csv(f, index_col=0, dtype=str); keys = json.load(open(f.replace("_display_matrix.csv", "_columns.json")))
            conv = cm.apply(lambda c: c.map(lambda v: np.nan if v == "-1" else ("U" if v == "0" else ("AMB" if v.startswith("(") else v))))
            sites = {}
            for j, k in enumerate(keys): sites.setdefault(k[:-2], []).append(conv.iloc[:, j])
            rows += intbc_rows(s, proj, "Ours · MJZ019" if s in mjz019 else "Ours · other libraries", cm.index, sites, p)
    a = pd.read_csv(p["kp_table"], sep="\t", index_col=0, dtype=str, keep_default_na=False)
    for t, g in a.groupby("Tumor"):
        cells = pd.Index(sorted(g.cellBC.unique()))
        if len(cells) < p["min_cells"]: continue
        g = g.drop_duplicates(["cellBC", "intBC"]).set_index("cellBC")
        sites = {}
        for ib, h in g.groupby("intBC"):
            ss = []
            for c in ["r1", "r2", "r3"]:
                v = h[c].replace("", np.nan).map(lambda x: "U" if isinstance(x, str) and "[None]" in x else x)
                ss.append(v.reindex(cells))
            sites[ib] = ss
        rows += intbc_rows(t, "KP-Tracer", "KP-Tracer (Yang 2022)", cells, sites, p)
    I = pd.DataFrame(rows); p["out"].mkdir(exist_ok=True); I.to_csv(p["out"] / "intbc_info_both.csv", index=False)
    S = []
    for ds, g in I.groupby("dataset"):
        g = g.sort_values("informative", ascending=False)
        kk = knee(g.informative.values, p["knee_min"])
        r = dict(dataset=ds, source=g.source.iloc[0], group=g.group.iloc[0], cells=int(g.cells.iloc[0]), n_intBCs=len(g), knee_k=kk,
                 missing_all=g.missing_sites.mean())
        for k in p["ks"] + ["knee"]:
            n = kk if k == "knee" else min(k, len(g))
            r[f"missing_best{k}"] = g.missing_sites.iloc[:n].mean(); r[f"informative_best{k}"] = g.informative.iloc[:n].mean()
        S.append(r)
    S = pd.DataFrame(S); S.to_csv(p["out"] / "selection_missing_both.csv", index=False)
    return {"output_path": str(p["out"]), "params_used": {k: str(v) for k, v in p.items()},
            "summary": S.groupby("group")[["n_intBCs", "knee_k", "missing_all", "missing_best10", "missing_best15", "missing_bestknee"]].median().round(3).to_dict()}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    r = run_pipeline(); print(pd.DataFrame(r["summary"]).to_string())

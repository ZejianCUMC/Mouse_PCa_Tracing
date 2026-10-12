# =============================================================================
# Per-sample barcode heatmap v3 Pipeline
# =============================================================================
# Description : v2 (gen_barcode_heatmaps_v2.py) + Zejian 2026-10-04: show only founder intBCs plus
#               genuine non-founder intBCs; DROP error intBCs = an intBC within Levenshtein <= 2 of a
#               more-detected intBC in the same mouse whose detected-cell ratio (child/parent) <= 0.30
#               (all 108 error intBCs found by scripts/intbc_error_scan_20261004 have ratio <= 0.298;
#               genuine non-founder JZ202 intBCs are >= 6 edits from any founder). Then the v2 10 %
#               coverage display rule. Colours unchanged (white = uncut, grey = missing).
# Input       : v2 inputs, with JZ135_n3000/JZ136/JZ137/JZ138 from the error-intBC rebuild; own _display_matrices (stage --matrices)
# Output      : <project>/figures/barcode_heatmaps_v3_<date>/*_barcode_heatmap.{png,pdf}, JOINS.csv,
#               ERROR_INTBCS_DROPPED.csv
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python
# Key deps    : ete3, matplotlib, pandas
# Validated   : 2026-10-04
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_barcode_heatmaps_v2 as v2

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
DEFAULTS = {
    **v2.DEFAULTS,
    "error_max_lev": 2,          # error intBC: <= 2 edits from a more-detected intBC in the same mouse
    "error_max_ratio": 0.30,     # ... and detected in <= 30 % as many cells as that intBC
    # rebuilt after error-intBC removal (2026-10-04): JZ137/JZ138/JZ135_n3000 trees+PhyloVelo, JZ136 tree/input/PhyloVelo
    "rb1_runs": [(s, v2.v1.RL / "intbcfix_20261004_seed42/cassiopeia", v2.v1.RL / "intbcfix_20261004_seed42")
                 if s in ("JZ137", "JZ138", "JZ135_n3000") else (s, c, t_) for s, c, t_ in v2.DEFAULTS["rb1_runs"]],
    "npp53_tree_override": {"JZ136": v2.v1.N / "lineage_tracing/JZ136_intbcfix_20261004/greedy/JZ136_tree_greedy.nwk"},
    "npp53_allele_override": {"JZ136": v2.v1.N / "lineage_tracing/JZ136_intbcfix_20261004/input/JZ136_allele_table.csv"},
    "npp53_pt_override": {"JZ136": v2.v1.N / "lineage_tracing/JZ136_intbcfix_20261004/phylovelo_out/JZ136_phylovelo_pseudotime.csv"},
    "npp53_out": v2.v1.N / "figures/barcode_heatmaps_v3_20261004",
    "rb1_out": v2.v1.R / "figures/barcode_heatmaps_v3_20261004",
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def lev(a, b):
    d = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        prev, d[0] = d[0], i
        for j, y in enumerate(b, 1):
            prev, d[j] = d[j], min(d[j] + 1, d[j - 1] + 1, prev + (x != y))
    return d[-1]

def error_intbcs(cm, keys, p):
    """{error intBC: (parent, distance, child cells, parent cells)} for one mouse."""
    grp = {}
    for j, k in enumerate(keys): grp.setdefault(k[:-2], []).append(j)
    obs = (cm >= 0) | cm.isna()   # detected = any observed state; NaN = ambiguous multi-allele call (observed), -1 = missing
    det = {g: int(obs.iloc[:, js].any(axis=1).sum()) for g, js in grp.items()}
    out = {}
    for g in grp:
        cand = [(h, lev(g, h)) for h in grp if det[h] > det[g]]
        cand = [(h, d) for h, d in cand if d <= p["error_max_lev"] and det[g] <= p["error_max_ratio"] * det[h]]
        if cand:
            h, d = min(cand, key=lambda x: (x[1], -det[x[0]]))
            out[g] = (h, d, det[g], det[h])
    return out

def patch(p, log):
    _plot = v2.plot_panel
    def plot_panel(tr, states, cm, keys, pt, pal, title, out_stem, pp, note=None):
        cm = cm.loc[[c for c in cm.index]]   # leaf rows (already restricted)
        err = error_intbcs(cm, keys, p)
        keep = [j for j, k in enumerate(keys) if k[:-2] not in err]
        for g, (h, d, n, m) in err.items():
            log.append(dict(sample=Path(out_stem).name, error_intBC=g, parent_intBC=h, distance=d, cells=n, parent_cells=m))
        if err:
            enote = f"{len(err)} error intBCs (1-2 base variants of a more-detected intBC) removed"
            note = f"{note} · {enote}" if note else enote
        return _plot(tr, states, cm.iloc[:, keep], [keys[j] for j in keep], pt, pal, title, out_stem, pp, note)
    v2.plot_panel = plot_panel

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    for k in ["npp53_out", "rb1_out"]: Path(p[k]).mkdir(parents=True, exist_ok=True)   # panels themselves refuse overwrite (v1 plot_panel)
    log = []; patch(p, log)
    res = v2.run_pipeline(**p)
    d = pd.DataFrame(log, columns=["sample", "error_intBC", "parent_intBC", "distance", "cells", "parent_cells"])
    for proj, k in [("Npp53", "npp53_out"), ("Np53Rb1", "rb1_out")]:
        d[d["sample"].str.startswith(proj + "_")].to_csv(Path(p[k]) / "ERROR_INTBCS_DROPPED.csv", index=False)
    return res

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--only", nargs="*"); ap.add_argument("--out_root"); ap.add_argument("--matrices", action="store_true")
    a = ap.parse_args(); kw = {}
    if a.only:
        kw["npp53_samples"] = [s for s in DEFAULTS["npp53_samples"] if s in a.only]
        kw["rb1_runs"] = [r for r in DEFAULTS["rb1_runs"] if r[0] in a.only]
    if a.out_root:
        kw["npp53_out"] = Path(a.out_root) / "Npp53"; kw["rb1_out"] = Path(a.out_root) / "Np53Rb1"
    ap2 = sys.argv
    if "--matrices" in ap2: v2.build_matrices({**DEFAULTS, **kw})   # stage 1 (Cassiopeia_v2 env)
    else: run_pipeline(**kw)

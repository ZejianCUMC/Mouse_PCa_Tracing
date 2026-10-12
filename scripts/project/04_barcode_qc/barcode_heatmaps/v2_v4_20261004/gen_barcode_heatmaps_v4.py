# =============================================================================
# Per-sample barcode heatmap v4 Pipeline
# =============================================================================
# Description : v3 (gen_barcode_heatmaps_v3.py) with Zejian 2026-10-05 changes:
#               1) colour key flipped: uncut (un-edited) = LIGHT GREY, missing (not detected / mis-sequenced) = WHITE;
#               2) JZ201-205 from the MJZ019 re-run (min_reads_per_umi=1): HPC-rebuilt Greedy trees, re-run inputs
#                  (error intBCs removed) and the new PhyloVelo (UMAP) pseudotime.
#               Display rules unchanged (error intBCs dropped, intBCs detected in >=10% of cells, all 3 cut sites).
# Input       : v3 inputs + Npp53_tracing/lineage_tracing/MJZ019_rerun_trees_20261005/{inputs,greedy,phylovelo_out}
# Output      : <project>/figures/barcode_heatmaps_v4_20261005/ (+ own _display_matrices, JOINS.csv, ERROR_INTBCS_DROPPED.csv)
# Conda env   : stage --matrices: /home/jrchristin/miniconda3/envs/Cassiopeia_v2/bin/python ; plot: hotspot env
# Validated   : 2026-10-05 (JZ201 -> JZ201_excl_outlier_20261006 on 2026-10-06)
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_barcode_heatmaps_v3 as v3

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
RR = v3.v2.v1.N / "lineage_tracing/MJZ019_rerun_trees_20261005"
M019 = ["JZ201", "JZ202", "JZ203", "JZ204", "JZ205"]
X201 = v3.v2.v1.N / "lineage_tracing/JZ201_excl_outlier_20261006"   # 2026-10-06: JZ201 re-run without isolated UMAP cell CACTGTCCACATGCAG
_rr = lambda s: X201 if s == "JZ201" else RR
DEFAULTS = {
    **v3.DEFAULTS,
    "uncut_color": "#D3D3D3",     # un-edited site = light grey (was white)
    "missing_color": "#FFFFFF",   # missing / mis-sequenced = white (was grey #BEBEBE)
    "npp53_tree_override": {**v3.DEFAULTS["npp53_tree_override"], **{s: _rr(s) / "greedy" / f"{s}_tree_greedy.nwk" for s in M019}},
    "npp53_allele_override": {**v3.DEFAULTS["npp53_allele_override"], **{s: _rr(s) / "inputs" / f"{s}_allele_table.csv" for s in M019}},
    "npp53_pt_override": {**v3.DEFAULTS["npp53_pt_override"], **{s: _rr(s) / "phylovelo_out" / f"{s}_phylovelo_pseudotime.csv" for s in M019}},
    "npp53_out": v3.v2.v1.N / "figures/barcode_heatmaps_v4_20261005",
    "rb1_out": v3.v2.v1.R / "figures/barcode_heatmaps_v4_20261005",
}
DEFAULTS.pop("display_matrix_dirs", None)   # build own display matrices (JZ201-205 inputs changed)
# 2026-10-06 (Zejian): JZ272 adopted from the state>=20-cells rebuild (isolated Mes-Chemokine cell stretched PhyloVelo pseudotime)
_X272 = v3.v2.v1.R / "lineage_tracing/JZ272_state_min20_test_20261006"
DEFAULTS["rb1_runs"] = [("JZ272", _X272 / "cassiopeia", _X272 / "viewer") if r[0] == "JZ272" else r for r in DEFAULTS["rb1_runs"]]
# 2026-10-06 (Zejian): same rule adopted for JZ267, JZ269, JZ274
_GM = v3.v2.v1.R / "lineage_tracing/state_min20_rebuild_20261006"
DEFAULTS["rb1_runs"] = [(r[0], _GM / r[0] / "cassiopeia", _GM / r[0] / "viewer") if r[0] in ("JZ267", "JZ269", "JZ274") else r for r in DEFAULTS["rb1_runs"]]

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def set_colours(p):
    """cm_to_rgb in gen_barcode_joint reads module-level UNCUT / MISSING at call time."""
    sys.path.insert(0, p["ref_dir"])
    import gen_barcode_joint as gbj
    gbj.UNCUT, gbj.MISSING = p["uncut_color"], p["missing_color"]

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    set_colours(p)
    return v3.run_pipeline(**p)

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--matrices", action="store_true"); ap.add_argument("--only", nargs="*"); ap.add_argument("--out_root")
    a = ap.parse_args(); kw = {}
    if a.only:
        kw["npp53_samples"] = [s for s in DEFAULTS["npp53_samples"] if s in a.only]
        kw["rb1_runs"] = [r for r in DEFAULTS["rb1_runs"] if r[0] in a.only]
    if a.out_root:
        kw["npp53_out"] = Path(a.out_root) / "Npp53"; kw["rb1_out"] = Path(a.out_root) / "Np53Rb1"
    if a.matrices: v3.v2.build_matrices({**DEFAULTS, **kw})
    else: run_pipeline(**kw)

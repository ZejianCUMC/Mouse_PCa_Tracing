# Barcode heatmaps (v3 rules, unchanged) for JZ201-205 on the MJZ019 re-run trees. Pseudotime: existing per-cell PhyloVelo
# values (v8 run) - a new PhyloVelo needs HPC. Stage 'matrices' in Cassiopeia_v2 env, 'plot' in hotspot env.
import sys
from pathlib import Path
sys.path.insert(0, "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/scripts/barcode_heatmaps_20261004")
import gen_barcode_heatmaps_v3 as v3
OUT = Path(__file__).resolve().parents[1]; S = ["JZ201", "JZ202", "JZ203", "JZ204", "JZ205"]
p = {**v3.DEFAULTS, "npp53_samples": S, "rb1_runs": [],
     "npp53_tree_override": {s: OUT / "greedy" / f"{s}_tree_greedy.nwk" for s in S},
     "npp53_allele_override": {s: OUT / "inputs" / f"{s}_allele_table.csv" for s in S},
     "npp53_pt_override": {},
     "npp53_out": OUT / "heatmaps", "rb1_out": OUT / "heatmaps_rb1_unused"}
if sys.argv[1] == "matrices": v3.v2.build_matrices(p)
else: v3.run_pipeline(**p)

# Report "before" panels redrawn with the v4 colour key (uncut light grey, missing white):
#  - original MJZ019 JZ201-205 (v8 trees/inputs, v3 rules)  - JZ136 with error intBCs (v2 rules, old tree)
import sys
from pathlib import Path
sys.path.insert(0, "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/scripts/barcode_heatmaps_20261004")
import gen_barcode_heatmaps_v4 as v4, gen_barcode_heatmaps_v3 as v3
O = Path(__file__).resolve().parent / "before_panels_v4colours"
v4.set_colours(v4.DEFAULTS)
if sys.argv[1] == "mjz019":     # original JZ201-205: v3 defaults (no MJZ019 overrides) + v4 colours, v3 display matrices
    v3.run_pipeline(npp53_samples=["JZ201", "JZ202", "JZ203", "JZ204", "JZ205"], rb1_runs=[], npp53_out=O / "MJZ019_original", rb1_out=O / "_unused1",
                    display_matrix_dirs={"Npp53": v3.v2.v1.N / "figures/barcode_heatmaps_v3_20261004/_display_matrices"})
else:                           # JZ136 with error intBCs: v2 rules on old tree, v2 display matrices
    v3.v2.run_pipeline(npp53_samples=["JZ136"], rb1_runs=[], npp53_out=O / "JZ136_v2_before", rb1_out=O / "_unused2",
                       display_matrix_dirs={"Npp53": v3.v2.v1.N / "figures/barcode_heatmaps_v2_20261004/_display_matrices"})

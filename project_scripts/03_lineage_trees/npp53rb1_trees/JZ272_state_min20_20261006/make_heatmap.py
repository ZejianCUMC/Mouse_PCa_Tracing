# Heatmap v4 (same rules/colours) for the JZ272 state>=20 test tree.  Stage 1: --matrices (Cassiopeia_v2 env); stage 2: plot (hotspot env).
import sys
from pathlib import Path
sys.path.insert(0, "/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/scripts/barcode_heatmaps_20261004")
import gen_barcode_heatmaps_v4 as v4
R = Path(__file__).resolve().parent; out = R / "heatmap"
kw = dict(npp53_samples=[], rb1_runs=[("JZ272", R / "cassiopeia", R / "viewer")], rb1_out=out, npp53_out=out / "_unused")
if "--matrices" in sys.argv: out.mkdir(exist_ok=True); (out / "_unused").mkdir(exist_ok=True); v4.v3.v2.build_matrices({**v4.DEFAULTS, **kw})
else: v4.run_pipeline(**kw)

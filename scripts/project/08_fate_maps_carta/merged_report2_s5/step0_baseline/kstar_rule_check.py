#!/usr/bin/env python
# Step 0b: which k-selection rule gives the k used by the published per-tumor maps (carta_extracted.json)?
# Chord-elbow (compute_best_k_v2.py logic: plateau filter, axes scaled to [0,1], max |x+y-1|) on k = 1..KMAX.
import json
from pathlib import Path
import numpy as np, pandas as pd
CF = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12")
TR = Path("/media/jrchristin/Backup10TB1/agent/Trash/JZ201_previous_runs_20261006/Npp53_tracing/analysis/evolution_features_v2/carta_fine12")
ext = json.loads((CF / "biology_review/carta_extracted.json").read_text())
summ = pd.read_csv(CF / "fine12_greedy_summary.tsv", sep="\t").set_index("sample")
def curve(s, kmax):
    d = CF / "results_fine12_greedy" / s
    if not d.exists(): d = TR / "results_fine12_greedy" / s
    rows = []
    for k in range(1, kmax + 1):
        f = d / f"{s}_results_{k}.txt"
        if f.exists():
            r = pd.read_csv(f, sep="\t").iloc[0]; rows.append((int(r["k"]), float(r["objective_score"])))
    c = pd.DataFrame(rows, columns=["k", "obj"]).drop_duplicates("k").sort_values("k").reset_index(drop=True)
    c["plateau"] = False
    for i in range(1, len(c)):
        if c.obj[i] >= c.obj[i - 1] - 1e-9: c.loc[i:, "plateau"] = True; break
    return c
def chord(c):
    u = c[~c.plateau]; xs, ys = u.k.values.astype(float), (u.obj / 1.0).values
    if len(xs) < 3: return int(xs[0])
    xn = (xs - xs.min()) / max(1e-12, xs.max() - xs.min()); yn = (ys - ys.min()) / max(1e-12, ys.max() - ys.min())
    return int(xs[np.argmax(np.abs(xn + yn - 1) / np.sqrt(2))])
rows = []
for s in sorted(ext):
    r = {"sample": s, "k_used_in_maps": ext[s]["k"], "n_states": int(summ.loc[s, "n_states"])}
    for kmax in (5, 6, 7, 8):
        r[f"chord_k1..{kmax}"] = chord(curve(s, kmax))
    rows.append(r)
t = pd.DataFrame(rows)
print(t.to_string(index=False))
for c in t.columns[3:]:
    print(c, "agrees with the maps' k in", int((t[c] == t.k_used_in_maps).sum()), "/ 25")
t.to_csv(Path(__file__).resolve().parent / "kstar_rule_check.tsv", sep="\t", index=False)

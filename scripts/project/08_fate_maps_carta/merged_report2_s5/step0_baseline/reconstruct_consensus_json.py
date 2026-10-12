#!/usr/bin/env python
# =============================================================================
# Step 0a: reconstruct the per-arm pooling behind Part 2 Fig 10.3.5/10.3.6
# =============================================================================
# Description : The script that wrote consensus_fate_per_arm.json was never saved.
#               For every tumor and every k (1..8) parse the CARTA progenitor sets,
#               derive terminal / direct-from-root / co-occurrence statistics under
#               several candidate rules, and search for the per-tumor k (and rule)
#               that reproduces the saved JSON exactly.
# Input       : Npp53 carta_fine12/results_fine12_greedy/*/_progenitors_k.txt,
#               inputs_fine12/*_labels_fine12.tsv, datadriven_fatemap/consensus_fate_per_arm.json,
#               k tables (best_k_per_sample*.csv, biology_review/carta_extracted.json)
# Output      : step0_baseline/reconstruct_consensus_json.{log,json}
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/scvi/bin/python
# Validated   : 2026-10-11
# =============================================================================
import ast, itertools, json, sys
from pathlib import Path
import pandas as pd

CF = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12")
RES = CF / "results_fine12_greedy"; LAB = CF / "inputs_fine12"
# JZ201's original CARTA inputs/results were moved aside on 2026-10-06 (JZ201 tree rebuilt); read them there, read-only
TRASH_CF = Path("/media/jrchristin/Backup10TB1/agent/Trash/JZ201_previous_runs_20261006/Npp53_tracing/analysis/evolution_features_v2/carta_fine12")
def res_dir(s): return RES / s if (RES / s).exists() else TRASH_CF / "results_fine12_greedy" / s
def lab_file(s): return LAB / f"{s}_labels_fine12.tsv" if (LAB / f"{s}_labels_fine12.tsv").exists() else TRASH_CF / "inputs_fine12" / f"{s}_labels_fine12.tsv"
JS = CF / "biology_review/datadriven_fatemap/consensus_fate_per_arm.json"
OUT = Path(__file__).resolve().parent
ARM = {"Intact": "Intact", "Castration": "Cas", "Castration + Enzalutamide": "Cas+Enz", "Cas": "Cas", "Cas+Enz": "Cas+Enz"}

J = json.loads(JS.read_text())
ren = lambda s: "Luminal-ARdep" if s == "Luminal-ARhi" else s
STATES = [ren(s) for s in J["states"]]
def tvec(d): return tuple(d.get(s, 0) for s in STATES)
TGT = {a: dict(term=tuple(J["in_terminal_by_arm"][a][s] for s in J["states"]),
               direct=tuple(J["direct_from_root_by_arm"][a][s] for s in J["states"]),
               cooc={tuple(ren(x) for x in k.split("|")): v for k, v in J["cooccur_by_arm"][a].items()})
       for a in J["n_by_arm"]}

summ = pd.read_csv(CF / "fine12_greedy_summary.tsv", sep="\t"); arm_of = {r.sample: ARM[r.arm] for r in summ.itertuples()}
present = {s: set(pd.read_csv(lab_file(s), sep="\t")["cell_state"]) for s in arm_of}

def parse(p):
    out = []
    for ln in p.read_text().splitlines():
        ln = ln.strip()
        if ln: out.append(frozenset(ren(x) for x in ast.literal_eval(ln)))
    return out
def stats(sets, keep, collapse=False):
    if collapse:   # carta_extracted.json-style DAG: restrict sets to present states, drop empties and repeats of an earlier set
        seen, fs = [], []
        for a in sets:
            b = frozenset(a & keep)
            if b and b not in seen: seen.append(b); fs.append(b)
        sets, keep = fs, set().union(*fs)
        direct_raw_override = None
    root, rest = sets[0], sets[1:]
    parents = set()
    for i, a in enumerate(rest):           # parent = smallest proper superset among earlier sets (root default)
        cands = [b for b in [root] + rest[:i] if b > a]
        if cands: parents.add(min(cands, key=len))
    leaves = [a for a in rest if a not in parents]
    leaves_f = [a & keep for a in leaves]
    term = set().union(*leaves_f) if leaves_f else set()
    inner = set().union(*rest) if rest else set()
    direct_raw = set(root) - inner
    pairs = set()
    for a in leaves_f:
        pairs |= {tuple(sorted(p, key=STATES.index)) for p in itertools.combinations(a, 2)}
    return term, direct_raw, pairs

ALL = set(STATES)
rules = {"leaves_raw": lambda s: ALL, "leaves_present": lambda s: present[s], "collapsed_present": lambda s: present[s],
         "smallest_set": lambda s: ALL}   # rule found in the 2026-09-07 session transcript (terminal = the single smallest raw set)
def smallest_stats(sets):
    srt = sorted(sets, key=len, reverse=True)
    root, inter = srt[0], srt[1:]
    term = min(srt, key=len)
    direct = {x for x in root if not any(x in q for q in inter)}
    pairs = {tuple(sorted(pq, key=STATES.index)) for pq in itertools.combinations(term, 2)}
    return set(term), direct, pairs
log = []
def P(*x):
    s = " ".join(str(y) for y in x); print(s); log.append(s)

result = {}
for rname, keepf in rules.items():
    P(f"\n=== rule: terminal = union of leaf potency sets, states kept = {rname}")
    sig = {}
    for s in arm_of:
        sig[s] = {}
        for k in range(1, 9):
            p = res_dir(s) / f"{s}_progenitors_{k}.txt"
            if not p.exists(): continue
            raw = parse(p)
            if rname == "smallest_set":
                t, d, pr = smallest_stats(raw)
            elif rname == "collapsed_present":
                t, _, pr = stats(raw, keepf(s), collapse=True)
                _, d, _ = stats(raw, ALL)            # direct-from-root stays on the raw sets (absent states count)
            else:
                t, d, pr = stats(raw, keepf(s))
            sig[s][k] = (tvec({x: 1 for x in t}), tvec({x: 1 for x in d}), frozenset(pr))
    sols_rule = {}
    for a, tg in TGT.items():
        samples = sorted([s for s in arm_of if arm_of[s] == a])
        tt, td = tg["term"], tg["direct"]
        sols = []
        def dfs(i, acc_t, acc_d, chosen):
            if any(x > y for x, y in zip(acc_t, tt)) or any(x > y for x, y in zip(acc_d, td)): return
            if i == len(samples):
                if acc_t == tt and acc_d == td: sols.append(dict(chosen))
                return
            s = samples[i]; seen = {}
            for k, (t, d, pr) in sig[s].items(): seen.setdefault((t, d, pr), []).append(k)
            for (t, d, pr), ks in seen.items():
                chosen[s] = ks
                dfs(i + 1, tuple(x + y for x, y in zip(acc_t, t)), tuple(x + y for x, y in zip(acc_d, d)), chosen)
            chosen.pop(s, None)
        dfs(0, (0,) * len(STATES), (0,) * len(STATES), {})
        good = []
        for so in sols:                   # check co-occurrence counts too (one count per tumor per pair)
            cnt = {}
            for s, ks in so.items():
                for pr in sig[s][ks[0]][2]: cnt[pr] = cnt.get(pr, 0) + 1
            if cnt == {k: v for k, v in tg["cooc"].items() if v}: good.append(so)
        P(f"  {a}: {len(sols)} k-assignments match terminal+direct counts; {len(good)} also match co-occurrence")
        for so in good[:5]: P("   ", {s: ks for s, ks in so.items()})
        sols_rule[a] = good
    result[rname] = {a: v for a, v in sols_rule.items()}

# compare with the saved k tables
ext = json.loads((CF / "biology_review/carta_extracted.json").read_text())
chord = pd.read_csv(CF / "best_k_per_sample_v2_chord.csv").set_index("sample")
best = pd.read_csv(CF / "best_k_per_sample.csv").set_index("sample")
tables = {"carta_extracted.json k": {s: ext[s]["k"] for s in ext},
          "v2_chord k_star_final": chord["k_star_final"].astype(int).to_dict(),
          "v2_chord k_star_chord": chord["k_star_chord"].astype(int).to_dict(),
          "best_k k_consensus": best["k_consensus"].astype(int).to_dict(),
          "summary elbow_k": summ.set_index("sample")["elbow_k"].astype(int).to_dict()}
P("\n=== which saved k table is consistent with the matching assignments")
for rname, per_arm in result.items():
    for tname, tab in tables.items():
        ok = all(any(all(tab[s] in ks for s, ks in so.items()) for so in sols) for a, sols in per_arm.items() if sols) and all(per_arm.values())
        bad = {}
        for a, sols in per_arm.items():
            if not sols: continue
            # samples whose table k is outside every matching assignment
            bad[a] = sorted(s for s in sols[0] if not any(tab[s] in so[s] for so in sols))
        P(f"  rule={rname:15s} table={tname:24s} reproduces JSON: {ok}   tumors off: {bad}")
(OUT / "reconstruct_consensus_json.log").write_text("\n".join(log) + "\n")
(OUT / "reconstruct_consensus_json.json").write_text(json.dumps({r: {a: s for a, s in v.items()} for r, v in result.items()}, indent=1, default=list))

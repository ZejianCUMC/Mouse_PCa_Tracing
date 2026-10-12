"""Final modules from the Hotspot modules (Zejian 2026-10-11: "start at hotspot modules result (10/11), then merge modules with closed
biological functions or <80 genes"; "let risk module picked as one of final module, do not split except sub-modules pass survival test").
Rule, same for both genotypes:
  0. The risk module (its two neighbouring Hotspot modules) is one fixed final module: neither of its parts passes the survival rule on
     its own (median split, HR > 1 and p < 0.05 in both SU2C and WCDT-JCO), so it is not split, and nothing is merged into it.
  1. Function: two Hotspot modules are merged when they share at least one enriched term (GO BP or MSigDB Hallmark, BH FDR < 0.05,
     universe = genes in that genotype's meta-analysis; 02_go_enrich.R). Merges are transitive.
  2. Size: while a final module has fewer than MIN_GENES genes, the smallest one is merged into the final module it is most co-expressed
     with (highest mean Hotspot local-correlation Z over all gene pairs between the two).
Read-only on Hotspot outputs. Env: moscot."""
import json, itertools, numpy as np, pandas as pd
from pathlib import Path
O = Path(__file__).resolve().parent; B = O.parents[2]
MIN_GENES = 80                                       # Zejian's size rule
RISK = {"Npp53": (7, 6), "Npp53Rb1": (6, 8)}
GO = pd.read_csv(O / "go_bp_per_module.csv"); HM = pd.read_csv(O / "hallmark_per_module.csv")
rows, steps = [], []
for g, risk in RISK.items():
    A = B / f"metaNorm_20261010_{g}_bcr" / "hotspot_out"; lid = f"META_{g}_FSU_top1000_n3000"
    S = json.load(open(A / f"{lid}_k30_t1000.json"))["sets"]["hs_top1000_min30"]; order = S["order"]
    Z = np.load(A / f"{lid}_k30_t1000_lcz.npy"); zg = (A / f"{lid}_k30_t1000_lcz_genes.txt").read_text().split(); pos = {x: i for i, x in enumerate(zg)}
    terms = {m: set(GO[(GO.genotype == g) & (GO.hotspot_module == m)].ID) | set(HM[(HM.genotype == g) & (HM.hotspot_module == m)].ID) for m in order}
    free = [m for m in order if m not in risk]
    parent = {m: m for m in free}
    def find(m):
        while parent[m] != m: m = parent[m]
        return m
    for a, b in itertools.combinations(free, 2):            # step 1: shared enriched term
        sh = terms[a] & terms[b]
        if sh:
            steps.append(dict(genotype=g, step="function", merged=f"M{a} + M{b}", reason="shared: " + "; ".join(sorted(sh)), mean_Z=np.nan))
            parent[find(b)] = find(a)
    grp = {}
    for m in free: grp.setdefault(find(m), []).append(m)
    units = {tuple(sorted(v, key=order.index)): sum(([pos[x] for x in S["modules"][str(m)]] for m in v), []) for v in grp.values()}
    rk = tuple(sorted(risk, key=order.index)); units[rk] = sum(([pos[x] for x in S["modules"][str(m)]] for m in risk), [])
    bz = lambda a, b: float(np.nanmean(Z[np.ix_(units[a], units[b])]))
    while True:                                              # step 2: size
        small = [u for u in units if u != rk and len(units[u]) < MIN_GENES]
        if not small: break
        u = min(small, key=lambda v: (len(units[v]), order.index(v[0])))
        cand = sorted(((bz(u, v), v) for v in units if v not in (u, rk)), reverse=True); s, v = cand[0]
        steps.append(dict(genotype=g, step="size", merged=f"{'+'.join('M%d' % m for m in u)} ({len(units[u])} genes) into {'+'.join('M%d' % m for m in v)}",
                          reason=f"< {MIN_GENES} genes; most co-expressed final module (vs risk module Z {bz(u, rk):.2f})", mean_Z=round(s, 2)))
        w = tuple(sorted(u + v, key=order.index)); units[w] = units.pop(u) + units.pop(v)
    fin = sorted(units, key=lambda u: min(order.index(m) for m in u)); k = 0
    for u in fin:
        if u == rk: name = "risk module"
        else: k += 1; name = f"module {k}"
        genes = sorted(zg[i] for i in units[u]); W = Z[np.ix_(units[u], units[u])]
        rows.append(dict(genotype=g, final_module=name, n_genes=len(genes), hotspot_modules=" ".join(f"M{m}" for m in u), n_hotspot_modules=len(u),
                         within_mean_Z=round(float(np.nanmean(W[~np.eye(len(W), dtype=bool)])), 2), genes=" ".join(genes)))
F = pd.DataFrame(rows); F.to_csv(O / "final_modules.csv", index=False); pd.DataFrame(steps).to_csv(O / "final_module_steps.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print(F.drop(columns="genes").to_string()); print(pd.DataFrame(steps).to_string())

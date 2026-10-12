"""opt5 step 06 — candidate gene sets with their fitness direction (usage: 06_candidates.py <arm>). As opt4 06:
  Hotspot (<arm>/hotspot_out/*_k30.json): each module + each dendrogram-adjacent pair (hs_* sets), merged-5 modules, ward-3/5 clusters;
  fitness-profile clusters (<arm>/profile_out/*.json): each cluster + each adjacent pair; whole lists (kept Npp53 / Npp53Rb1 n3000 lists);
  bidirectional whole signatures (bidir_lists.json: up = FSU list, down = FSD list).
Identical mouse sets of the same direction merged. Human orthologs via orthologs.py (audit fixes: ALL human orthologs of a homology
class; outdated 10x mouse symbols rescued through HGNC); kept if 5-400 human genes (BI: each side 5-1500).
Output <arm>/candidates_human.json, candidates_meta.csv, ortholog_map_list_genes.csv, pool_human.txt (human orthologs of every gene
in a primary list = random-set pool of step 09). Conda env: moscot."""
import sys, json, glob, pandas as pd
from pathlib import Path
O = Path(__file__).resolve().parent; sys.path.insert(0, str(O)); import orthologs as ORT
ARM = sys.argv[1]; A = O / ARM
L = json.load(open(A / "fsu_lists.json")); C = {}
def add(prov, genes): C.setdefault((L[prov.split("|")[0]]["direction"], frozenset(genes)), []).append(prov)
for src, pat in [("hs", "hotspot_out/*.json"), ("prof", "profile_out/*.json")]:
    for f in sorted(glob.glob(str(A / pat))):
        r = json.load(open(f)); base = f"{r['list_id']}|" + (f"k{r['k']}" + ("" if r.get("top", 500) == 500 else f"t{r['top']}") if src == "hs" else src)
        for sname, s in r["sets"].items():
            mods = {int(k): v for k, v in s["modules"].items()}; order = [int(x) for x in s["order"]]
            for m in order: add(f"{base}|{sname}|M{m}", mods[m])
            if sname.startswith(("hs_", "prof")):
                for a, b in zip(order[:-1], order[1:]): add(f"{base}|{sname}|M{a}+M{b}", mods[a] + mods[b])
prim = {lid: d for lid, d in L.items() if d["cells"] in ("Npp53", "Npp53Rb1", "all")}
for lid, d in prim.items(): add(f"{lid}|whole|list|all", d["genes"])
rows, human = [], {}
for i, ((dr, key), prov) in enumerate(sorted(C.items(), key=lambda kv: (kv[0][0], len(kv[0][1]), sorted(kv[0][1])))):
    h = ORT.to_human(key)
    if not 5 <= len(h) <= 400: continue
    cid = f"Q{i:05d}"; human[cid] = h; genos = sorted({L[p.split("|")[0]]["genotype"] for p in prov})
    rows.append(dict(cand=cid, direction=dr, genotype="+".join(genos), n_mouse=len(key), n_human=len(h), n_configs=len(prov), n_lists=len({p.split("|")[0] for p in prov}),
                     first_config=prov[0], configs=";".join(prov), genes=" ".join(sorted(key))))
B = json.load(open(A / "bidir_lists.json")); seen = {}
for j, (bid, d) in enumerate(sorted(B.items())):
    up, dn = ORT.to_human(d["up"]), ORT.to_human(d["dn"]); k = (tuple(up), tuple(dn))
    if not (5 <= len(up) <= 1500 and 5 <= len(dn) <= 1500): continue
    if k in seen: rows[seen[k]]["configs"] += f";{bid}|whole|bidir|all"; rows[seen[k]]["n_configs"] += 1; continue
    cid = f"B{j:04d}"; human[cid] = dict(up=up, dn=dn); seen[k] = len(rows)
    rows.append(dict(cand=cid, direction="BI", genotype=bid.split("_")[1], n_mouse=len(d["up"]) + len(d["dn"]), n_human=len(up) + len(dn), n_configs=1, n_lists=1,
                     first_config=f"{bid}|whole|bidir|all", configs=f"{bid}|whole|bidir|all", genes="UP: " + " ".join(d["up"]) + " | DN: " + " ".join(d["dn"])))
json.dump(human, open(A / "candidates_human.json", "w")); M = pd.DataFrame(rows); M.to_csv(A / "candidates_meta.csv", index=False)
allg = sorted({g for d in prim.values() for g in d["genes"]}); R = ORT.rescues(allg); R.to_csv(A / "ortholog_map_list_genes.csv", index=False)
(A / "pool_human.txt").write_text("\n".join(ORT.to_human(allg)) + "\n")
print(f"{ARM}: hotspot outputs {len(glob.glob(str(A / 'hotspot_out/*.json')))} | profile outputs {len(glob.glob(str(A / 'profile_out/*.json')))} | candidates {len(M)}")
print(M.groupby(["genotype", "direction"]).size().to_string()); print("list-gene ortholog mapping:", R.how.value_counts().to_dict())

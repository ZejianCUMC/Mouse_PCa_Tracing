"""opt5 step 05c — fitness-PROFILE modules, no Hotspot (usage: 05c_profile_modules.py <arm>). Same as opt4 05c: for each kept FSU /
FSD list of the primary cell sets (Npp53; Npp53Rb1 with JZ135 tree n3000), the genes (500 with the largest |meta-z| when longer;
genes tested in < half the units rank last) are clustered by the Pearson correlation of their signed per-tumour fitness z across that
genotype's tumours (untested = 0), average linkage on 1 - r, cut into 3 / 5 / 8; clusters >= 8 genes kept. Deterministic.
Output <arm>/profile_out/<list_id>.json. Conda env: moscot."""
import sys, json, numpy as np, pandas as pd
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform
O = Path(__file__).resolve().parent; A = O / sys.argv[1]; OUT = A / "profile_out"; OUT.mkdir(exist_ok=True)
L = json.load(open(A / "fsu_lists.json")); UZ = {}; n = 0
def absz(d, x): v = d["meta_z"].get(x); return 0.0 if v is None else abs(v)
for lid, d in L.items():
    if d["cells"] not in ("Npp53", "Npp53Rb1", "all"): continue
    if d["method"] not in UZ: UZ[d["method"]] = pd.read_csv(A / f"unit_z_{d['method']}.csv.gz", index_col=0)
    Z = UZ[d["method"]].loc[d["units"]]; g = [x for x in sorted(d["genes"], key=lambda x: (-absz(d, x), x)) if x in Z.columns][:500]
    P = Z[g].fillna(0).values.T; R = np.nan_to_num(np.corrcoef(P)); np.fill_diagonal(R, 1)
    Lk = linkage(squareform(1 - R, checks=False), "average"); order = [g[i] for i in leaves_list(Lk)]; sets = {}
    for kc in (3, 5, 8):
        cl = pd.Series(fcluster(Lk, kc, "maxclust"), index=g); mo = [int(c) for c in pd.unique(cl.reindex(order).values) if (cl == c).sum() >= 8]
        sets[f"prof{kc}"] = dict(order=mo, modules={c: sorted(cl.index[cl == c]) for c in mo})
    (OUT / f"{lid}.json").write_text(json.dumps(dict(list_id=lid, n_genes=len(g), sets=sets))); n += 1
print(sys.argv[1], "profile module sets written:", n)

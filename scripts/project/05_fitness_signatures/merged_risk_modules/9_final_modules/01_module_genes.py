"""Export the Hotspot modules (min-30 set, k = 30, top 1,000) and the meta-analysis universe per genotype for GO enrichment. Env: moscot."""
import json, pandas as pd
from pathlib import Path
O = Path(__file__).resolve().parent; B = O.parents[2]
rows, uni = [], []
for g in ["Npp53", "Npp53Rb1"]:
    A = B / f"metaNorm_20261010_{g}_bcr"; S = json.load(open(A / "hotspot_out" / f"META_{g}_FSU_top1000_n3000_k30_t1000.json"))["sets"]["hs_top1000_min30"]
    for i, m in enumerate(S["order"]):
        for x in S["modules"][str(m)]: rows.append(dict(genotype=g, hotspot_module=m, dendro_pos=i, gene=x))
    M = pd.read_csv(A / f"meta_{g}.csv.gz", usecols=["gene"]); uni += [dict(genotype=g, gene=x) for x in M.gene]
pd.DataFrame(rows).to_csv(O / "hotspot_module_genes.csv", index=False); pd.DataFrame(uni).to_csv(O / "universe_genes.csv", index=False)
print(pd.DataFrame(rows).groupby(["genotype", "hotspot_module"]).size().to_string())

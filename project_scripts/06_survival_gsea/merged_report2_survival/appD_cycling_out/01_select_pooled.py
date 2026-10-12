"""Appendix D (cycling cells excluded) (Zejian 2026-10-11: "Appendix C show best candidates risk factor for Npp53-Npp53Rb1 merged"): the same pipeline as the
per-genotype risk modules, run on all 47 tumors pooled (arm opt5_corrected_20261009/metaNorm_20261011_all_cycout_bcr: log1p CP10k, Blom fitness,
CR1 SE, random-effects meta-analysis over 47 tumors, FSU top 1,000, Hotspot k 30 / top 1,000 / min 30 on all 43,082 cells). Candidates =
each Hotspot module and each dendrogram-adjacent pair (5-400 human genes); rule = median split, HR > 1 and log-rank p < 0.05 in SU2C and
WCDT-JCO; expected by chance = sum of the random-set pass rate at the nearest size (09_null.R). Top 5 by the worse selection p.
Writes pooled_candidates.csv (all), pooled_top.csv, pooled_top_human.json, pooled_summary.json. Env: moscot."""
import json, numpy as np, pandas as pd
from pathlib import Path
O = Path(__file__).resolve().parent; B = O.parents[2]; A = B / "metaNorm_20261011_all_cycout_bcr"; SU, JC, TOP = "SU2C-mCRPC-union", "WCDT-JCO", 5
M = pd.read_csv(A / "candidates_meta.csv"); M = M[(M.direction == "FSU") & M.configs.str.contains(r"META_all_FSU_top1000_n3000\|k30t1000\|hs_top1000_min30\|")].set_index("cand")
K = pd.read_csv(A / "km_scan.csv"); K = K[(K.split == "median") & K.cand.isin(M.index) & K.cohort.isin([SU, JC])].pivot_table(index="cand", columns="cohort", values=["HR", "p"])
K.columns = [f"{a}|{b}" for a, b in K.columns]; D = M[["n_mouse", "n_human", "genes"]].join(K)
D["hotspot"] = M.configs.str.extract(r"hs_top1000_min30\|([^;]+)")[0]; D["worst_p"] = D[[f"p|{SU}", f"p|{JC}"]].max(1)
D["meets_rule"] = (D[f"HR|{SU}"] > 1) & (D[f"HR|{JC}"] > 1) & (D.worst_p < .05)
N = pd.read_csv(A / "km_null.csv"); N = N[N.split == "median"].pivot_table(index="cand", columns="cohort", values=["HR", "p"]); N.columns = [f"{a}|{b}" for a, b in N.columns]
ok = (N[f"HR|{SU}"] > 1) & (N[f"p|{SU}"] < .05) & (N[f"HR|{JC}"] > 1) & (N[f"p|{JC}"] < .05)
rate = pd.Series(ok.values, index=N.index.str.extract(r"rand_s(\d+)_")[0].astype(int).values).groupby(level=0).mean(); sz = rate.index.values
near = sz[np.abs(np.log(D.n_human.values[:, None] / sz[None, :])).argmin(1)]
S = json.load(open(A / "hotspot_out" / "META_all_FSU_top1000_n3000_k30_t1000.json")); hs = S["sets"]["hs_top1000_min30"]
RISK = {g: set(pd.read_csv(B / "pair_q03389_q03123_20261010" / "report2_update_20261010" / "fig45_function" / "final_modules.csv").query("genotype == @g and final_module == 'risk module'").genes.iloc[0].split()) for g in ["Npp53", "Npp53Rb1"]}
D["shared_Npp53_risk"] = D.genes.apply(lambda s: len(set(s.split()) & RISK["Npp53"])); D["shared_Npp53Rb1_risk"] = D.genes.apply(lambda s: len(set(s.split()) & RISK["Npp53Rb1"]))
D = D.sort_values("worst_p"); D.drop(columns="genes").to_csv(O / "pooled_candidates.csv")
summ = dict(hotspot_modules=len(hs["order"]), module_sizes={str(m): len(hs["modules"][str(m)]) for m in hs["order"]}, n_autocorr_fdr05=S["n_autocorr_fdr05"],
            candidates=len(D), passers=int(D.meets_rule.sum()), expected_by_chance=round(float(rate.loc[near].sum()), 3), null_rate_by_size={int(k): round(float(v), 4) for k, v in rate.items()})
json.dump(summ, open(O / "pooled_summary.json", "w"), indent=1)
top = D.head(TOP).copy(); top["rank"] = range(1, len(top) + 1); top.drop(columns="genes").to_csv(O / "pooled_top.csv")
HJ = json.load(open(A / "candidates_human.json")); json.dump({f"all|{c}": sorted(set(HJ[c])) for c in top.index}, open(O / "pooled_top_human.json", "w"))
json.dump({c: top.loc[c, "genes"].split() for c in top.index}, open(O / "pooled_top_mouse.json", "w"))
pd.set_option("display.width", 220); print(json.dumps(summ)); print(D.drop(columns="genes").round(3).head(10).to_string())

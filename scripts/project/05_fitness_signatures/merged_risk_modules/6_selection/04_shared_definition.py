"""Shared definition across genotypes: bcr arm, list FSU_top1000, Hotspot k30 / top 1000. Per genotype: all module types, and Q03123's
grouping only (create_modules min 30: modules + adjacent pairs). Counts, passers (median, SU2C + WCDT-JCO), nearest-size chance. Env: moscot."""
import json, numpy as np, pandas as pd
from pathlib import Path
O = Path(__file__).resolve().parent; SU, JC, WM = "SU2C-mCRPC-union", "WCDT-JCO", "WCDT-merged"; out, rows = {}, []
ok = lambda X: (X[f"HR|{SU}"] > 1) & (X[f"p|{SU}"] < .05) & (X[f"HR|{JC}"] > 1) & (X[f"p|{JC}"] < .05)
wide = lambda K: (lambda W: W.set_axis([f"{a}|{b}" for a, b in W.columns], axis=1))(K[K.split == "median"].pivot_table(index="cand", columns="cohort", values=["HR", "p"]))
for g in ["Npp53", "Npp53Rb1"]:
    A = O.parent / f"metaNorm_20261010_{g}_bcr"; M = pd.read_csv(A / "candidates_meta.csv"); M = M[M.direction == "FSU"].set_index("cand"); W = wide(pd.read_csv(A / "km_scan.csv"))
    Rw = wide(pd.read_csv(A / "km_null.csv")); rate = pd.Series(ok(Rw).values, index=Rw.index.str.extract(r"rand_s(\d+)_")[0].astype(int).values).groupby(level=0).mean(); sz = rate.index.values
    for lab, pat in [("same_run_any_type", r"FSU_top1000_n3000\|k30t1000\|"), ("hotspot_min30_modules_and_pairs", r"FSU_top1000_n3000\|k30t1000\|hs_top1000_min30\|")]:
        c = M[M.configs.str.contains(pat)]; w = W.loc[c.index]; p = ok(w); near = sz[np.abs(np.log(c.n_human.values[:, None] / sz[None, :])).argmin(1)]
        out[f"{g}|{lab}"] = dict(fsu_candidates=len(c), passers=int(p.sum()), expected=round(float(rate.loc[near].sum()), 3))
        for cand in w.index[p]:
            rows.append(dict(genotype=g, set=lab, cand=cand, n_mouse=int(c.loc[cand, "n_mouse"]), module=";".join(x.split("|", 2)[2] for x in c.loc[cand, "configs"].split(";") if "k30t1000" in x),
                             SU2C=f"{w.loc[cand, f'HR|{SU}']:.2f} ({w.loc[cand, f'p|{SU}']:.3f})", JCO=f"{w.loc[cand, f'HR|{JC}']:.2f} ({w.loc[cand, f'p|{JC}']:.3f})",
                             merged=f"{w.loc[cand, f'HR|{WM}']:.2f} ({w.loc[cand, f'p|{WM}']:.3f})", worst_p=round(float(max(w.loc[cand, f'p|{SU}'], w.loc[cand, f'p|{JC}'])), 4)))
    out[f"{g}|autocorr_fdr05"] = json.load(open(A / "hotspot_out" / f"META_{g}_FSU_top1000_n3000_k30_t1000.json"))["n_autocorr_fdr05"]
P = pd.DataFrame(rows).drop_duplicates(["genotype", "set", "cand"]); P.to_csv(O / "shared_definition_passers.csv", index=False); json.dump(out, open(O / "shared_definition.json", "w"), indent=1)
print(json.dumps(out, indent=1)); print(P.to_string(index=False))

"""opt5 step 08 — direction-aware selection (usage: 08_select.py <arm>) on SU2C + WCDT-merged (primary; WCDT-JCO reported).
Risk reading: FSU / BI high score = adverse (HR_adv = HR); FSD low score = adverse (HR_adv = 1 / HR).
Rules, same split, HR_adv > 1 in both: strict = SU2C p < 0.05 & WCDT-merged p < 0.05; relaxed = SU2C p < 0.05 & WCDT-merged p < 0.1;
strict_jco / relaxed_jco = the same with WCDT-JCO. Splits: median, mean, tertile, q1q4, cont, optimal.
Output <arm>/selection.csv (rows passing any rule), <arm>/param_diagnostics.csv (median-split risk-z by parameter). Conda env: moscot."""
import sys, numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import norm
O = Path(__file__).resolve().parent; ARM = sys.argv[1]; A = O / ARM
K = pd.read_csv(A / "km_scan.csv"); M = pd.read_csv(A / "candidates_meta.csv")
K = K.merge(M[["cand", "direction"]], on="cand"); K["HR_adv"] = np.where(K.direction == "FSD", 1 / K.HR, K.HR)
W = K.pivot_table(index=["cand", "split"], columns="cohort", values=["HR_adv", "p"]); W.columns = [f"{a}|{b}" for a, b in W.columns]; W = W.reset_index()
su, wm, wj = "SU2C-mCRPC-union", "WCDT-merged", "WCDT-JCO"
adv = lambda c, a: (W[f"HR_adv|{c}"] > 1) & (W[f"p|{c}"] < a)
W["strict"] = adv(su, .05) & adv(wm, .05); W["relaxed"] = adv(su, .05) & adv(wm, .1); W["strict_jco"] = adv(su, .05) & adv(wj, .05); W["relaxed_jco"] = adv(su, .05) & adv(wj, .1)
S = W[W.relaxed | W.relaxed_jco].merge(M[["cand", "direction", "genotype", "n_mouse", "n_human", "n_configs", "n_lists", "first_config", "genes"]], on="cand")
S["split"] = pd.Categorical(S.split, ["median", "mean", "tertile", "q1q4", "cont", "optimal"], ordered=True); S = S.sort_values(["direction", "split", f"p|{su}"])
S.to_csv(A / "selection.csv", index=False)
print(f"{ARM}: candidates {M.cand.nunique()} {M.groupby(['genotype', 'direction']).size().to_dict()}")
for r in ["strict", "relaxed", "strict_jco", "relaxed_jco"]:
    X = W[W[r]].merge(M[["cand", "genotype", "direction"]], on="cand"); mm = X[X.split.isin(["median", "mean"])]
    print(f"  {r:12s} any split {X.cand.nunique():3d} {X.groupby(['genotype', 'direction']).cand.nunique().to_dict()} | median/mean {mm.cand.nunique()} {mm.groupby(['genotype', 'direction']).cand.nunique().to_dict()}")
K["z"] = np.sign(np.log(K.HR_adv)) * norm.isf(K.p.clip(1e-300) / 2); Z = K[K.split == "median"].pivot_table(index="cand", columns="cohort", values="z")
X = M[["cand", "n_mouse", "configs", "direction"]].assign(cfg=lambda d: d.configs.str.split(";")).explode("cfg"); p = X.cfg.str.split("|", expand=True)
X["list"], X["source"] = p[0], p[1]; q = X.list.str.split("_", expand=True); X["method"], X["genotype"], X["rule"] = q[0], q[1], q[3]
X["size"] = pd.cut(X.n_mouse, [0, 15, 30, 60, 120, 100000], labels=["<=15", "16-30", "31-60", "61-120", ">120"]); X = X.merge(Z, left_on="cand", right_index=True); out = []
for par in ["genotype", "direction", "method", "rule", "source", "size"]:
    for lev, g in X.groupby(par, observed=True):
        out.append(dict(param=par, level=lev, n=len(g), mean_z_su2c=g[su].mean(), mean_z_wcdt_merged=g[wm].mean(), mean_z_wcdt_jco=g[wj].mean(), r_su2c_wcdtm=g[su].corr(g[wm]),
                        frac_both_adverse=((g[su] > 0) & (g[wm] > 0)).mean()))
pd.DataFrame(out).to_csv(A / "param_diagnostics.csv", index=False)

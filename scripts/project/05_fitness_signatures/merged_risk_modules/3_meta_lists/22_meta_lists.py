"""opt5 step 22 — meta-analysis workflow, step 3 of 3: random-effects meta-analysis across tumours -> FSU / FSD lists, replacing the vote
(usage: 22_meta_lists.py <arm> [se: cr | ols]). Per gene and tumour set (Npp53 25 tumours; Npp53Rb1 22 incl. JZ135; all 47):
  input        beta and its subclone-robust SE per tumour (step 21); gene must be tested in >= max(5, half) of the set's tumours
  pooling      DerSimonian-Laird tau^2; random-effects mean mu with weights 1 / (SE^2 + tau^2)
  test         modified Hartung-Knapp-Sidik-Jonkman: Var(mu) = max(1, q) / sum(w*), q = sum w*(beta - mu)^2 / (k - 1); t with k - 1 df
               (recommended small-sample random-effects test); BH FDR over genes within the set
  also         I^2 heterogeneity; 95% prediction interval mu +- t_{k-2} sqrt(tau^2 + Var(mu)) (effect expected in a new tumour)
  stability    leave-one-tumour-out: refit without each tumour; stable = same sign and p < 0.05 in every refit
Noise genes removed before listing (Zejian): §9 non-coding (Rik, Gm#, mt-, accession-like), tracing transgenes, cytosolic / mitochondrial
ribosomal proteins, haemoglobins, Malat1 / Neat1. Lists (FSU = mu > 0, FSD = mu < 0), unsupervised fixed rules:
  q05loo   FDR < 0.05 and leave-one-out stable          (round 1: empty in every set — no gene passes FDR < 0.1)
  q10loo   FDR < 0.10 and leave-one-out stable
  q05loopi q05loo and the 95% prediction interval excludes 0
  p01loo / p05loo   nominal p < 0.01 / 0.05 and the sign stable in every leave-one-out refit
  top150 / top300 / top500 / top1000   genes with the largest |t| whose sign is stable in every leave-one-out refit (ranked, unsupervised)
Lists of 30-1500 genes go to Hotspot (same list format as opt5 step 03; method "META"). Output <arm>/meta_<set>.csv.gz, fsu_lists.json,
bidir_lists.json, fsu_grid.csv, unit_z_META.csv.gz (beta / SE per tumour, for the fitness-profile modules). Conda env: moscot."""
import re, sys, json, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
O = Path(__file__).resolve().parent; ARM = sys.argv[1]; A = O / ARM; SEC = "se_" + (sys.argv[2] if len(sys.argv) > 2 else "cr")
EXCL = [re.compile(p) for p in [r"Rik$", r"^Gm\d+", r"^mt-", r"^[A-Z]{2}\d{6}$", r"^Rp[ls]\d", r"^Rplp\d", r"^Rpsa$", r"^Mrp[ls]\d", r"^Hb[ab]-", r"^(Malat1|Neat1)$"]]
TG = {"mCherry", "Cas9", "GFP", "Cre"}; excl = lambda g: g in TG or any(p.search(g) for p in EXCL)
E = pd.read_csv(A / "effects.csv.gz"); E = E[~E.gene.map(excl)]
B = E.pivot_table(index="unit", columns="gene", values="beta"); S = E.pivot_table(index="unit", columns="gene", values=SEC).reindex_like(B)
(B / S).to_csv(A / "unit_z_META.csv.gz", compression="gzip")
T = pd.read_csv(A / "cells_n3000_PRIVATE.csv", index_col=0).groupby("sample_id").model.first()
SETS = {"Npp53": list(T.index[T == "Npp53"]), "Npp53Rb1": list(T.index[T == "Npp53Rb1"]), "all": list(T.index)}
import os
if os.environ.get("META_SETS"): SETS = {k: v for k, v in SETS.items() if k in os.environ["META_SETS"].split(",")}   # round 2+: restrict sets
def remeta(b, v):
    """b, v: (k, G) with NaN for untested; returns mu, se_hk, t, p, tau2, I2, k, pi_lo, pi_hi."""
    m = ~np.isnan(b); k = m.sum(0); w = np.where(m, 1 / v, 0); bz = np.where(m, b, 0)
    sw = w.sum(0); mu_fe = (w * bz).sum(0) / sw; Q = (w * (bz - mu_fe) ** 2 * m).sum(0); Cc = sw - (w ** 2).sum(0) / sw
    with np.errstate(divide="ignore", invalid="ignore"):
        tau2 = np.maximum(0, (Q - (k - 1)) / Cc); ws = np.where(m, 1 / (v + tau2), 0); sws = ws.sum(0); mu = (ws * bz).sum(0) / sws
        q = (ws * (bz - mu) ** 2 * m).sum(0) / (k - 1); se = np.sqrt(np.maximum(1, q) / sws); t = mu / se; p = 2 * stats.t.sf(np.abs(t), k - 1)
        I2 = np.where(Q > 0, np.maximum(0, (Q - (k - 1)) / Q), 0); half = stats.t.ppf(.975, np.maximum(k - 2, 1)) * np.sqrt(tau2 + 1 / sws)
    return mu, se, t, p, tau2, I2, k, mu - half, mu + half
rows, lists, kept, bidir = [], {}, {}, {}
for sname, units in SETS.items():
    b = B.loc[units].values; v = (S.loc[units].values) ** 2; kmin = max(5, int(np.ceil(len(units) / 2)))
    mu, se, t, p, tau2, I2, k, lo, hi = remeta(b, v); ok = k >= kmin
    loo_sign, loo_p = np.ones(b.shape[1], bool), np.zeros(b.shape[1])
    for u in range(len(units)):
        bb, vv = b.copy(), v.copy(); bb[u] = np.nan; vv[u] = np.nan; m2, _, _, p2, *_ = remeta(bb, vv)
        loo_sign &= np.sign(m2) == np.sign(mu); loo_p = np.maximum(loo_p, np.nan_to_num(p2, nan=1.0))
    M = pd.DataFrame(dict(gene=B.columns, k=k, mu=mu, se_hk=se, t=t, p=p, tau2=tau2, I2=I2, pi_lo=lo, pi_hi=hi, loo_sign_stable=loo_sign, loo_max_p=loo_p))[ok].copy()
    M["q"] = stats.false_discovery_control(M.p.fillna(1).values); M["loo_stable"] = M.loo_sign_stable & (M.loo_max_p < .05)
    M.sort_values("t", ascending=False).to_csv(A / f"meta_{sname}.csv.gz", index=False, compression="gzip")
    got = {}
    for dr, sg in [("FSU", 1), ("FSD", -1)]:
        d = M[np.sign(M.mu) == sg]
        got[(dr, "q05loo")] = d[(d.q < .05) & d.loo_stable]; got[(dr, "q10loo")] = d[(d.q < .10) & d.loo_stable]
        pi_ok = (d.pi_lo > 0) if sg == 1 else (d.pi_hi < 0); got[(dr, "q05loopi")] = d[(d.q < .05) & d.loo_stable & pi_ok]
        st = d[d.loo_sign_stable]; st = st.reindex(st.t.abs().sort_values(ascending=False).index)
        got[(dr, "p01loo")] = st[st.p < .01]; got[(dr, "p05loo")] = st[st.p < .05]
        for N in (150, 300, 500, 1000, 2000): got[(dr, f"top{N}")] = st.head(N)
    print(f"{sname}: genes tested {len(M)} | " + " | ".join(f"{dr} {r}: {len(x)}" for (dr, r), x in got.items()))
    for (dr, rule), d in got.items():
        lid = f"META_{sname}_{dr}_{rule}_n3000"; genes = sorted(d.gene); key = (sname, dr); kept.setdefault(key, [])
        r = dict(list_id=lid, method="META", genotype=sname, direction=dr, jz135="n3000", rule=rule, n_genes=len(genes), status="")
        if not 30 <= len(genes) <= 2000: r["status"] = "size out of range"
        else:
            s = set(genes); dup = next((kk for kk in kept[key] if len(s & lists[kk]["set"]) / len(s | lists[kk]["set"]) >= .9), None)
            if dup: r["status"] = f"alias of {dup}"; lists[dup]["aliases"].append(lid)
            else:
                r["status"] = "kept"; kept[key].append(lid)
                lists[lid] = dict(method="META", genotype=sname, direction=dr, jz135="n3000", cells=sname, rule=rule, genes=genes, set=s, aliases=[],
                                  meta_z={g: round(float(x), 3) for g, x in zip(d.gene, d.t)}, units=units)
        rows.append(r)
    for rule in ["q05loo", "q10loo", "q05loopi", "p01loo", "p05loo", "top150", "top300", "top500", "top1000", "top2000"]:
        up, dn = sorted(got[("FSU", rule)].gene), sorted(got[("FSD", rule)].gene)
        if 30 <= len(up) <= 1500 and 30 <= len(dn) <= 1500: bidir[f"META_{sname}_BI_{rule}_n3000"] = dict(up=up, dn=dn)
G = pd.DataFrame(rows); G.to_csv(A / "fsu_grid.csv", index=False)
json.dump({k: {kk: vv for kk, vv in d.items() if kk != "set"} for k, d in lists.items()}, open(A / "fsu_lists.json", "w"), allow_nan=False); json.dump(bidir, open(A / "bidir_lists.json", "w"))
print(G.to_string(index=False)); print(f"{ARM}: kept lists for Hotspot {(G.status == 'kept').sum()}")

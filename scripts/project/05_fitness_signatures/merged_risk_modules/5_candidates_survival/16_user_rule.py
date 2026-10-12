"""opt5 step 16 — Zejian's own selection rule only (usage: 16_user_rule.py <arm>). Zejian 2026-10-09: the family / JZ135-tree /
external-validation filters "are not my rules". Rule:
  survival : FSU module; SU2C HR > 1 with p < 0.05 AND a WCDT version (WCDT-merged or WCDT-JCO; "any WCDT version is fine") HR > 1 with
             p < 0.1, under the same split (median / mean preferred; other splits reported);
  mouse    : involved in Castration and/or Castration+ENZA (not necessarily specific): per-cell score (scanpy score_genes, log1p CP10k,
             random_state 42) on the arm's cells (cycOUT: non-cycling), per-tumour mean, z-scored; OLS z ~ treatment (ref Intact) + model
             over the 47 tumours, and z ~ treatment within the module's own genotype. Involved = Castration or Castration+ENZA higher
             than Intact with p < 0.05 in either model.
Output <arm>/user_rule_passers.csv (one row per passing module, best split = median/mean first, then smallest SU2C p). Env: moscot."""
import sys, numpy as np, pandas as pd, scanpy as sc, anndata as ad, statsmodels.formula.api as smf
from pathlib import Path
O = Path(__file__).resolve().parent; ARM = sys.argv[1]; A = O / ARM; MERGED = O.parents[5]
su, wm, wj = "SU2C-mCRPC-union", "WCDT-merged", "WCDT-JCO"
M = pd.read_csv(A / "candidates_meta.csv"); K = pd.read_csv(A / "km_scan.csv"); K = K[K.cand.isin(M[M.direction == "FSU"].cand)]
W = K.pivot_table(index=["cand", "split"], columns="cohort", values=["HR", "p"]); W.columns = [f"{a}|{b}" for a, b in W.columns]; W = W.reset_index()
ok = lambda c, a: (W[f"HR|{c}"] > 1) & (W[f"p|{c}"] < a)
W["pass_merged"] = ok(su, .05) & ok(wm, .1); W["pass_jco"] = ok(su, .05) & ok(wj, .1); P = W[W.pass_merged | W.pass_jco].copy()
P["pref"] = P.split.map({"median": 0, "mean": 0}).fillna(1); P = P.sort_values(["cand", "pref", f"p|{su}"])
best = P.groupby("cand").head(1).set_index("cand"); allsp = P.groupby("cand").split.apply(lambda s: ",".join(s))
print(f"{ARM}: FSU modules passing the survival rule: {len(best)} (median/mean: {(best.pref == 0).sum()})")
cells = pd.read_csv(A / "cells_n3000_PRIVATE.csv", index_col=0)
sh = ad.read_h5ad(MERGED / "shared_counts.h5ad"); sh = sh[cells.index.values].copy(); sh.var_names = sh.var.gene_symbol.astype(str).values; sh.var_names_make_unique()
sc.pp.normalize_total(sh, target_sum=1e4); sc.pp.log1p(sh); rows = []
meta = M.set_index("cand")
for c, b in best.iterrows():
    genes = [g for g in meta.loc[c, "genes"].split() if g in sh.var_names]
    sc.tl.score_genes(sh, genes, score_name="s", random_state=42)
    d = cells.assign(s=sh.obs["s"].values).groupby("sample_id").agg(s=("s", "mean"), treatment=("treatment", "first"), model=("model", "first")).reset_index()
    d["trt"] = pd.Categorical(d.treatment, categories=["Intact", "Castration", "Castration+ENZA"]); d["z"] = (d.s - d.s.mean()) / d.s.std()
    f = smf.ols("z ~ trt + model", data=d).fit(); g = meta.loc[c, "genotype"]; dg = d[d.model == g]
    fg = smf.ols("z ~ trt", data=dg).fit() if len(dg) > 5 else f                  # pooled ("all") modules: own-genotype model = all-tumour model
    co = lambda fit, t: (fit.params.get(f"trt[T.{t}]", np.nan), fit.pvalues.get(f"trt[T.{t}]", np.nan))
    (bc, pc), (be, pe), (bcg, pcg), (beg, peg) = co(f, "Castration"), co(f, "Castration+ENZA"), co(fg, "Castration"), co(fg, "Castration+ENZA")
    inv = [lab for lab, bb, pp in [("Cas (all)", bc, pc), ("Cas+ENZA (all)", be, pe), (f"Cas ({g})", bcg, pcg), (f"Cas+ENZA ({g})", beg, peg)] if bb > 0 and pp < .05]
    rows.append(dict(cand=c, genotype=g, n_genes=meta.loc[c, "n_mouse"], best_split=b.split, splits_passing=allsp[c],
                     SU2C=f"{b[f'HR|{su}']:.2f} ({b[f'p|{su}']:.3f})", WCDT_merged=f"{b[f'HR|{wm}']:.2f} ({b[f'p|{wm}']:.3f})", WCDT_JCO=f"{b[f'HR|{wj}']:.2f} ({b[f'p|{wj}']:.3f})",
                     Cas_all=f"{bc:+.2f} (p {pc:.3f})", CasENZA_all=f"{be:+.2f} (p {pe:.3f})", Cas_own=f"{bcg:+.2f} (p {pcg:.3f})", CasENZA_own=f"{beg:+.2f} (p {peg:.3f})",
                     involved_Cas_or_CasENZA="; ".join(inv) if inv else "no", genes=" ".join(sorted(meta.loc[c, "genes"].split()))))
R = pd.DataFrame(rows); R["rule_met"] = R.involved_Cas_or_CasENZA != "no"; R.to_csv(A / "user_rule_passers.csv", index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_colwidth", 60)
print(R.drop(columns=["genes"]).to_string(index=False))
print(f"{ARM}: FSU modules meeting the full rule (survival + Cas/Cas+ENZA involvement): {int(R.rule_met.sum())}")

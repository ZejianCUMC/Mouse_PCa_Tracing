"""opt5 step 25 — one fixed survival rule for both genotypes (Zejian 2026-10-10: "find risk modules without strict castration enriched;
method to get modules in Npp53 and Npp53Rb1 should be very similar; survival splitting method and final picked cohorts must be the
same, like split patients by median and pick SU2C and WCDT-JCO").
Usage: 25_fixed_rule.py <out_dir> <arm> [<arm> ...]     env SPLIT (default median), COHORTS (default "SU2C-mCRPC-union,WCDT-JCO"), ALPHA (0.05)
Rule : FSU module; in EVERY listed cohort HR > 1 and log-rank p < ALPHA under SPLIT (high vs low score). No castration requirement;
       castration effects (16_user_rule.py; computed for modules passing the earlier p < 0.1 rule, a superset) are information only.
Chance: the arm's 4,000 random gene sets (km_null.csv; 1,000 each of 15/30/60/120 human genes from the arm's list-gene ortholog pool)
       read with the same rule; expected = sum over the arm's FSU candidates (per genotype) of the pass rate at the nearest size.
Modules: identical gene sets across arms merged; families = connected components at Jaccard >= 0.5; ranking = larger of the cohort p.
Output <out_dir>/funnel.csv, passers.csv (one row per distinct gene set), families.csv, cumulative.csv. Env: moscot."""
import os, sys, numpy as np, pandas as pd
from pathlib import Path
O = Path(__file__).resolve().parent; OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True); ARMS = sys.argv[2:]
SPLIT = os.environ.get("SPLIT", "median"); CO = os.environ.get("COHORTS", "SU2C-mCRPC-union,WCDT-JCO").split(","); AL = float(os.environ.get("ALPHA", .05))
SE = {"metaR1": "cr", "metaR2b": "cr", "metaR2z": "cr", "metaR3b": "ols", "metaR3z": "ols"}; FIT = {"metaR1": "rank", "metaR2b": "rank", "metaR2z": "z", "metaR3b": "rank", "metaR3z": "z"}
GRID = {"metaR1": "k30 x top500"}
for a in ARMS:
    if a.startswith(("metaRb1_", "metaNorm_")): s = a.rsplit("_", 1)[1]; FIT[a] = "rank" if s[0] == "b" else "z"; SE[a] = s[1:]
def passes(K):
    K = K[K.split == SPLIT]; W = K.pivot_table(index="cand", columns="cohort", values=["HR", "p"]); W.columns = [f"{x}|{y}" for x, y in W.columns]
    ok = np.logical_and.reduce([(W[f"HR|{c}"] > 1) & (W[f"p|{c}"] < AL) for c in CO]); W["worst_p"] = W[[f"p|{c}" for c in CO]].max(1); return W, ok
fun, rows = [], []
for arm in ARMS:
    A = O / arm; M = pd.read_csv(A / "candidates_meta.csv").set_index("cand"); M = M[M.direction == "FSU"]
    K = pd.read_csv(A / "km_scan.csv"); W, ok = passes(K[K.cand.isin(M.index)]); Wp = W[ok]
    R = pd.read_csv(A / "km_null.csv"); Rw, rok = passes(R); size = Rw.index.str.extract(r"rand_s(\d+)_")[0].astype(int).values
    rate = pd.Series(rok, index=size).groupby(level=0).mean(); sizes = rate.index.values
    U = pd.read_csv(A / "user_rule_passers.csv").set_index("cand"); P = pd.read_csv(A / "prolif_corr.csv").set_index("cand")
    K2 = K.set_index(["cand", "split", "cohort"])
    for geno, g in M.groupby("genotype"):
        near = sizes[np.abs(np.log(g.n_human.values[:, None] / sizes[None, :])).argmin(1)]
        fun.append(dict(arm=arm, genotype=geno, fitness=FIT[arm], se=SE[arm], fsu_candidates=len(g), passers=int(Wp.index.isin(g.index).sum()),
                        expected_by_chance=round(float(rate.loc[near].sum()), 2), null_rate_by_size=";".join(f"{s}:{r:.4f}" for s, r in rate.items())))
    for c in Wp.index:
        hr = lambda co, sp="median": f"{K2.loc[(c, sp, co), 'HR']:.2f} ({K2.loc[(c, sp, co), 'p']:.3f})" if (c, sp, co) in K2.index else ""
        rows.append(dict(arm=arm, fitness=FIT[arm], se=SE[arm], cand=c, genotype=M.loc[c, "genotype"], n_genes=int(M.loc[c, "n_mouse"]), n_configs=int(M.loc[c, "n_configs"]),
                         worst_p=Wp.loc[c, "worst_p"], SU2C=hr("SU2C-mCRPC-union", SPLIT), WCDT_JCO=hr("WCDT-JCO", SPLIT), WCDT_merged=hr("WCDT-merged", SPLIT),
                         SU2C_mean=hr("SU2C-mCRPC-union", "mean"), JCO_mean=hr("WCDT-JCO", "mean"), SU2C_cont=hr("SU2C-mCRPC-union", "cont"), JCO_cont=hr("WCDT-JCO", "cont"),
                         castration_info=U.loc[c, "involved_Cas_or_CasENZA"] if c in U.index else "not tested",
                         Cas_all=U.loc[c, "Cas_all"] if c in U.index else "", CasENZA_all=U.loc[c, "CasENZA_all"] if c in U.index else "",
                         Cas_own=U.loc[c, "Cas_own"] if c in U.index else "", CasENZA_own=U.loc[c, "CasENZA_own"] if c in U.index else "",
                         rho_E2F_SU2C=round(P.loc[c, "rho_SU2C"], 2), padj_E2F_SU2C=round(P.loc[c, "padj_SU2C"], 3), padj_E2F_WCDTm=round(P.loc[c, "padj_WCDTm"], 3),
                         genes=M.loc[c, "genes"]))
F = pd.DataFrame(fun); F.to_csv(OUT / "funnel.csv", index=False)
C = F.groupby("genotype").agg(arms=("arm", "nunique"), fsu_candidates=("fsu_candidates", "sum"), passers=("passers", "sum"), expected_by_chance=("expected_by_chance", "sum")).reset_index()
C.to_csv(OUT / "cumulative.csv", index=False)
R = pd.DataFrame(rows); R["gset"] = R.genes.map(lambda x: frozenset(x.split())); out, fam = [], []
for geno, g in R.groupby("genotype"):
    agg = []
    for s, d in g.groupby("gset"):
        d = d.sort_values("worst_p"); b = d.iloc[0].to_dict()
        b.update(ids=";".join(f"{a}:{c}" for a, c in zip(d.arm, d.cand)), n_arms=d.arm.nunique(), configs_total=int(d.n_configs.sum()), se_types=",".join(sorted(set(d.se)))); agg.append(b)
    G = pd.DataFrame(agg).sort_values("worst_p").reset_index(drop=True); S = list(G.gset); par = list(range(len(G)))
    def find(i):
        while par[i] != i: i = par[i]
        return i
    for i in range(len(G)):
        for j in range(i + 1, len(G)):
            if len(S[i] & S[j]) / len(S[i] | S[j]) >= .5: par[find(j)] = find(i)
    G["family"] = [find(i) for i in range(len(G))]; fm = {f: f"{geno}-F{k+1}" for k, f in enumerate(pd.unique(G.family))}; G["family"] = G.family.map(fm)
    for f, d in G.groupby("family", sort=False):
        arms = sorted({x.split(":")[0] for i in d.ids for x in i.split(";")}); core = frozenset.intersection(*d.gset)
        fam.append(dict(genotype=geno, family=f, members=len(d), lead=d.ids.iloc[0], lead_genes=int(d.n_genes.iloc[0]), lead_worst_p=round(d.worst_p.iloc[0], 4),
                        arms=";".join(arms), n_arms=len(arms), se_types=",".join(sorted({SE[a] for a in arms})), configs_total=int(d.configs_total.sum()),
                        core_n=len(core), core_genes=" ".join(sorted(core)) if len(core) <= 40 else ""))
    G["family_members"] = G.family.map(G.family.value_counts()); out.append(G)
T = pd.concat(out).drop(columns=["gset"]) if out else pd.DataFrame(); T.to_csv(OUT / "passers.csv", index=False); Fa = pd.DataFrame(fam); Fa.to_csv(OUT / "families.csv", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 60)
print(f"rule: {SPLIT} split; HR > 1 and p < {AL} in {CO}"); print(F.drop(columns=["null_rate_by_size"]).to_string(index=False)); print(C.to_string(index=False))
if len(Fa): print(Fa.to_string(index=False))
for geno, d in (T.groupby("genotype") if len(T) else []):
    print(f"\n== {geno}: {len(d)} distinct passers")
    print(d.head(12)[["ids", "family", "n_genes", "SU2C", "WCDT_JCO", "WCDT_merged", "rho_E2F_SU2C", "padj_E2F_SU2C", "configs_total", "n_arms", "se_types", "castration_info"]].to_string(index=False))

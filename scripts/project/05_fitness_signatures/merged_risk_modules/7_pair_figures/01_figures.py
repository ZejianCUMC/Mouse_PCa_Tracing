"""Figures and tables for the matched pair Npp53 Q03389 / Npp53Rb1 Q03123 (identical setting: CP10k bcr arm, list FSU_top1000, Hotspot k30 /
top 1000, create_modules min 30, adjacent pair). Read-only on saved outputs. Env: moscot."""
import json, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
O = Path(__file__).resolve().parent; B = O.parent; sys.path.insert(0, str(B)); import orthologs as ORT
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
PAIR = [("Npp53", "Q03389", (7, 6)), ("Npp53Rb1", "Q03123", (6, 8))]; K = []; F = {}
for g, cand, mods in PAIR:
    A = B / f"metaNorm_20261010_{g}_bcr"; lid = f"META_{g}_FSU_top1000_n3000"; M = pd.read_csv(A / f"meta_{g}.csv.gz")
    cm = pd.read_csv(A / "candidates_meta.csv").set_index("cand"); genes = cm.loc[cand, "genes"].split()
    H = json.load(open(A / "hotspot_out" / f"{lid}_k30_t1000.json")); S = H["sets"]["hs_top1000_min30"]
    Z = np.load(A / "hotspot_out" / f"{lid}_k30_t1000_lcz.npy"); zg = (A / "hotspot_out" / f"{lid}_k30_t1000_lcz_genes.txt").read_text().split(); pos = {x: i for i, x in enumerate(zg)}
    ordg, bounds = [], []
    for m in S["order"]:
        mg = sorted(S["modules"][str(m)]); bounds.append((int(m), len(ordg), len(ordg) + len(mg))); ordg += mg
    fig, ax = plt.subplots(figsize=(7.0, 5.8)); im = ax.imshow(Z[np.ix_([pos[x] for x in ordg], [pos[x] for x in ordg])], cmap="RdBu_r", vmin=-10, vmax=10, interpolation="nearest")
    for m, a, b in bounds:
        hit = m in mods; ax.add_patch(plt.Rectangle((a - .5, a - .5), b - a, b - a, fill=False, ec="k" if hit else "#888", lw=1.6 if hit else .4))
        if hit: ax.text(a - 4, (a + b) / 2, f"M{m}", ha="right", fontsize=8, va="center")
    ax.set_xticks([]); ax.set_yticks([]); fig.colorbar(im, ax=ax, shrink=.7, label="local-correlation Z")
    ax.set_title(f"{g}: Hotspot k30, top {len(zg)} autocorrelated genes of the top-1000 list\n{len(S['order'])} modules (min 30) in dendrogram order; {cand} = M{mods[0]} + M{mods[1]} ({len(genes)} genes)", fontsize=8)
    fig.tight_layout(); fig.savefig(O / f"fig_heatmap_{cand}.png", dpi=150); plt.close(fig)
    lp = -np.log10(M.p.clip(lower=1e-300)); fig, ax = plt.subplots(figsize=(6, 4.2)); ax.scatter(M.mu, lp, s=2, c="#c8c8c8", rasterized=True, label=f"all tested genes ({len(M):,})")
    t1 = M.gene.isin(json.load(open(A / "fsu_lists.json"))[lid]["genes"]); ax.scatter(M.mu[t1], lp[t1], s=3, c="#9ecae1", rasterized=True, label="top-1000 fitness-up list")
    q = M.gene.isin(genes); ax.scatter(M.mu[q], lp[q], s=10, c="#c0392b", label=f"{cand} ({int(q.sum())})"); ax.axhline(2, color="grey", lw=.5, ls="--")
    ax.set_xlabel("pooled effect per 1 SD fitness (log1p CP10k)"); ax.set_ylabel("-log10 p (Hartung-Knapp)")
    ax.set_title(f"{g} meta-analysis across {M.k.max()} tumours (lowest q {M.q.min():.3f})", fontsize=9); ax.legend(fontsize=7, frameon=False, markerscale=2)
    fig.tight_layout(); fig.savefig(O / f"fig_volcano_{g}.png", dpi=150); plt.close(fig)
    Ms = M.set_index("gene"); sub = {x: f"M{m}" for m in mods for x in S["modules"][str(m)]}
    rows = [dict(gene=x, module=sub[x], human=";".join(ORT.mouse_to_human(x)[0]), pct_per_SD=round(100 * np.expm1(Ms.loc[x, "mu"]), 3), p=Ms.loc[x, "p"]) for x in genes]
    G = pd.DataFrame(rows).sort_values(["module", "gene"]); G.to_csv(O / f"{cand}_genes.csv", index=False)
    U = pd.read_csv(A / "user_rule_passers.csv").set_index("cand").loc[cand]; k = pd.read_csv(A / "km_scan.csv"); K.append(k[k.cand == cand])
    F[cand] = dict(genotype=g, n_mouse=len(genes), n_human=int(cm.loc[cand, "n_human"]), sizes={f"M{m}": len(S["modules"][str(m)]) for m in mods}, n_modules=len(S["order"]), autocorr=H["n_autocorr_fdr05"],
                   genes_p01=int((G.p < .01).sum()), genes_p05=int((G.p < .05).sum()), castration={x: U[x] for x in ["Cas_all", "CasENZA_all", "Cas_own", "CasENZA_own"]})
pd.concat(K).to_csv(O / "km_all_splits.csv", index=False); json.dump(F, open(O / "facts.json", "w"), indent=1); print(json.dumps(F, indent=1))

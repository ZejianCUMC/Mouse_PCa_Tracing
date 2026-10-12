"""Report-2 update (2026-10-10, thread 'Optimize result report'): redraw Fig 4.3 so per-tumor differences are visible (linear diverging
bars per tumor with the chance level, plus how many tumors share each call) and redraw Fig 4.2b with the genotype spelled Npp53Rb1.
Read-only on all inputs. Env: moscot."""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
O = Path(__file__).resolve().parent; B = O.parents[1]
plt.rcParams.update({"font.family": "DejaVu Sans"})
GENOS = ["Npp53", "Npp53Rb1"]; A = B / "metaNorm_20261010_Npp53_bcr"   # effects/cells tables cover both genotypes
UP, DN = "#cc1a1a", "#2b6cb0"; GCOL = {"Npp53": "#4a7fc1", "Npp53Rb1": "#c44e52"}; ALPHA = 0.01
cells = pd.read_csv(A / "cells_n3000_PRIVATE.csv", usecols=["cell", "model", "sample_id", "cycling"])
nc = cells.groupby(["model", "sample_id"]).agg(n=("cell", "size"), cyc=("cycling", "sum")).reset_index()
order = {g: nc[nc.model == g].sort_values("n", ascending=False).sample_id.tolist() for g in GENOS}   # same order as Fig 4.2b
E = pd.read_csv(A / "effects.csv.gz"); z = E.beta / E.se_cr; th = stats.norm.isf(ALPHA / 2)
E["up"] = z > th; E["dn"] = z < -th
S = E.groupby("unit").agg(tested=("gene", "size"), up=("up", "sum"), dn=("dn", "sum"))
S["chance"] = S.tested * ALPHA / 2
out = []
# ---- Fig 4.3 ----
fig = plt.figure(figsize=(15, 6.6)); gs = fig.add_gridspec(2, 2, width_ratios=[4.3, 1], hspace=.75, wspace=.12, top=.9)
ymax = max(S.up.max(), S.dn.max()) * 1.08
for r, g in enumerate(GENOS):
    ax = fig.add_subplot(gs[r, 0]); t = order[g]; d = S.loc[t]; x = np.arange(len(t))
    ax.bar(x, d.up, color=UP, width=.78); ax.bar(x, -d.dn, color=DN, width=.78)
    for i, ch in enumerate(d.chance): ax.hlines([ch, -ch], i - .39, i + .39, color="k", lw=1.1, ls=(0, (2, 1.2)))
    ax.axhline(0, color="k", lw=.8); ax.set_xlim(-.6, len(t) - .4); ax.set_ylim(-ymax, ymax)
    ax.set_xticks(x); ax.set_xticklabels(t, rotation=90, fontsize=7.5)
    ax.set_yticks(np.arange(-800, 801, 200)); ax.set_yticklabels([f"{abs(v)}" for v in np.arange(-800, 801, 200)], fontsize=8)
    ax.set_ylabel("genes at p < 0.01\n(down | up)", fontsize=9)
    ax.set_title(f"{'ab'[r]}  {g}: fitness genes per tumor ({len(t)} tumors, ordered by cells with LBI)", loc="left", fontsize=10.5)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    for _, row in d.assign(tumor=t).iterrows(): out.append(dict(genotype=g, tumor=row.tumor, genes_tested=int(row.tested), up=int(row.up), down=int(row.dn), chance_per_direction=round(row.chance, 1)))
fig.legend(handles=[Patch(color=UP, label="fitness-up"), Patch(color=DN, label="fitness-down"), plt.Line2D([], [], color="k", ls=(0, (2, 1.2)), label="expected by chance (0.5% of tested genes)")],
                   fontsize=8.5, frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(.5, .985))
share = []
for r, g in enumerate(GENOS):
    ax = fig.add_subplot(gs[r, 1]); e = E[E.unit.isin(order[g])]; k = np.arange(1, 7)
    for j, (col, c) in enumerate([("up", UP), ("dn", DN)]):
        n = e[e[col]].groupby("gene").size(); f = np.array([(n == i).mean() * 100 for i in k[:-1]] + [(n >= 6).mean() * 100])
        ax.bar(k + (j - .5) * .38, f, width=.38, color=c)
        share.append(dict(genotype=g, direction={"up": "fitness-up", "dn": "fitness-down"}[col], genes_called=len(n), **{f"pct_in_{i}": round(v, 1) for i, v in zip(["1", "2", "3", "4", "5", "6plus"], f)}))
    ax.set_xticks(k); ax.set_xticklabels(["1", "2", "3", "4", "5", "≥6"], fontsize=8); ax.set_ylim(0, 100); ax.tick_params(labelsize=8)
    ax.set_xlabel("tumors in which the gene is called", fontsize=8.5); ax.set_ylabel("% of called genes", fontsize=8.5)
    ax.set_title(f"{'cd'[r]}  {g}: calls shared across tumors", loc="left", fontsize=10.5)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
fig.savefig(O / "fig4_3_fitness_genes_per_tumor.png", dpi=160, bbox_inches="tight"); plt.close(fig); print("saved fig4_3")
pd.DataFrame(out).to_csv(O / "fig4_3_per_tumor_counts.csv", index=False); pd.DataFrame(share).to_csv(O / "fig4_3_sharing.csv", index=False)
# ---- Fig 4.2b ----
fig, ax = plt.subplots(figsize=(15, 4.2)); x0 = 0; ticks, labs = [], []
for g in GENOS:
    d = nc[nc.model == g].set_index("sample_id").loc[order[g]]; x = np.arange(x0, x0 + len(d))
    ax.bar(x, d.n, color=GCOL[g], width=.8, label=f"{g} ({len(d)} tumors)"); ax.bar(x, d.cyc, color="k", alpha=.45, width=.8)
    ticks += list(x); labs += list(d.index); x0 += len(d)
h, _ = ax.get_legend_handles_labels(); ax.legend(handles=h + [Patch(fc="k", alpha=.45, label="cycling cells")], fontsize=9, frameon=False)
ax.set_xticks(ticks); ax.set_xticklabels(labs, rotation=90, fontsize=7.5); ax.set_ylabel("cells with LBI")
ax.set_title(f"Cells with lineage fitness per tumor ({len(cells):,} cells; {int(cells.cycling.sum()):,} cycling, kept)", loc="left")
for s in ["top", "right"]: ax.spines[s].set_visible(False)
fig.savefig(O / "fig4_2b_cells_per_tumor.png", dpi=160, bbox_inches="tight"); plt.close(fig); print("saved fig4_2b")
print(pd.DataFrame(share).to_string())

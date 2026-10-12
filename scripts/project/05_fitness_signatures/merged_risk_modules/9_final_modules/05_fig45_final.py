"""Fig 4.5 with the final modules (03_merge_final.py): Hotspot local-correlation heatmap per genotype, genes ordered by final module
(Hotspot modules in dendrogram order inside each final module), thin bar = Hotspot modules, thick bar / side bar / boxes = final modules,
black box = risk module. Legend gives each final module's size and its top enriched term (04_go_final.R). Env: moscot."""
import json, numpy as np, pandas as pd
from pathlib import Path
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent)); from hallmark_names import hallmark   # local helper
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
O = Path(__file__).resolve().parent; B = O.parents[2]
plt.rcParams.update({"font.family": "DejaVu Sans"})
RCOL = "#cc1a1a"; FCOL = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2"]
F = pd.read_csv(O / "final_modules.csv"); GO = pd.read_csv(O / "go_bp_final_modules.csv"); HM = pd.read_csv(O / "hallmark_final_modules.csv")
def top_term(g, m):
    d = GO[(GO.genotype == g) & (GO.final_module == m)].sort_values("p.adjust")
    if len(d): return d.Description.iloc[0]
    h = HM[(HM.genotype == g) & (HM.final_module == m)].sort_values("p.adjust")
    if len(h): return hallmark(h.ID.iloc[0])
    return "no enriched term"
F["top_term"] = [top_term(r.genotype, r.final_module) for r in F.itertuples()]; F.drop(columns="genes").to_csv(O / "final_modules_annotated.csv", index=False)
fig = plt.figure(figsize=(20, 11.4)); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, .03], wspace=.3)
for col, (g, lab) in enumerate(zip(["Npp53", "Npp53Rb1"], "ab")):
    A = B / f"metaNorm_20261010_{g}_bcr" / "hotspot_out"; lid = f"META_{g}_FSU_top1000_n3000"
    S = json.load(open(A / f"{lid}_k30_t1000.json"))["sets"]["hs_top1000_min30"]; order = S["order"]
    Z = np.load(A / f"{lid}_k30_t1000_lcz.npy"); zg = (A / f"{lid}_k30_t1000_lcz_genes.txt").read_text().split(); pos = {x: i for i, x in enumerate(zg)}
    d = F[F.genotype == g].reset_index(drop=True); hmem = {x: m for m in order for x in S["modules"][str(m)]}
    ordg, mem = [], {}
    for r in d.itertuples():
        for m in [int(t[1:]) for t in r.hotspot_modules.split()]:
            for x in sorted(S["modules"][str(m)]): ordg.append(x); mem[x] = r.final_module
    ordg += sorted(set(zg) - set(mem)); n = len(ordg)
    fcol = {}; k = 0
    for r in d.itertuples():
        if r.final_module == "risk module": fcol[r.final_module] = RCOL
        else: fcol[r.final_module] = FCOL[k]; k += 1
    hcm = plt.get_cmap("tab20"); hcol = {m: hcm(i % 20) for i, m in enumerate(order)}
    ax = fig.add_subplot(gs[col]); ax.imshow(np.clip(Z[np.ix_([pos[x] for x in ordg], [pos[x] for x in ordg])], -8, 8), cmap="RdBu_r", vmin=-8, vmax=8, interpolation="nearest")
    for i, x in enumerate(ordg):
        ax.add_patch(Rectangle((i - .5, -60), 1, 18, color=hcol[hmem[x]] if x in hmem else (.88, .88, .88), clip_on=False))
        c = fcol[mem[x]] if x in mem else (.88, .88, .88)
        ax.add_patch(Rectangle((i - .5, -38), 1, 30, color=c, clip_on=False)); ax.add_patch(Rectangle((-32, i - .5), 24, 1, color=c, clip_on=False))
    ax.text(-40, -51, f"Hotspot ({len(order)})", ha="right", va="center", fontsize=8.5, color="#555"); ax.text(-40, -23, f"final ({len(d)})", ha="right", va="center", fontsize=10, weight="bold")
    for r in d.itertuples():
        ii = [i for i, x in enumerate(ordg) if mem.get(x) == r.final_module]; a0, a1 = min(ii), max(ii); risky = r.final_module == "risk module"
        ax.add_patch(Rectangle((a0 - .5, a0 - .5), a1 - a0 + 1, a1 - a0 + 1, fill=False, ec="black" if risky else fcol[r.final_module], lw=2.6 if risky else 2.2))
        if risky: ax.annotate(f"{g} risk module\n({r.n_genes} genes)", xy=(a0 - 4, (a0 + a1) / 2), xytext=(max(a0 - 420, 20), a1 - 40), fontsize=11, weight="bold", bbox=dict(fc="white", ec="black", lw=.8), arrowprops=dict(arrowstyle="->", lw=1.3))
    ax.legend(handles=[Patch(fc=fcol[r.final_module], label=f"{r.final_module} ({r.n_genes} genes): {r.top_term}") for r in d.itertuples()] + [Patch(fc=(.88, .88, .88), label=f"unassigned ({n - len(mem)} genes)")],
              loc="upper left", bbox_to_anchor=(0, -.01), ncol=1, frameon=False, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(-33, n - .5); ax.set_ylim(n - .5, -61)
    ax.set_title(f"{lab}  {g}: Hotspot local correlations of the top-1,000 fitness-up genes, ordered by final module\n    {len(order)} Hotspot modules merged into {len(d)} final modules; {n - len(mem)} genes unassigned (grey)", fontsize=10.5, loc="left", pad=8)
cax = fig.add_subplot(gs[2]); fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-8, 8), cmap="RdBu_r"), cax=cax).set_label("local-correlation Z (clipped ±8)")
fig.savefig(O / "fig4_5_final_modules.png", dpi=130, bbox_inches="tight"); plt.close(fig)
pd.set_option("display.width", 220); print(F.drop(columns="genes").to_string())

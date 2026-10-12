# =============================================================================
# Npp53 tree viewer: swap in the error-intBC-free JZ136 results
# =============================================================================
# Description : The Npp53 viewer has no builder (patched in place, ANALYSIS_RECORD). Replace ONLY JZ136:
#               trees (Greedy/NJ/Spectral/MaxCut/Percolation; Lineage Trajectory = Greedy), per-cell LBI 'f'
#               (new jungle), 'pt' (new v10-UMAP PhyloVelo, rank-normalised within sample = pt_rank rule of
#               combined_pseudotime_percell), 'plast' (Fitch per-cell plasticity, v10 function, min clade 20,
#               on the new Greedy tree). Clade tooltip data (plasticityFitness) for JZ136 is recomputed on the
#               new tree with viewer node ids (preorder, assignIds) so it matches the displayed tree.
#               All other samples untouched. Previous HTML backed up next to this run.
# Input       : ../{greedy,nj,maxcut,spectral,percolation,jungle_fitness,phylovelo_out}/, Npp53 viewer HTML
# Output      : reports/Npp53_tree_viewer_interactive.html (+ backup, PATCH_RECEIPT.json)
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python
# Validated   : 2026-10-04
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import json, re, shutil, hashlib, datetime, sys
from pathlib import Path
import numpy as np, pandas as pd
from ete3 import Tree

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
RUN = Path(__file__).resolve().parents[1]
N = RUN.parents[1]
DEFAULTS = {
    "sample": "JZ136",
    "viewer": N / "reports/Npp53_tree_viewer_interactive.html",
    "solvers": {"Greedy": "greedy", "NJ": "nj", "Spectral": "spectral", "MaxCut": "maxcut", "Percolation": "percolation"},
    "plast_min_clade": 20,      # v10 compute_fitch_plasticity_per_cell default
    "pt_round": 4,
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

def fitch(tree, state_map):
    """Fitch node states (v10 rebuild_sc_fitness_plasticity_v10.py logic)."""
    node_set = {}
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            s = state_map.get(n.name); node_set[n] = {s} if s else set()
        else:
            cs = [node_set[c] for c in n.children if node_set[c]]
            if not cs: node_set[n] = set(); continue
            inter = set.intersection(*cs); node_set[n] = inter if inter else set.union(*cs)
    st = {}
    for n in tree.traverse("preorder"):
        if n.is_root(): st[n] = next(iter(node_set[n])) if node_set[n] else None
        else:
            ps = st.get(n.up)
            st[n] = ps if ps and ps in node_set[n] else (next(iter(node_set[n])) if node_set[n] else ps)
    return st

def clade_plasticity(node, st):
    trans = edges = 0
    for sub in node.traverse():
        if sub is node or sub.up is None: continue
        if st[sub] is not None and st[sub.up] is not None:
            edges += 1; trans += st[sub] != st[sub.up]
    return trans / edges if edges else None

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}; s = p["sample"]
    html = p["viewer"].read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", html, re.S); D = json.loads(m.group(1))
    old_cells = set(re.findall(r"[(,]([ACGT]{16}):", D["trees"][s]["Greedy"]))
    trees = {}
    for k, d in p["solvers"].items():
        t = Tree(open(RUN / d / f"{s}_tree_{d}.nwk").read(), format=1)
        assert set(t.get_leaf_names()) == old_cells, f"{k}: leaf set differs from current viewer"
        trees[k] = t
    newick = {k: t.write(format=1) for k, t in trees.items()}
    lbi = pd.read_csv(RUN / "jungle_fitness" / f"{s}_lbi_fitness.csv")
    f = lbi[lbi.is_leaf].set_index("node").fitness
    pt_raw = pd.read_csv(RUN / "phylovelo_out" / f"{s}_phylovelo_pseudotime.csv", index_col=0).iloc[:, 0]
    pt = pt_raw.rank(pct=True).round(p["pt_round"])          # pt_rank rule
    assert set(f.index) == old_cells and np.isfinite(f).all()
    states = {c: D["meta"][c]["s"] for c in old_cells}
    g = trees["Greedy"]; st = fitch(g, {c: re.sub(r"^c\d+_", "", v) for c, v in states.items()})
    psum, pn = {}, {}
    for node in g.traverse("postorder"):
        if node.is_leaf() or len(node) < p["plast_min_clade"]: continue
        v = clade_plasticity(node, st)
        if v is None: continue
        for l in node.get_leaves(): psum[l.name] = psum.get(l.name, 0) + v; pn[l.name] = pn.get(l.name, 0) + 1
    plast = {c: psum[c] / pn[c] for c in psum}
    # clade tooltip data keyed by the viewer's preorder ids (assignIds) on the displayed Greedy tree
    ids = {n: i for i, n in enumerate(Tree(newick["Greedy"], format=1).traverse("preorder"))}
    gg = Tree(newick["Greedy"], format=1); st2 = fitch(gg, {c: re.sub(r"^c\d+_", "", v) for c, v in states.items()})
    cp, cs = {}, {}
    for i, n in enumerate(gg.traverse("preorder")):
        if n.is_leaf(): continue
        cs[str(i)] = len(n)
        v = clade_plasticity(n, st2)
        if v is not None: cp[str(i)] = round(v, 4)
    # apply
    for k in p["solvers"]: D["trees"][s][k] = newick[k]
    D["trees"][s]["Lineage Trajectory"] = newick["Greedy"]
    for c in old_cells:
        D["meta"][c]["f"] = float(f[c])
        if c in pt.index: D["meta"][c]["pt"] = float(pt[c])
        else: D["meta"][c].pop("pt", None)
        if c in plast: D["meta"][c]["plast"] = float(plast[c])
        else: D["meta"][c].pop("plast", None)
        if c in D.get("cellPlasticity", {}) and c in plast: D["cellPlasticity"][c] = float(plast[c])
    D["plasticityFitness"][s] = dict(cellFitness={c: round(float(f[c]), 4) for c in old_cells}, cladePlasticity=cp, cladeSize=cs)
    new_html = html[:m.start(1)] + json.dumps(D, separators=(",", ":"), allow_nan=False) + html[m.end(1):]
    bk = RUN / f"previous_{p['viewer'].name}"; assert not bk.exists(); shutil.copy2(p["viewer"], bk)
    p["viewer"].write_text(new_html)
    receipt = dict(patched_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), sample=s, viewer_sha256=sha(p["viewer"]),
                   backup=str(bk), backup_sha256=sha(bk), cells=len(old_cells), pt_cells=int(pt.reindex(list(old_cells)).notna().sum()),
                   plast_cells=len(plast), sources={str(x.relative_to(RUN)): sha(x) for x in
                   [RUN / d / f"{s}_tree_{d}.nwk" for d in p["solvers"].values()] + [RUN / "jungle_fitness" / f"{s}_lbi_fitness.csv", RUN / "phylovelo_out" / f"{s}_phylovelo_pseudotime.csv"]},
                   note="Only JZ136 replaced. Other samples' plasticityFitness predates the v8 trees (ids do not match their displayed trees) - pre-existing, untouched.")
    (RUN / "PATCH_RECEIPT.json").write_text(json.dumps(receipt, indent=2))
    return {"output_path": str(p["viewer"]), "params_used": {k: str(v) for k, v in p.items()}, "summary": receipt}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    print(json.dumps(run_pipeline()["summary"], indent=2))

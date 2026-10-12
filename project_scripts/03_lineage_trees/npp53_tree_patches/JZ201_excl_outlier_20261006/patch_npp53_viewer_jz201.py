# =============================================================================
# Npp53 tree viewer: swap in JZ201 re-run without isolated UMAP cell CACTGTCCACATGCAG (2026-10-06)
# =============================================================================
# Description : Same replacement as JZ136_intbcfix_20261004/scripts/patch_npp53_viewer_jz136.py, for JZ201-205:
#               trees (Greedy/NJ/Spectral/MaxCut/Percolation; Lineage Trajectory = Greedy), per-cell LBI 'f', 'pt'
#               (new UMAP PhyloVelo, rank-normalised within sample), 'plast' (Fitch per-cell plasticity, min clade 20),
#               clade tooltip data (viewer preorder ids). The re-run trees contain a slightly different cell set:
#               cells new to a tree get meta from the RNA metadata; dropped cells keep their unused meta entries.
# Input       : ../{greedy,nj,maxcut,spectral,percolation,jungle_fitness,phylovelo_out}/, Npp53 viewer HTML
# Output      : reports/Npp53_tree_viewer_interactive.html (+ previous_*.html backup, PATCH_RECEIPT.json)
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/hotspot/bin/python
# Validated   : 2026-10-05
# =============================================================================
import json, re, shutil, hashlib, datetime, sys
from pathlib import Path
import pandas as pd
from ete3 import Tree
RUN = Path(__file__).resolve().parents[1]
N = RUN.parents[1]
sys.path.insert(0, str(N / "lineage_tracing/JZ136_intbcfix_20261004/scripts"))
from patch_npp53_viewer_jz136 import fitch, clade_plasticity     # identical helpers
DEFAULTS = {
    "samples": ["JZ201"],
    "viewer": N / "reports/Npp53_tree_viewer_interactive.html",
    "obs": N / "data/npp53_merge/obs_metadata.csv",
    "solvers": {"Greedy": "greedy", "NJ": "nj", "Spectral": "spectral", "MaxCut": "maxcut", "Percolation": "percolation"},
    "plast_min_clade": 20, "pt_round": 4,
}
sha = lambda x: hashlib.sha256(Path(x).read_bytes()).hexdigest()

def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    html = p["viewer"].read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", html, re.S); D = json.loads(m.group(1))
    obs = pd.read_csv(p["obs"], index_col=0); obs = obs[obs.run == "MJZ019"].copy()
    obs["bc"] = obs.index.str.replace(r"^MJZ019_", "", regex=True).str.replace(r"-\d+$", "", regex=True)
    pal = D["palettes"]["state"]
    norm = lambda x: re.sub(r"^c\d+_", "", x).replace("Lum-", "Luminal-")
    name = {norm(k): k for k in pal}
    receipt = dict(patched_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), samples={})
    for s in p["samples"]:
        old_cells = set(re.findall(r"[(,]([ACGT]{16}):", D["trees"][s]["Greedy"]))
        tmpl = D["meta"][next(iter(old_cells))]
        trees = {k: Tree(open(RUN / d / f"{s}_tree_{d}.nwk").read(), format=1) for k, d in p["solvers"].items()}
        cells = set(trees["Greedy"].get_leaf_names())
        for k, t in trees.items(): assert set(t.get_leaf_names()) == cells, (s, k)
        o = obs[obs.sample_id == s].set_index("bc")
        added, dropped = cells - old_cells, old_cells - cells
        no_meta = {c for c in cells if c not in D["meta"]}   # incl. old-tree cells that never had a meta entry (pre-existing gap)
        for c in added | no_meta:
            assert c in o.index, (s, c)
            D["meta"][c] = dict(s=name[o.at[c, "putative_state"]], t=tmpl["t"], p=tmpl.get("p", ""), sample=s)
        newick = {k: t.write(format=1) for k, t in trees.items()}
        lbi = pd.read_csv(RUN / "jungle_fitness" / f"{s}_lbi_fitness.csv"); f = lbi[lbi.is_leaf].set_index("node").fitness
        assert set(f.index) == cells
        pt = pd.read_csv(RUN / "phylovelo_out" / f"{s}_phylovelo_pseudotime.csv", index_col=0).iloc[:, 0].rank(pct=True).round(p["pt_round"])
        states = {c: D["meta"][c]["s"] for c in cells}
        strip = lambda: {c: re.sub(r"^c\d+_", "", v) for c, v in states.items()}
        g = trees["Greedy"]; st = fitch(g, strip()); psum, pn = {}, {}
        for node in g.traverse("postorder"):
            if node.is_leaf() or len(node) < p["plast_min_clade"]: continue
            v = clade_plasticity(node, st)
            if v is None: continue
            for l in node.get_leaves(): psum[l.name] = psum.get(l.name, 0) + v; pn[l.name] = pn.get(l.name, 0) + 1
        plast = {c: psum[c] / pn[c] for c in psum}
        gg = Tree(newick["Greedy"], format=1); st2 = fitch(gg, strip()); cp, cs = {}, {}
        for i, n in enumerate(gg.traverse("preorder")):
            if n.is_leaf(): continue
            cs[str(i)] = len(n); v = clade_plasticity(n, st2)
            if v is not None: cp[str(i)] = round(v, 4)
        for k in p["solvers"]: D["trees"][s][k] = newick[k]
        D["trees"][s]["Lineage Trajectory"] = newick["Greedy"]
        for c in cells:
            D["meta"][c]["f"] = float(f[c])
            if c in pt.index: D["meta"][c]["pt"] = float(pt[c])
            else: D["meta"][c].pop("pt", None)
            if c in plast: D["meta"][c]["plast"] = float(plast[c])
            else: D["meta"][c].pop("plast", None)
            if c in D.get("cellPlasticity", {}) or c in added:
                if c in plast: D.setdefault("cellPlasticity", {})[c] = float(plast[c])
        D["plasticityFitness"][s] = dict(cellFitness={c: round(float(f[c]), 4) for c in cells}, cladePlasticity=cp, cladeSize=cs)
        receipt["samples"][s] = dict(cells=len(cells), added=len(added), meta_filled=len(no_meta), dropped_from_tree=len(dropped), pt_cells=int(pt.reindex(list(cells)).notna().sum()))
        print(s, receipt["samples"][s], flush=True)
    new_html = html[:m.start(1)] + json.dumps(D, separators=(",", ":"), allow_nan=False) + html[m.end(1):]
    bk = RUN / f"previous_{p['viewer'].name}"; assert not bk.exists(); shutil.copy2(p["viewer"], bk)
    p["viewer"].write_text(new_html)
    receipt.update(viewer_sha256=sha(p["viewer"]), backup=str(bk), backup_sha256=sha(bk),
                   note="Dropped cells keep their unused meta entries; CARTA panel and other samples untouched.")
    (RUN / "PATCH_RECEIPT.json").write_text(json.dumps(receipt, indent=2))
    return {"output_path": str(p["viewer"]), "summary": receipt}

if __name__ == "__main__":
    print(json.dumps(run_pipeline()["summary"]["samples"], indent=1))

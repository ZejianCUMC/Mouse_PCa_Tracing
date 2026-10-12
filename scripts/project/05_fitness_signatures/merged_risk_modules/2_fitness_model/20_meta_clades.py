"""opt5 step 20 — meta-analysis workflow (Zejian 2026-10-09: replace the vote with a random-effects meta-analysis; Greedy-tree LBI as is;
all cells incl. cycling; unsupervised). Step 1 of 3: subclones for cluster-robust standard errors (usage: 20_meta_clades.py <arm> [max_frac]).
For every tumour, the Greedy tree that produced its LBI is located and verified: by the tree sha256 recorded in the LBI run's
COMPLETE.json when available, otherwise by an identical leaf set (tree leaves == LBI leaves). The tree is then cut top-down into
subclones: any clade with more than max(10, ceil(max_frac x n_cells)) cells (default max_frac = 1/40) is split into its children, so
each tumour gets >= ~40 subclones (polytomy leaves become singletons). Cells in one subclone share most of their lineage and LBI.
Input <arm>/cells_n3000_PRIVATE.csv (all cells, Greedy LBI; JZ135 = n3000 tree). Output <arm>/clades_PRIVATE.csv (cell, sample_id,
clade, clade_size), <arm>/tree_manifest.csv. Conda env: moscot (Bio.Phylo)."""
import sys, json, glob, hashlib, math, numpy as np, pandas as pd
from pathlib import Path
from Bio import Phylo
O = Path(__file__).resolve().parent; A = O / sys.argv[1]; MAXF = float(sys.argv[2]) if len(sys.argv) > 2 else 1 / 40
MB = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
cells = pd.read_csv(A / "cells_n3000_PRIVATE.csv", index_col=0); sys.setrecursionlimit(100000)
bc = lambda s: str(s).split("_")[-1].split("-")[0]
CAND = {}
for root in [MB / "Npp53_tracing/lineage_tracing", MB / "Npp53Rb1_Tracing/lineage_tracing"]:
    for f in glob.glob(str(root / "**" / "*greedy*.nwk"), recursive=True):
        if "Cassiopeia_QC" in f or "/Trash/" in f: continue
        CAND.setdefault(Path(f).name.split("_tree_greedy")[0], []).append(f)
sha = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
rows, man = [], []
for t, src in cells.groupby("sample_id").lbi_source.first().items():
    p = Path(src) if src.startswith("/") else MB / src
    L = pd.read_csv(p); leaves_lbi = {bc(n) for n, il in zip(L.node, L.is_leaf.astype(str).str.lower()) if il == "true"}
    want = None; cj = p.parent / "COMPLETE.json"
    if cj.exists():
        try: want = json.load(open(cj)).get("tree_sha256")
        except Exception: want = None
    key = "JZ135_n3000" if t == "JZ135" else t; found, how = None, None
    for f in CAND.get(key, []):
        if want and sha(f) == want: found, how = f, "sha256"; break
    if not found:
        for f in sorted(CAND.get(key, []), key=lambda x: (str(p.parent.parent.parent) not in x, x)):
            tr = Phylo.read(f, "newick")
            if {bc(c.name) for c in tr.get_terminals()} == leaves_lbi: found, how = f, "leaf set"; break
    assert found, f"{t}: no Greedy tree matching its LBI"
    tr = Phylo.read(found, "newick"); n_leaf = len(tr.get_terminals()); cap = max(10, math.ceil(MAXF * n_leaf)); clades = []
    stack = [tr.root]
    while stack:                                    # split oversized clades; a node's direct leaf children (cells with identical
        c = stack.pop(); term = c.get_terminals()   # barcodes in a polytomy = one clone) stay together as one subclone
        if len(term) <= cap or not c.clades: clades.append([bc(x.name) for x in term]); continue
        direct = [bc(x.name) for x in c.clades if x.is_terminal()]
        if direct: clades.append(direct)
        stack.extend(x for x in c.clades if not x.is_terminal())
    lab = {b: f"{t}_c{i}" for i, cl in enumerate(clades) for b in cl}; sz = {f"{t}_c{i}": len(cl) for i, cl in enumerate(clades)}
    sub = cells[cells.sample_id == t]; cl = sub.bc16.map(lab); assert cl.notna().all(), f"{t}: cells missing from tree"
    rows.append(pd.DataFrame({"cell": sub.index, "sample_id": t, "clade": cl.values, "clade_size": cl.map(sz).values}))
    k = pd.Series([sz[c] for c in set(cl)]); man.append(dict(tumour=t, tree=found.replace(str(MB) + "/", ""), matched_by=how, tree_leaves=n_leaf, cells=len(sub), cap=cap,
                                                                    n_clades=len(set(cl)), singletons=int((k == 1).sum()), median_clade=float(k.median()), max_clade=int(k.max())))
R = pd.concat(rows); R.to_csv(A / "clades_PRIVATE.csv", index=False); M = pd.DataFrame(man); M.to_csv(A / "tree_manifest.csv", index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 90); print(M.to_string(index=False))

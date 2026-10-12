#!/usr/bin/env python
"""v10 rebuild: per-sample rings.csv (Pseudotime = v10 PhyloVelo UMAP pseudotime) +
   sidecar.h5ad taken directly from cassiopeia_trees_v8_rebuilt/phylovelo_out/{sample}_phylovelo_out.h5ad
   (which has the RECOMPUTED SHARED v10 X_umap + obs['pseudotime_umap']).
   Fallback: if a tip has no v10 pseudotime (rare), use normalized topological depth."""
from pathlib import Path
import pandas as pd, numpy as np
from ete3 import Tree
import anndata as ad
import shutil, warnings; warnings.filterwarnings('ignore')

M   = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
OUT = M/'1-mouse_based/Npp53_tracing/analysis/videos_out'
SAMPLES = ['JZ201', 'JZ202', 'JZ203', 'JZ205', 'JZ206']

def to_major(s):
    if s == 'Basal-like':   return 'Basal'
    if s == 'Partial-EMT':  return 'EMT'
    if s.startswith('Luminal-'): return 'Luminal'
    if s.startswith('Mes-'):     return 'Mesenchymal'
    return 'Other'

for sample in SAMPLES:
    print(f'\n=== {sample} ===')
    tf = M/f'cassiopeia_trees_v8_rebuilt/greedy/{sample}_tree_greedy.nwk'
    v10_h5 = M/f'cassiopeia_trees_v8_rebuilt/phylovelo_out/{sample}_phylovelo_out.h5ad'
    if not tf.exists() or not v10_h5.exists():
        print(f'  missing tree or v10 h5ad; skip'); continue
    tr = Tree(str(tf), format=1)
    # Load v10 h5ad and pull pseudotime + state
    a = ad.read_h5ad(v10_h5)
    if 'pseudotime_umap' not in a.obs.columns:
        print(f'  ! no pseudotime_umap col in v10 h5ad; skip'); continue
    pt_v10 = a.obs['pseudotime_umap'].to_dict()
    st_v10 = a.obs['putative_state'].to_dict() if 'putative_state' in a.obs.columns else {}
    # Topological-depth fallback for tips without v10 pt (rare — cells filtered out during PhyloVelo)
    for n in tr.traverse():
        if n.up is not None: n.dist = 1.0
    for l in tr.get_leaves():
        d = 0; cur = l
        while cur.up is not None:
            d += 1; cur = cur.up
        l.add_feature('_topo_depth', float(d))
    _tdmax = max(l._topo_depth for l in tr.get_leaves())
    n_v10, n_fb, rows = 0, 0, []
    for l in tr.get_leaves():
        bc = l.name
        st = st_v10.get(bc)
        if st is None or st in ('Ionocyte','Cycling'):
            continue
        pt = pt_v10.get(bc)
        if pt is None or (isinstance(pt, float) and np.isnan(pt)):
            pt = l._topo_depth / (_tdmax + 1e-9)   # fallback = normalized topological depth
            n_fb += 1
        else:
            n_v10 += 1
        rows.append(dict(TipLabel=bc, CellType=to_major(st), State=st, Pseudotime=float(pt)))
    df = pd.DataFrame(rows)
    ring_fp = OUT/'anno'/f'{sample}_rings.csv'
    df.to_csv(ring_fp, index=False)
    print(f'  {len(df)} tips → {ring_fp.name}   (v10 pt: {n_v10}, topo-depth fallback: {n_fb})')
    print(f'  CellType dist: {df["CellType"].value_counts().to_dict()}')
    # copy tree (unchanged from v9)
    (OUT/'nwk_pt').mkdir(exist_ok=True)
    dest = OUT/'nwk_pt'/f'{sample}_cassiopeia_pt0.nwk'
    shutil.copy(tf, dest)
    print(f'  copied tree → {dest.name}')
    # v10 sidecar = phylovelo_out.h5ad directly (has X_umap + obs.pseudotime_umap)
    (OUT/'sidecar').mkdir(exist_ok=True)
    sd = OUT/'sidecar'/f'{sample}_sidecar.h5ad'
    shutil.copy(v10_h5, sd)
    print(f'  sidecar (v10) → {sd.name} ({sd.stat().st_size/1e6:.1f} MB)')
print('\n[done]')

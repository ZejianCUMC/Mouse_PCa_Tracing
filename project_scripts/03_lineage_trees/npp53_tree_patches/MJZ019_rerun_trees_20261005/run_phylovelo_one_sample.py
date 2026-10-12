#!/usr/bin/env python
"""[MJZ019 rerun 2026-10-05] copy of analysis/hpc_tree_rebuild/run_phylovelo_one_sample.py (v8); ONLY input/output paths changed: input = MJZ019 re-run (min_reads_per_umi=1) per-sample tables, error intBCs removed; outputs -> lineage_tracing/MJZ019_rerun_trees_20261005/. Methods/params identical.
Run PhyloVelo on ONE sample: new Greedy tree + phylovelo h5ad.
   Outputs: {sample}_phylovelo_pseudotime.csv + {sample}_phylovelo_out.h5ad
   Uses the correct API: pv.scData + velocity_inference + calc_phylo_pseudotime
   Usage: run_phylovelo_one_sample.py <SAMPLE>"""
import os, sys, copy, warnings
warnings.filterwarnings('ignore')

LOCAL_ROOT = '/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing'
HPC_ROOT   = '/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing'
ROOT = LOCAL_ROOT if os.path.isdir(LOCAL_ROOT) else HPC_ROOT
EXCLUDE = {'Ionocyte', 'Cycling'}
SHARED_EMB_CSV = f'{ROOT}/npp53_merge/shared_embeddings_v10.csv'          # v10: RECOMPUTED UMAP+tSNE (sample.cellBC index)
SHARED_DIFFMAP_CSV = f'{ROOT}/npp53_merge/shared_diffmap_v10.csv'         # v10: RECOMPUTED diffmap
SHARED_SCVI_CSV = f'{ROOT}/npp53_merge/shared_scvi_v10.csv'               # v10: RECOMPUTED scVI (30-dim)

def main():
    sample = sys.argv[1]
    h5f     = f'{ROOT}/phylovelo_input_h5ad/{sample}.h5ad'
    tree_fp = f'/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/lineage_tracing/MJZ019_rerun_trees_20261005/greedy/{sample}_tree_greedy.nwk'
    outdir  = f'/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/lineage_tracing/MJZ019_rerun_trees_20261005/phylovelo_out'
    os.makedirs(outdir, exist_ok=True)
    print(f'[env] ROOT={ROOT}  SAMPLE={sample}', flush=True)

    import numpy as np, pandas as pd, scanpy as sc
    from Bio import Phylo
    import phylovelo as pv
    from phylovelo import pseudo_time as _pv_pt

    # Monkey-patch pv.pseudo_time.prim to handle disconnected NN graphs gracefully.
    # Original crashes with `KeyError(int)` on `nodes.remove(next)` when a visited-set
    # has no reachable unvisited neighbour. Break instead; PhyloVelo's built-in fallback
    # (calc_phylo_pseudotime line ~47) then assigns pseudotime to unvisited cells via
    # NearestNeighbors on Xdr.
    from tqdm import tqdm as _tqdm
    def _prim_patched(graph, root):
        nodes = set(graph.keys()); nodes.remove(root)
        visited = [root]; path = []
        with _tqdm(total=len(nodes)) as pbar:
            while nodes:
                distance = float("inf"); pre = None; nxt = None
                for s in visited:
                    for d in graph[s]:
                        if d in visited or s == d: continue
                        if graph[s][d] < distance:
                            distance = graph[s][d]; pre = s; nxt = d
                if nxt is None: break  # disconnected — fallback handles rest
                path.append((pre, nxt)); visited.append(nxt); nodes.remove(nxt); pbar.update(1)
        return path
    _pv_pt.prim = _prim_patched

    # Load
    a = sc.read_h5ad(h5f)
    if 'putative_state' in a.obs.columns:
        a = a[~a.obs['putative_state'].isin(EXCLUDE)].copy()
    print(f'  h5ad cells: {a.n_obs}', flush=True)

    sc.pp.filter_cells(a, min_genes=100)
    sc.pp.filter_genes(a, min_cells=10)
    a.layers['counts'] = a.X.copy()
    sc.pp.normalize_total(a); sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000)

    # --- v10 (Zejian 2026-09-02): inject RECOMPUTED shared scVI + UMAP + diffmap
    #     (computed fresh on the 18,241 no-cyc-no-ion cell set — matches §6.5b of the report). ---
    shared_emb  = pd.read_csv(SHARED_EMB_CSV,  index_col=0)
    shared_dm   = pd.read_csv(SHARED_DIFFMAP_CSV, index_col=0)
    shared_scvi = pd.read_csv(SHARED_SCVI_CSV, index_col=0)
    a.obs['_sk'] = sample + '.' + a.obs_names.astype(str)
    in_shared   = a.obs['_sk'].isin(shared_emb.index) & a.obs['_sk'].isin(shared_dm.index) & a.obs['_sk'].isin(shared_scvi.index)
    n_before    = a.n_obs
    a = a[in_shared.values].copy()
    lookup_keys = a.obs['_sk'].values
    print(f'  shared-embedding intersect (v10): {n_before} → {a.n_obs} cells', flush=True)
    a.obsm['X_umap']    = shared_emb.loc[lookup_keys, ['UMAP_1','UMAP_2']].to_numpy()
    a.obsm['X_diffmap'] = shared_dm.loc[lookup_keys, ['DC2','DC3','DC4','DC5']].to_numpy()
    a.obsm['X_scVI']    = shared_scvi.loc[lookup_keys].to_numpy()   # 30-dim RECOMPUTED scVI

    # Load tree; intersect cells that are tree tips
    base_tree = Phylo.read(tree_fp, 'newick')
    tree_tips = np.array([t.name for t in base_tree.get_terminals()])
    a = a[a.obs_names.isin(tree_tips)].copy()
    print(f'  after tree intersect: {a.n_obs} cells', flush=True)
    if a.n_obs < 30:
        print(f'[abort] too few cells ({a.n_obs})', file=sys.stderr); sys.exit(2)

    # Run PhyloVelo on UMAP (primary/default) + diffmap (secondary if valid).
    # DC0 always skipped (captures technical/library-depth artifact, not biology — Coifman 2005;
    # confirmed by Zejian 2026-09-01). Use DC1+DC2 = 2D diffmap.
    # When diffmap is degenerate (e.g. all-zero DC1/DC2, or too few cells for stable eigenvectors),
    # a *_diffmap_INVALID.flag file is written to outdir with the reason.
    def _diffmap_valid(dm_slice, min_var=1e-8):
        """Check DC1+DC2 have real spread — reject if either DC is constant/zero-variance."""
        import numpy as np
        for j in range(dm_slice.shape[1]):
            v = float(np.nanvar(dm_slice[:, j]))
            if not np.isfinite(v) or v < min_var:
                return False, f'DC{j+1} variance {v:.2e} < {min_var:.0e} (degenerate)'
        if np.isnan(dm_slice).any():
            return False, 'NaN entries in DC1/DC2'
        return True, 'OK'

    dm_slice = a.obsm['X_diffmap'][:, 1:3]
    dm_ok, dm_reason = _diffmap_valid(dm_slice)
    print(f'[diffmap validity] {dm_reason}', flush=True)
    if not dm_ok:
        flag_fp = f'{outdir}/{sample}_diffmap_INVALID.flag'
        with open(flag_fp, 'w') as f:
            f.write(f'sample={sample}\nreason={dm_reason}\ndiffmap_shape={a.obsm["X_diffmap"].shape}\n')
        print(f'  → wrote {flag_fp}; skipping diffmap branch, UMAP-only', flush=True)
        EMBEDS = {'umap': a.obsm['X_umap']}
    else:
        EMBEDS = {
            'umap':    a.obsm['X_umap'],                # PRIMARY — 2D, always robust
            'diffmap': dm_slice,                        # secondary — DC1+DC2 (DC0 skipped)
        }
    a.obs['pseudotime_umap']    = np.nan
    a.obs['pseudotime_diffmap'] = np.nan
    a.obs['pseudotime']         = np.nan            # unified default = UMAP result
    for space, emb in EMBEDS.items():
        print(f'\n[{space}] running PhyloVelo', flush=True)
        xdr = pd.DataFrame(emb, index=a.obs_names)
        sel = tree_tips[np.isin(tree_tips, xdr.index)]
        sel_data = a[sel].copy()
        print(f'  scData: {sel_data.n_obs} cells × {sel_data.n_vars} genes', flush=True)
        temp_tree = copy.deepcopy(base_tree)
        sd = pv.scData(x_normed=sel_data.to_df(), phylo_tree=temp_tree)
        depths = [len(sd.phylo_tree.get_path(sd.phylo_tree.find_any(name=i)))
                  for i in sd.x_normed.index]
        sd.Xdr = xdr.loc[sd.x_normed.index]
        try:
            pv.velocity_inference(sd, depths, cutoff=0.97, target='x_normed')
            pv.velocity_embedding(sd, target='x_normed')
            # With the prim() monkey-patch (top of main), disconnected NN graphs no longer crash;
            # PhyloVelo's built-in fallback (calc_phylo_pseudotime line ~47) assigns pseudotime to
            # unvisited cells via NearestNeighbors. Use default n_neighbors=100.
            k = min(100, sel_data.n_obs - 1)
            print(f'  [{space}] NN k={k} (prim-patched for disconnected graphs)', flush=True)
            pv.calc_phylo_pseudotime(sd, n_neighbors=k, r_sample=1)
        except Exception as e:
            print(f'  [{space}] FAILED at runtime: {e}', file=sys.stderr, flush=True)
            # Flag runtime failure per space so downstream can skip / warn
            flag_fp = f'{outdir}/{sample}_{space}_RUNTIME_FAIL.flag'
            with open(flag_fp, 'w') as f:
                f.write(f'sample={sample}\nspace={space}\nreason={e}\n')
            print(f'  → wrote {flag_fp}', flush=True)
            continue
        pt = pd.Series(sd.phylo_pseudotime, index=sd.x_normed.index, name='pseudotime')
        pt.to_csv(f'{outdir}/{sample}_phylovelo_pseudotime_{space}.csv')
        col = f'pseudotime_{space}'
        a.obs.loc[pt.index, col] = pt.values
        # copy UMAP result to the unified 'pseudotime' column (default)
        if space == 'umap':
            a.obs.loc[pt.index, 'pseudotime'] = pt.values
            pt.to_csv(f'{outdir}/{sample}_phylovelo_pseudotime.csv')   # default file
        a.obsm[f'X_phylovelo_emb_{space}'] = sd.Xdr.values
        if hasattr(sd, 'velocity_embeded'):
            a.obsm[f'X_phylovelo_velocity_{space}'] = sd.velocity_embeded
        print(f'  [{space}] pseudotime range: [{pt.min():.4f}, {pt.max():.4f}]', flush=True)

    # save combined h5ad
    a.write_h5ad(f'{outdir}/{sample}_phylovelo_out.h5ad')
    print(f'\n[saved] {sample}_phylovelo_out.h5ad + _pseudotime_{{diffmap,umap}}.csv', flush=True)

if __name__ == '__main__':
    main()

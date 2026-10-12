import sys,json,time
from pathlib import Path
import scanpy as sc
import numpy as np
out=Path(sys.argv[1]);start=time.time()
variant=sys.argv[2] if len(sys.argv)>2 else '';assert variant in ['', 'endogenous']
suffix='_endogenous' if variant else ''
assert (out/'PREPARED.json').exists() and not (out/f'HVG{suffix}_COMPLETE.json').exists()
a=sc.read_h5ad(out/'shared_counts.h5ad')
if variant:a=a[:,~a.var.gene_symbol.isin(['GFP','mCherry','Cas9','Cre','Cas9-full','Cre-full'])].copy()
# Identical feature universe and HVGs for every integration configuration.
sc.pp.filter_genes(a,min_cells=3)
sc.pp.highly_variable_genes(a,n_top_genes=3000,flavor='seurat_v3',batch_key='run_id')
a.var.to_csv(out/f'hvg_selection{suffix}.csv')
a=a[:,a.var.highly_variable].copy();assert a.n_vars==3000
a.write_h5ad(out/f'hvg_counts{suffix}.h5ad',compression='gzip')
(out/f'HVG{suffix}_COMPLETE.json').write_text(json.dumps({'cells':a.n_obs,'genes':a.n_vars,'batch_key':'run_id','flavor':'seurat_v3','variant':variant,'min_cells':3,'elapsed_s':time.time()-start,'scanpy':sc.__version__},indent=2))
print('HVG COMPLETE',a.shape,flush=True)

"""Final corrected full-gene counts for local loupeR; no external uploads."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import time,json,hashlib
import anndata as ad,numpy as np,pandas as pd,h5py
p=Path(__file__).resolve().parents[1]
out=p/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/exclude_current_c17_c19_seed42'
dest=out/'loupe_inputs';dest.mkdir(exist_ok=True)
assert not (dest/'COUNTS_COMPLETE.json').exists()
start=time.time();m=pd.read_csv(out/'metadata_PRIVATE.csv',index_col=0);labels=pd.read_csv(out/'state_labels_PRIVATE.csv',index_col=0)
assert m.index.equals(labels.index) and len(m)==40141
source=p/'analysis/min400_20260921/merged_counts.h5ad'
print('Loading corrected full counts',flush=True)
a=ad.read_h5ad(source);pos=a.obs_names.get_indexer(m.index);assert (pos>=0).all()
x=a.layers['counts'][pos].T.tocsc();x.eliminate_zeros()
assert x.shape==(33702,40141) and np.allclose(x.data,np.rint(x.data)) and np.all(x.data>=0)
assert a.var_names.is_unique and a.var.gene_id.notna().all()
assert a.var.gene_id.is_unique and {'GFP','mCherry'}<=set(a.var_names)
# Verify the prior gene-filtered tumor matrix is exactly preserved on overlap.
old=ad.read_h5ad(p/'analysis/min400_20260921/tumor_umap_20260921/tumor_full_counts.h5ad')
ix=a.var_names.get_indexer(old.var_names);assert (ix>=0).all()
assert (x[ix].T.tocsr()!=old.X[old.obs_names.get_indexer(m.index)].tocsr()).nnz==0
added=~a.var_names.isin(old.var_names)
added_detected=int(np.count_nonzero(x[added].getnnz(axis=1)))
del old
with h5py.File(dest/'final_counts.h5','w') as h:
 g=h.create_group('matrix')
 for name,arr in [('data',x.data.astype(np.int32)),('indices',x.indices.astype(np.int32)),('indptr',x.indptr.astype(np.int64)),('shape',np.array(x.shape,dtype=np.int64))]:g.create_dataset(name,data=arr,compression='gzip')
 g.create_dataset('barcodes',data=np.asarray(m.index,dtype='S'))
 f=g.create_group('features')
 for key,vals in [('id',a.var.gene_id.astype(str)),('name',a.var_names),('feature_type',['Gene Expression']*len(a.var)),('genome',['GRCm39']*len(a.var))]:f.create_dataset(key,data=np.asarray(vals,dtype='S'))
 f.create_dataset('_all_tag_keys',data=np.asarray(['genome'],dtype='S'))
export=m[['sample_id','run_id','treatment','hashtag']].copy();export['putative_state']=labels.inferred_state;export['original_r1_1_cluster']=labels.original_cluster.astype(str)
export.to_csv(dest/'metadata_PRIVATE.csv')
pd.DataFrame({'gene_id':a.var.gene_id.astype(str).values,'symbol':a.var_names}).to_csv(dest/'features.csv',index=False)
pd.DataFrame({'cell':m.index,'sum_corrected_counts':np.asarray(x.sum(axis=0)).ravel()}).to_csv(dest/'cell_count_checks_PRIVATE.csv',index=False)
digest=hashlib.sha256()
with source.open('rb') as handle:
 for block in iter(lambda:handle.read(8*1024*1024),b''):digest.update(block)
assert digest.hexdigest()=='f6c9e39f00dbcc5015ca4832e3846742bfd21a43c5aa1b866ad7070cfe246892', 'HPC source checksum mismatch'
(dest/'COUNTS_COMPLETE.json').write_text(json.dumps(dict(cells=40141,genes=a.n_vars,nnz=int(x.nnz),sum_counts=int(x.sum()),matrix_shape_genes_by_cells=list(x.shape),all_retained_genes=True,no_gene_filter=True,added_genes=int(added.sum()),added_genes_detected=added_detected,prior_25310_gene_counts_exact=True,corrected_counts=True,source=str(source),source_sha256=digest.hexdigest(),metadata_sha256=hashlib.sha256((out/'metadata_PRIVATE.csv').read_bytes()).hexdigest(),elapsed_seconds=time.time()-start),indent=2))
print(f'PASS: all {a.n_vars} genes exported for40141 exact retained cells; {added.sum()} genes restored ({added_detected} detected).',flush=True)

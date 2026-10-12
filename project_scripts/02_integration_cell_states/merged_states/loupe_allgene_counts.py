"""Lossless union of source features, without inventing ambiguous biological IDs."""
from pathlib import Path
import json,time,hashlib,sys
import h5py,numpy as np,pandas as pd
from anndata._io.specs import read_elem
from scipy import sparse
p=Path(__file__).resolve().parents[1];o=p/'analysis/merged_20260923_seed42'
c=o/'filtered_no_c17_c18_seed42/resolution_1p0/final_cleanup'
if '--merged' in sys.argv:c=o/'restored_old_c18_seed42/recluster_r1_seed42/approved_merges'
d=c/'loupe_export';d.mkdir(exist_ok=True)
assert not (d/'COUNTS_COMPLETE.json').exists()
start=time.time();meta=pd.read_csv(c/'final_metadata_PRIVATE.csv',index_col=0)
mapping=pd.read_csv(o/'feature_mapping.csv');shared=dict(zip(mapping.gene_symbol,mapping.gene_id))
sources={'Npp53':p.parent/'Npp53_tracing/data/npp53_merge/npp53_edited.h5ad','Npp53Rb1':p.parent/'Npp53Rb1_Tracing/analysis/min400_20260921/merged_counts.h5ad'}
features=[];feature_pos={};audit=[];blocks=[];order=[];source_info=[]
for model,path in sources.items():
 print('Reading full source counts:',model,flush=True)
 with h5py.File(path) as h:
  obs=read_elem(h['obs']);var=read_elem(h['var']);m=meta[meta.model.eq(model)]
  ix=obs.index.get_indexer(m.source_cell_id);assert (ix>=0).all()
  x=read_elem(h['layers/counts'])[ix].tocsr()
 assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.rint(x.data)).all()
 np.testing.assert_array_equal(np.asarray(x.sum(1)).ravel(),m.total_counts)
 cols=[]
 for j,source_feature in enumerate(var.index):
  symbol=str(source_feature) if model=='Npp53' else str(var.iloc[j].original_symbol)
  if symbol in shared:
   key=shared[symbol];display=symbol;scope='shared_unique_mapping'
  elif model=='Npp53Rb1':
   key=str(var.iloc[j].gene_id);display=f'{symbol} [Npp53Rb1:{key}]';scope='Npp53Rb1_only_source_resolved'
  else:
   key='UNRESOLVED_Npp53:'+str(source_feature);display=f'{symbol} [Npp53 unresolved]';scope='Npp53_only_unresolved_source_feature'
  if key not in feature_pos:
   feature_pos[key]=len(features);features.append({'gene_id':key,'symbol':display,'original_symbol':symbol,'availability':scope})
  cols.append(feature_pos[key]);audit.append(dict(model=model,source_feature=source_feature,source_column=j,export_column=feature_pos[key],export_id=key,export_name=display,scope=scope))
 assert len(set(cols))==len(cols)
 blocks.append((x,np.asarray(cols)));order.extend(m.index)
 source_info.append(dict(model=model,path=str(path),features=len(var),cells=len(m),counts=int(x.sum()),sha256=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()))
union=[]
for x,cols in blocks:
 z=sparse.csr_matrix((x.data.astype(np.int32),cols[x.indices],x.indptr),shape=(x.shape[0],len(features)));z.sort_indices();union.append(z)
x=sparse.vstack(union,format='csr')[pd.Index(order).get_indexer(meta.index)]
np.testing.assert_array_equal(np.asarray(x.sum(1)).ravel(),meta.total_counts)
f=pd.DataFrame(features);assert f.gene_id.is_unique and f.symbol.is_unique
assert {'GFP','mCherry'}<=set(f.symbol)
for reporter in ['GFP','mCherry']:np.testing.assert_array_equal(x[:,f.index[f.symbol.eq(reporter)][0]].toarray().ravel(),meta['raw_'+reporter])
# Exact full-shared-feature counts check, not just total counts.
with h5py.File(o/'shared_counts.h5ad') as h:
 sm=read_elem(h['obs']);sv=read_elem(h['var']);sx=read_elem(h['X'])[sm.index.get_indexer(meta.index)]
j=pd.Index(f.gene_id).get_indexer(sv.index);assert (j>=0).all() and (x[:,j]!=sx).nnz==0
x=x.T.tocsc();x.eliminate_zeros();x.sort_indices()
with h5py.File(d/'all_counts_PRIVATE.h5','w') as h:
 g=h.create_group('matrix')
 for key,values in [('data',x.data),('indices',x.indices.astype(np.int32)),('indptr',x.indptr.astype(np.int64)),('shape',np.array(x.shape,dtype=np.int64))]:g.create_dataset(key,data=values,compression='gzip')
 g.create_dataset('barcodes',data=np.asarray(meta.index,dtype='S'));fg=g.create_group('features')
 for key,values in [('id',f.gene_id),('name',f.symbol),('feature_type',['Gene Expression']*len(f)),('genome',['GRCm39']*len(f))]:fg.create_dataset(key,data=np.asarray(values,dtype='S'))
 fg.create_dataset('_all_tag_keys',data=np.asarray(['genome'],dtype='S'))
f.to_csv(d/'features.csv',index=False);pd.DataFrame(audit).to_csv(d/'source_feature_mapping.csv',index=False)
meta.to_csv(d/'metadata_PRIVATE.csv')
record=dict(status='PASS',cells=len(meta),features=len(f),all_source_features_preserved=True,no_gene_filter=True,exact_shared_counts=True,exact_source_library_totals=True,exact_reporter_counts=True,counts=int(x.sum()),nnz=int(x.nnz),sources=source_info,source_specific_features=int((f.availability!='shared_unique_mapping').sum()),unavailable_features_note='Source-specific features are structurally zero in the other cohort: not evidence of biological absence; do not use those features for cross-model DE.',seconds=time.time()-start)
(d/'COUNTS_COMPLETE.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2),flush=True)

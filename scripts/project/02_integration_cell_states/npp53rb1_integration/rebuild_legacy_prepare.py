"""Fresh corrected-count merge; no reuse of legacy embeddings or tumor gates."""
import os, json
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
import anndata as ad, scanpy as sc
from feature_name_mapping import map_corrected_features
from sample_retention import retain_sample
out=Path(os.environ['RB1_RUN'])
minimum=int(os.environ.get('RB1_MIN_VALID_CELLS','0'))
assert minimum >= 0
assert not (out/'merged_counts.h5ad').exists(), 'Refusing to overwrite an existing prepared cohort'
libs=['MJZ005','MJZ008','MJZ025','MJZ026','MJZ027','MJZ028']
objects=[]; funnel=[]
with h5py.File(out/'cellbender/MJZ005_filtered.h5') as f:
 g=f['matrix']; cb=sp.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:]))
 cbsymbols=g['features/name'][:].astype(str).tolist()
 cbids=g['features/id'][:].astype(str).tolist()
 cbcells=pd.Index(g['barcodes'][:].astype(str))
 assert cbcells.is_unique, 'Corrected cell barcodes must be unique'
raw_reference='/groups/ms3625_gp/zw2994/1-project/12-Jingbo/Nov13/counts/MJZ005_S51/outs/raw_feature_bc_matrix.h5'
with h5py.File(raw_reference) as f:
 rawids=f['matrix/features/id'][:].astype(str).tolist()
 rawsymbols=f['matrix/features/name'][:].astype(str).tolist()
exported=(out/'inputs/MJZ005/genes.txt').read_text().splitlines()
cbgenes=pd.Index(map_corrected_features(rawids,rawsymbols,exported,cbids,cbsymbols))
feature_map=pd.DataFrame({'gene_id':cbids,'original_symbol':cbsymbols,'seurat_name':cbgenes})
feature_map.to_csv(out/'MJZ005_feature_name_map.csv',index=False)
print(f'Mapped {len(cbgenes)} corrected feature rows by stable gene ID; no rows aggregated or removed.',flush=True)
for lib in libs:
 p=out/'inputs'/lib; md=pd.read_csv(p/'metadata.csv',dtype=str).set_index('cell')
 genes=pd.Index((p/'genes.txt').read_text().splitlines())
 m=sio.mmread(p/'counts.mtx').tocsc(); before=len(md)
 if lib=='MJZ005':
  good=md.index.isin(cbcells); md=md.loc[good].copy()
  m=cb[:,cbcells.get_indexer(md.index)]; genes=cbgenes
  md['count_source']='CellBender_fresh_MJZ005'
 assert m.shape==(len(genes),len(md))
 a=ad.AnnData(m.T.tocsr(),obs=md.copy(),var=pd.DataFrame(index=genes))
 a.var['mt']=a.var_names.str.startswith('mt-')
 sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
 # Apply the reference's basic RNA QC consistently on corrected full counts.
 a=a[(a.obs.n_genes_by_counts>200)&(a.obs.pct_counts_mt<10)].copy()
 a.obs['source_barcode']=a.obs_names; a.obs_names=lib+'_'+a.obs_names
 a.obs['batch']=lib; a.obs['cohort']='legacy' if lib in ['MJZ005','MJZ008'] else 'new'
 funnel.append(dict(library=lib,source_singlets=before,corrected_matrix_overlap=len(md),post_corrected_RNA_QC=a.n_obs))
 objects.append(a)
a=ad.concat(objects,join='inner',merge='same')
assert a.obs_names.is_unique and not a.obs.sample_id.isin(['JZ136','JZ282']).any()
assert set(['JZ135','JZ137','JZ138']).issubset(a.obs.sample_id)
assert (a.obs.loc[a.obs.sample_id.isin(['JZ135','JZ137','JZ138']),'treatment']=='Intact').all()
assert a.n_vars>10000
feature_map=feature_map.set_index('seurat_name')
assert a.var_names.isin(feature_map.index).all()
a.var['gene_id']=feature_map.loc[a.var_names,'gene_id'].to_numpy()
a.var['original_symbol']=feature_map.loc[a.var_names,'original_symbol'].to_numpy()
a.layers['counts']=a.X.copy()
for gene in ['GFP','mCherry']:
 assert gene in a.var_names
 a.obs['raw_'+gene]=np.asarray(a[:,gene].X.sum(axis=1)).ravel()
a.obs['reporter_positive']=(a.obs.raw_GFP>0)|(a.obs.raw_mCherry>0)
a.obs['reporter_ge2']=(a.obs.raw_GFP>=2)|(a.obs.raw_mCherry>=2)
sample_counts=a.obs.groupby(['library','sample_id','hashtag','treatment'],observed=True).size().rename('cells').reset_index()
sample_counts['included_downstream']=sample_counts.cells.map(lambda n:retain_sample(n,minimum))
sample_counts.to_csv(out/'sample_retention.csv',index=False)
sample_counts.to_csv(out/'sample_counts_before_threshold.csv',index=False)
excluded=sample_counts.loc[~sample_counts.included_downstream,'sample_id'].tolist()
before_threshold=a.n_obs
a=a[~a.obs.sample_id.isin(excluded)].copy()
assert a.n_obs>0 and not a.obs.sample_id.isin(excluded).any()
if minimum==400:
 assert set(excluded)=={'JZ257','JZ275','JZ278','JZ281'}
 assert a.n_obs==56686 and a.obs.sample_id.nunique()==22
(out/'sample_retention_summary.json').write_text(json.dumps({'minimum_valid_cells':minimum,'rule':'retain cells >= threshold per sample, measured after corrected-RNA QC','before_cells':before_threshold,'retained_cells':a.n_obs,'retained_samples':int(a.obs.sample_id.nunique()),'excluded_samples':excluded,'excluded_cells':before_threshold-a.n_obs},indent=2))
tg=['GFP','mCherry','Cas9','Cas9_full','Cas9-full','Cre','Cre_full','Cre-full']
eligible=a[:,~a.var_names.isin(tg)].copy()
sc.pp.highly_variable_genes(eligible,flavor='seurat_v3',n_top_genes=2000,batch_key='batch',layer='counts')
hvg=eligible.var_names[eligible.var.highly_variable]
a.var['highly_variable']=a.var_names.isin(hvg)
a.write(out/'merged_counts.h5ad')
sio.mmwrite(out/'merged_counts.mtx',a.X.T)
(out/'genes.txt').write_text('\n'.join(a.var_names)+'\n')
(out/'hvg.txt').write_text('\n'.join(hvg)+'\n')
a.obs.to_csv(out/'metadata.csv')
pd.DataFrame(funnel).to_csv(out/'QC_funnel.csv',index=False)
a.obs.groupby(['library','sample_id','hashtag','treatment'],observed=True).size().rename('cells').reset_index().to_csv(out/'sample_counts.csv',index=False)
(out/'cohort_summary.json').write_text(json.dumps({'cells':a.n_obs,'genes':a.n_vars,'mice':int(a.obs.sample_id.nunique()),'sample_minimum_applied':minimum>0,'minimum_valid_cells':minimum,'excluded_samples':excluded,'input_singlets':sum(x['source_singlets'] for x in funnel)},indent=2))
print(a); print(pd.DataFrame(funnel).to_string(index=False))

"""Freeze retained cohorts and align only uniquely mapped shared genes."""
import argparse, hashlib, json, platform, time
from pathlib import Path
import anndata as ad
from anndata._io.specs import read_elem
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

p=argparse.ArgumentParser();p.add_argument('--parent',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
args=p.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
assert not (out/'PREPARED.json').exists(), 'Refusing duplicate preparation'
start=time.time()
n=args.parent/'Npp53_tracing';r=args.parent/'Npp53Rb1_Tracing'
final=r/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/exclude_current_c17_c19_seed42'
sources=[n/'data/npp53_merge/npp53_edited.h5ad',r/'analysis/min400_20260921/merged_counts.h5ad']
cohorts=[pd.read_csv(n/'data/npp53_merge/obs_no_cycling_no_ionocyte.csv',index_col=0),pd.read_csv(final/'metadata_PRIVATE.csv',index_col=0)]
labels=pd.read_csv(final/'state_labels_PRIVATE.csv',index_col=0)
arrays=[];variables=[];observations=[];manifest=[]
for model,path,cohort,expected in zip(['Npp53','Npp53Rb1'],sources,cohorts,[18241,40141]):
 print('Reading',model,path,flush=True)
 with h5py.File(path) as f:
  obs=read_elem(f['obs']);var=read_elem(f['var'])
  assert len(cohort)==expected and cohort.index.is_unique and obs.index.is_unique
  ix=obs.index.get_indexer(cohort.index);assert (ix>=0).all()
  counts=read_elem(f['layers/counts'])[ix,:].tocsr()
 assert np.isfinite(counts.data).all() and (counts.data>=0).all()
 assert np.allclose(counts.data,np.rint(counts.data),atol=1e-6,rtol=0), 'Noninteger counts'
 counts=counts.astype(np.float32)
 meta=pd.DataFrame(index=cohort.index)
 meta['source_cell_id']=cohort.index;meta['model']=model
 meta['sample_id']=cohort.sample_id.astype(str)
 meta['run_id']=cohort['run' if model=='Npp53' else 'run_id'].astype(str)
 meta['treatment']=cohort.treatment.replace({'Castration + Enzalutamide':'Castration+ENZA','Castration + ENZA':'Castration+ENZA'})
 assert set(meta.treatment)<=set(['Intact','Castration','Castration+ENZA'])
 meta['original_state']=cohort.putative_state if model=='Npp53' else labels.loc[cohort.index,'inferred_state']
 meta['total_counts']=np.asarray(counts.sum(1)).ravel()
 meta['n_genes_by_counts']=np.asarray((counts>0).sum(1)).ravel()
 sym=var.index.to_series() if model=='Npp53' else var.original_symbol
 for reporter in ['GFP','mCherry']:
  cols=np.flatnonzero(sym.values==reporter);assert len(cols)==1, (model,reporter,cols)
  meta['raw_'+reporter]=np.asarray(counts[:,cols].sum(1)).ravel()
 mt=np.asarray(sym.str.lower().str.startswith('mt-'))
 meta['pct_counts_mt']=100*np.asarray(counts[:,mt].sum(1)).ravel()/meta.total_counts.values
 meta.index.name='cell'
 arrays.append(counts);variables.append(var);observations.append(meta)
 manifest.append({'model':model,'path':str(path.resolve()),'size':path.stat().st_size,'mtime_ns':path.stat().st_mtime_ns,'cells':len(meta),'features':len(var),'counts_layer':'counts','integer_nonnegative_finite':True})
 print('Validated',model,counts.shape,flush=True)

nv,rv=variables
ambiguous=set(rv.loc[rv.original_symbol.duplicated(False),'original_symbol'])
mapping=rv.loc[~rv.original_symbol.duplicated(False)].set_index('original_symbol')['gene_id']
shared=[g for g in nv.index if g in mapping.index];shared_set=set(shared)
assert len(shared)>33000
var=pd.DataFrame({'gene_symbol':shared,'gene_id':mapping.loc[shared].values},index=mapping.loc[shared].values)
var.index.name='feature_id';assert var.index.is_unique
var['mapping_provenance']='Npp53 symbol uniquely matches Npp53Rb1 reference ID; original Npp53 IDs absent'
excluded=[]
for model,v in zip(['Npp53','Npp53Rb1'],variables):
 symbols=v.index if model=='Npp53' else v.original_symbol
 for feature,symbol in zip(v.index,symbols):
  if symbol not in shared_set:excluded.append({'model':model,'source_feature':feature,'symbol':symbol,'reason':'ambiguous multi-ID symbol' if symbol in ambiguous else 'not shared/construct naming differs'})
pd.DataFrame(excluded).to_csv(out/'excluded_integration_features.csv',index=False)
rsymbols=rv.original_symbol.tolist()
ridx={s:i for i,s in enumerate(rsymbols) if s not in ambiguous}
X=sparse.vstack([arrays[0][:,nv.index.get_indexer(shared)],arrays[1][:,[ridx[g] for g in shared]]],format='csr')
obs=pd.concat(observations);assert len(obs)==58382 and obs.index.is_unique, 'Unexpected cohort overlap'
assert obs.groupby('sample_id').model.nunique().max()==1
obs['shared_counts']=np.asarray(X.sum(1)).ravel()
a=ad.AnnData(X=X,obs=obs,var=var)
a.uns['counts_description']='X is unnormalized, integer-valued source counts. No log transform.'
a.uns['source_manifest_json']=json.dumps(manifest)
print('Writing shared full-count object',a.shape,flush=True)
a.write_h5ad(out/'shared_counts.h5ad',compression='gzip')
obs.to_csv(out/'metadata_PRIVATE.csv')
obs.groupby(['model','sample_id','run_id','treatment']).size().rename('n_cells').reset_index().to_csv(out/'cohort_summary.csv',index=False)
var.to_csv(out/'feature_mapping.csv')
record={'host':platform.node(),'elapsed_s':time.time()-start,'n_cells':len(obs),'n_shared_genes':len(var),'n_samples':obs.sample_id.nunique(),'n_runs':obs.run_id.nunique(),'counts_by_model':obs.model.value_counts().to_dict(),'sources':manifest,'excluded_ambiguous_symbols':sorted(ambiguous),'seed':42,'status':'PREPARED','counts_sha256':hashlib.file_digest(open(out/'shared_counts.h5ad','rb'),'sha256').hexdigest()}
(out/'PREPARED.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record,indent=2),flush=True)

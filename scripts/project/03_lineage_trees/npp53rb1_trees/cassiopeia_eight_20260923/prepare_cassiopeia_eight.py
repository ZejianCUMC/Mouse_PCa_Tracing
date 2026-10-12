"""Freeze the eight reviewed RNA/TSL intersections; never process held runs."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json,hashlib,time,platform
from pathlib import Path
import pandas as pd

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    project=Path(__file__).resolve().parents[1];mouse=project.parent
    final=project/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/exclude_current_c17_c19_seed42'
    out=project/'lineage_tracing/cassiopeia_eight_20260923_seed42'
    assert not (out/'PREPARED.json').exists(), 'Frozen inputs already exist'
    out.mkdir(parents=True,exist_ok=True);(out/'inputs').mkdir(exist_ok=True)
    meta=pd.read_csv(final/'metadata_PRIVATE.csv',index_col=0)
    labels=pd.read_csv(final/'state_labels_PRIVATE.csv',index_col=0)
    assert meta.index.equals(labels.index) and len(meta)==40141
    meta['putative_state']=labels.inferred_state
    meta['cellBC']=meta.source_barcode.str.replace(r'-[0-9]+$','',regex=True)
    expected=pd.read_csv(project/'metadata/tsl_availability.csv')
    assert expected.sample_id.tolist()==['JZ135','JZ137','JZ138','JZ256','JZ259','JZ261','JZ262','JZ263']
    sources={
        'JZ135':Path('/media/jrchristin/BackUp8T2/project/12-Jingbo/MJZ005/MJZ005_allele_table_min_umi_per_cell_10_with_annotation.csv'),
        'JZ137':mouse/'MJZ008/MJZ008_allele_table_min_umi_per_cell_10_JZ137.csv',
        'JZ138':mouse/'MJZ008/MJZ008_allele_table_min_umi_per_cell_10_JZ138.csv',
    }
    shared=mouse/'Npp53_tracing/lineage_tracing/allele_tables/per_library/MJZ025_allele_table.csv'
    columns=['cellBC','intBC','allele','r1','r2','r3','lineageGrp','UMI','readCount']
    cache={};records=[];started=time.time()
    for r in expected.itertuples(index=False):
        source=sources.get(r.sample_id,shared)
        if source not in cache:
            print('Read',source,flush=True)
            # Only one full library table is held at a time.
            cache={source:pd.read_csv(source,usecols=columns,dtype={'cellBC':str,'intBC':str})}
        a=cache[source];m=meta[meta.sample_id==r.sample_id].copy()
        assert len(m)==r.final_tumor_cells and set(m.run_id)=={r.run_id} and m.cellBC.is_unique
        assert r.run_id in ['MJZ005','MJZ008','MJZ025']
        z=a[a.cellBC.isin(m.cellBC)].copy()
        assert z.cellBC.nunique()==r.tsl_matched_cells
        assert z[['cellBC','intBC','r1','r2','r3']].notna().all().all()
        assert (z.UMI>0).all() and (z.readCount>0).all()
        conflicts=z.groupby(['cellBC','intBC'],sort=False).allele.nunique()
        duplicate_rows=int(z.duplicated(columns).sum())
        assert duplicate_rows==0,'Unexpected exact duplicate source rows: inspect before solving'
        matched=m[m.cellBC.isin(z.cellBC)].copy()
        mapping=dict(zip(matched.cellBC,matched.index))
        # Prefix every barcode with run_id; never join cells across libraries.
        z['cellBC']=z.cellBC.map(mapping)
        z=z.sort_values(['cellBC','intBC','r1','r2','r3','UMI','readCount'],kind='stable').reset_index(drop=True)
        dest=out/'inputs'/r.sample_id;dest.mkdir(exist_ok=False)
        z.to_csv(dest/'allele_table.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        matched=matched.sort_index()
        matched[['sample_id','run_id','treatment','hashtag','source_barcode','putative_state']].to_csv(dest/'cell_metadata.csv',index_label='cell')
        m.assign(has_TSL=m.cellBC.isin(mapping))[['sample_id','run_id','source_barcode','has_TSL']].to_csv(dest/'rna_tsl_overlap.csv',index_label='cell')
        qc=dict(sample_id=r.sample_id,run_id=r.run_id,final_tumor_cells=int(len(m)),matched_cells=int(len(matched)),
            allele_rows=int(len(z)),cell_intbc_pairs=int(len(conflicts)),multi_allele_pairs=int((conflicts>1).sum()),
            median_intbcs=float(z.groupby('cellBC').intBC.nunique().median()),median_umis=float(z.groupby('cellBC').UMI.sum().median()),
            source=str(source),source_sha256=sha(source),allele_table_sha256=sha(dest/'allele_table.csv.gz'),
            metadata_sha256=sha(dest/'cell_metadata.csv'),ambiguity_policy='Preserve all distinct allele states as Cassiopeia tuples; no arbitrary dominant-allele call',
            no_new_cell_exclusions=True)
        (dest/'INPUT_QC.json').write_text(json.dumps(qc,indent=2));records.append(qc)
        print('Prepared',r.sample_id,len(matched),'cells;',qc['multi_allele_pairs'],'multi-allele pairs',flush=True)
    pd.DataFrame([{k:v for k,v in r.items() if k not in ['source','ambiguity_policy']} for r in records]).to_csv(out/'input_summary.csv',index=False)
    manifest=dict(samples=[r['sample_id'] for r in records],seed=42,solver='VanillaGreedySolver',allele_rep_thresh=.99,
        collapse_duplicates=True,mutation_priors=None,collapse_mutationless_edges=False,
        held_runs=['MJZ026','MJZ027','MJZ028'],matched_cells=sum(r['matched_cells'] for r in records),
        records=records,source_metadata_sha256=sha(final/'metadata_PRIVATE.csv'),source_state_labels_sha256=sha(final/'state_labels_PRIVATE.csv'),
        host=platform.node(),elapsed_seconds=time.time()-started)
    (out/'PREPARED.json').write_text(json.dumps(manifest,indent=2))
    print('PASS: exact eight-sample RNA/TSL inputs frozen; no held runs.',flush=True)

if __name__=='__main__':main()

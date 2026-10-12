"""Freeze RNA/TSL intersections for the MJZ026-028 Np53Rb1 mice (hold lifted by Zejian 2026-09-23).
Adapted from scripts/prepare_cassiopeia_eight.py: identical join/format/QC logic; differences are
(1) samples = every final-cohort mouse in runs MJZ026/027/028, (2) each run maps to ITS OWN library
allele table (no MJZ025 fallback), (3) matched-cell counts are recorded as observed, not pre-asserted."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json,hashlib,time,platform
from pathlib import Path
import pandas as pd

RUNS=['MJZ026','MJZ027','MJZ028']

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    out=Path(__file__).resolve().parents[1]
    project=out.parents[1]
    final=project/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/exclude_current_c17_c19_seed42'
    tables=project/'lineage_tracing/allele_tables_MJZ026_028_20260925'
    sources={r:tables/f'{r}_allele_table.csv' for r in RUNS}
    assert not (out/'PREPARED.json').exists(), 'Frozen inputs already exist'
    (out/'inputs').mkdir(exist_ok=True)
    meta=pd.read_csv(final/'metadata_PRIVATE.csv',index_col=0)
    labels=pd.read_csv(final/'state_labels_PRIVATE.csv',index_col=0)
    assert meta.index.equals(labels.index) and len(meta)==40141
    meta['putative_state']=labels.inferred_state
    meta['cellBC']=meta.source_barcode.str.replace(r'-[0-9]+$','',regex=True)
    samples=(meta[meta.run_id.isin(RUNS)].groupby('sample_id').run_id.first().sort_index())
    columns=['cellBC','intBC','allele','r1','r2','r3','lineageGrp','UMI','readCount']
    cache={};records=[];started=time.time()
    for sample_id,run_id in samples.items():
        source=sources[run_id]
        if source not in cache:
            print('Read',source,flush=True)
            full=pd.read_csv(source,usecols=columns,keep_default_na=False,
                dtype={c:str for c in ['cellBC','intBC','allele','r1','r2','r3','lineageGrp']})
            # Incomplete reads: empty allele or any uncovered cut site (''). Uninterrupted Cassiopeia runs keep
            # these rows; checkpoint-resumed runs (e.g. MJZ025) lose them. Zejian 2026-09-25: drop them, so the
            # new mice match MJZ025's input rule. Counts recorded per sample below.
            full['incomplete']=(full[['allele','r1','r2','r3']]=='').any(axis=1)
            cache={source:full}
        a=cache[source];m=meta[meta.sample_id==sample_id].copy()
        assert set(m.run_id)=={run_id} and m.cellBC.is_unique
        z_all=a[a.cellBC.isin(m.cellBC)]
        n_incomplete=int(z_all.incomplete.sum());cells_all=z_all.cellBC.nunique()
        z=z_all[~z_all.incomplete].drop(columns='incomplete').copy()
        z['lineageGrp']=pd.to_numeric(z.lineageGrp)
        assert z[['cellBC','intBC','r1','r2','r3']].notna().all().all(), f'{sample_id}: missing cellBC/intBC/site values'
        assert (z.UMI>0).all() and (z.readCount>0).all()
        conflicts=z.groupby(['cellBC','intBC'],sort=False).allele.nunique()
        duplicate_rows=int(z.duplicated(columns).sum())
        assert duplicate_rows==0,'Unexpected exact duplicate source rows: inspect before solving'
        matched=m[m.cellBC.isin(z.cellBC)].copy()
        mapping=dict(zip(matched.cellBC,matched.index))
        z['cellBC']=z.cellBC.map(mapping)          # run_id-prefixed cell ids; never join across libraries
        z=z.sort_values(['cellBC','intBC','r1','r2','r3','UMI','readCount'],kind='stable').reset_index(drop=True)
        dest=out/'inputs'/sample_id;dest.mkdir(exist_ok=False)
        z.to_csv(dest/'allele_table.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        matched=matched.sort_index()
        matched[['sample_id','run_id','treatment','hashtag','source_barcode','putative_state']].to_csv(dest/'cell_metadata.csv',index_label='cell')
        m.assign(has_TSL=m.cellBC.isin(mapping))[['sample_id','run_id','source_barcode','has_TSL']].to_csv(dest/'rna_tsl_overlap.csv',index_label='cell')
        qc=dict(sample_id=sample_id,run_id=run_id,final_tumor_cells=int(len(m)),matched_cells=int(len(matched)),
            pct_matched=round(100*len(matched)/len(m),1),
            incomplete_rows_dropped=n_incomplete,cells_before_incomplete_drop=int(cells_all),
            incomplete_row_policy='drop rows with empty allele or any empty r1/r2/r3 (match MJZ025 input rule)',
            allele_rows=int(len(z)),cell_intbc_pairs=int(len(conflicts)),multi_allele_pairs=int((conflicts>1).sum()),
            median_intbcs=float(z.groupby('cellBC').intBC.nunique().median()),median_umis=float(z.groupby('cellBC').UMI.sum().median()),
            source=str(source),source_sha256=sha(source),allele_table_sha256=sha(dest/'allele_table.csv.gz'),
            metadata_sha256=sha(dest/'cell_metadata.csv'),ambiguity_policy='Preserve all distinct allele states as Cassiopeia tuples; no arbitrary dominant-allele call',
            no_new_cell_exclusions=True)
        (dest/'INPUT_QC.json').write_text(json.dumps(qc,indent=2));records.append(qc)
        print('Prepared',sample_id,run_id,len(matched),'/',len(m),'cells;',qc['multi_allele_pairs'],'multi-allele pairs',flush=True)
    pd.DataFrame([{k:v for k,v in r.items() if k not in ['source','ambiguity_policy']} for r in records]).to_csv(out/'input_summary.csv',index=False)
    manifest=dict(samples=[r['sample_id'] for r in records],seed=42,solver='VanillaGreedySolver',allele_rep_thresh=.99,
        collapse_duplicates=True,mutation_priors=None,collapse_mutationless_edges=False,
        held_runs=[],hold_lifted='Zejian 2026-09-23: build trees for MJZ026-028 after allele-table QC',
        matched_cells=sum(r['matched_cells'] for r in records),
        records=records,source_metadata_sha256=sha(final/'metadata_PRIVATE.csv'),source_state_labels_sha256=sha(final/'state_labels_PRIVATE.csv'),
        host=platform.node(),elapsed_seconds=time.time()-started)
    (out/'PREPARED.json').write_text(json.dumps(manifest,indent=2))
    print(f'PASS: {len(records)} MJZ026-028 RNA/TSL inputs frozen.',flush=True)

if __name__=='__main__':main()

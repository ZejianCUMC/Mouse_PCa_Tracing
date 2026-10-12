"""Reconcile corrected CHOIR/AMI with exact matching completed sc-SHC tests.

Preserve rejected reference artifacts. No HTML publication or state naming.
Run only after both persistent services finish and manual correction review.
"""
import datetime,hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
import pandas as pd

project=Path(__file__).resolve().parents[1]
out=project/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42'
fix=out/'choir_label_fix_seed42'
for unit in ['rb1-filtered-states-20260922-resume1.service','rb1-choir-label-fix-20260922.service']:
    query=subprocess.run(['systemctl','--user','show',unit,'-p','MainPID','-p','ActiveState','-p','LoadState','-p','ExecMainStatus'],capture_output=True,text=True)
    state=dict(line.split('=',1) for line in query.stdout.splitlines() if '=' in line)
    assert state.get('MainPID','0')=='0' and state.get('ActiveState') in ['inactive','failed'],(unit,state,query.stderr)
    assert state.get('LoadState') in ['loaded','not-found'],(unit,state)
    if 'label-fix' in unit and state.get('LoadState')=='loaded':assert state.get('ExecMainStatus')=='0'
assert (fix/'CHOIR_LABEL_FIX_VERIFIED.json').exists()
assert json.loads((fix/'CHOIR_LABEL_AUDIT.json').read_text())['valid']
assert json.loads((fix/'TREE_COMPARISON.json').read_text())['same_initial_tree']
assert not (out/'ADOPTED_CORRECTED_REFERENCE.json').exists()
assert not (out/'publication_acceptance.json').exists()
original_manifest=json.loads((out/'RESUME_STARTED.json').read_text())
for name,digest in original_manifest['retained_input_sha256'].items():
    assert hashlib.sha256((out/name).read_bytes()).hexdigest()==digest,name
plan=json.loads((fix/'CORRECTION_PLAN.json').read_text())
for name,digest in plan['input_sha256'].items():
    assert hashlib.sha256((fix/name).read_bytes()).hexdigest()==digest,name
    assert hashlib.sha256((out/name).read_bytes()).hexdigest()==digest,name
for name,digest in plan['script_sha256'].items():
    assert hashlib.sha256((project/'scripts'/name).read_bytes()).hexdigest()==digest,name
sc_manifest=json.loads((out/'scSHC_COMPLETE.json').read_text())
assert sc_manifest['cells']==42806 and sc_manifest['seed']==42
assert sc_manifest['batch']=='run_id' and sc_manifest['PCs']==30 and sc_manifest['alpha']==.05
assert sc_manifest['RNG_kind'][0]=="L'Ecuyer-CMRG"
rank=pd.read_csv(fix/'AMI_ranking.csv')
old_summary=pd.read_csv(out/'scSHC_summary.csv')
required=rank.head(5).resolution.tolist()
missing=[r for r in required if not np.isclose(old_summary.resolution,r).any()]
assert not missing,f'Newly selected candidates need genuine sc-SHC tests before adoption: {missing}'
assign=pd.read_csv(out/'resolution_assignments_PRIVATE.csv',index_col=0,dtype=str)
rows=[];ledger=[]
for item in rank.head(5).itertuples():
    filename=f'scSHC_r{item.resolution:.1f}_labels_PRIVATE.csv'
    labels=pd.read_csv(out/filename,index_col=0,dtype=str)
    expected=assign[f'leiden_r{item.resolution:.1f}_seed42']
    assert labels.index.equals(assign.index) and len(labels)==42806 and not labels.isna().any().any()
    assert labels.leiden_cluster.equals(expected.rename('leiden_cluster'))
    assert labels.groupby('leiden_cluster').scSHC_cluster.nunique().eq(1).all()
    summary=old_summary.loc[np.isclose(old_summary.resolution,item.resolution)].iloc[0].to_dict()
    assert summary['input_clusters']==item.clusters==expected.nunique()
    assert summary['defended_clusters']==labels.scSHC_cluster.nunique()
    rows.append(summary)
    ledger.append(dict(resolution=item.resolution,label_sha256=hashlib.sha256((out/filename).read_bytes()).hexdigest(),
                       reason='Exact same retained counts, barcode order, seed42 partition, batch variable, RNG and test settings; CHOIR affects selection only'))
summary=pd.DataFrame(rows)
summary[['input_clusters','defended_clusters']]=summary[['input_clusters','defended_clusters']].astype(int)
evidence=rank.head(5)[['rank','resolution','AMI','clusters']].merge(summary,on='resolution',validate='one_to_one')
evidence=evidence[['rank','resolution','AMI','input_clusters','defended_clusters','retained_fraction']]

archive=out/'rejected_choir_reference_20260922'
archive.mkdir(exist_ok=False)
old_files=['choir_HVG_object_PRIVATE.RDS','choir_labels_PRIVATE.csv','choir_cluster_sizes.csv',
    'CHOIR_COMPLETE.json','CHOIR_LABEL_AUDIT.json','choir_invalid_internal_labels_PRIVATE.csv',
    'AMI_ranking.csv','RANKED.json','scSHC_candidates.csv','scSHC_summary.csv','scSHC_COMPLETE.json',
    'top5_evidence.csv','COMPLETE.json','COMPUTATION_VERIFIED.json','PLOTS_COMPLETE.json',
    'state_resolution_review.png','state_seed_sensitivity.png']
for name in old_files:
    if (out/name).exists():shutil.move(str(out/name),str(archive/name))
for name in ['choir_HVG_object_PRIVATE.RDS','choir_labels_PRIVATE.csv','choir_cluster_sizes.csv',
             'CHOIR_COMPLETE.json','CHOIR_LABEL_AUDIT.json','choir_invalid_internal_labels_PRIVATE.csv',
             'CHOIR_LABEL_FIX_VERIFIED.json','AMI_ranking.csv','RANKED.json','scSHC_candidates.csv']:
    shutil.copy2(fix/name,out/name)
summary.to_csv(out/'scSHC_summary.csv',index=False)
evidence.to_csv(out/'top5_evidence.csv',index=False)
sc_manifest.update(candidates=required,reused_candidate_tests=ledger,
                   selection_reference='Corrected and internally audited CHOIR; original candidate tests unchanged')
(out/'scSHC_COMPLETE.json').write_text(json.dumps(sc_manifest,indent=2))
record=dict(cells=42806,seed=42,completed=datetime.datetime.now().astimezone().isoformat(),
    correction='Signed-level CHOIR label parser; exactly matching initial tree verified',
    source='choir_label_fix_seed42',rejected_reference_archive=archive.name,
    candidate_test_reuse=ledger,status='Corrected computations assembled; validation required')
(out/'COMPLETE.json').write_text(json.dumps(record,indent=2))
(out/'ADOPTION_STARTED.json').write_text(json.dumps(record,indent=2))
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='4'
os.environ['MPLCONFIGDIR']=tempfile.mkdtemp(prefix='rb1-corrected-plot-')
for script,arguments,logname in [
    ('plot_top5_evidence.py',[str(out)],'corrected_plots.log'),
    ('test_filtered_states.py',['--complete'],'corrected_validation.log')]:
    with (out/logname).open('x') as handle:
        subprocess.run([sys.executable,str(project/'scripts'/script),*arguments],stdout=handle,stderr=subprocess.STDOUT,check=True)
record['status']='Corrected computational validation passed; manual interpretation/publication pending'
(out/'COMPUTATION_VERIFIED.json').write_text(json.dumps(record,indent=2))
(out/'ADOPTED_CORRECTED_REFERENCE.json').write_text(json.dumps(record,indent=2))
(out/'STATUS.json').write_text(json.dumps({**record,'stage':record['status'],'updated':record['completed']},indent=2))
print(evidence.to_string(index=False))
print('PASS: rejected reference preserved; corrected CHOIR and AMI adopted; exact matching sc-SHC tests reconciled; no HTML published.')

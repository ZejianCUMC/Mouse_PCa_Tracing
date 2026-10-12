"""Submit the error-intBC rebuild (Zejian 2026-10-04) once; no automatic retries.
Same methods/resources as cassiopeia_eight_20260923 + tree_viewer_20260923 (JZ137/JZ138/JZ135 subsets):
  Greedy array (5)  -> verify          [1CPU/32GB/6h, max4]
  variants (afterok greedy): JZ137/JZ138 x NJ/Spectral/MaxCut/Percolation + JZ135_n* NJ only   [1CPU/16GB/6h, max4]
  LBI      (afterok greedy): 5        [1CPU/8GB/2h, max4]
  PhyloVelo(afterok greedy): 5        [1CPU/32GB/8h, max2]"""
from pathlib import Path
import os,json,subprocess,shlex,time,hashlib
ROOT=Path(__file__).resolve().parents[1]
assert str(ROOT).startswith('/groups/')
assert not (ROOT/'SUBMISSION.json').exists() and not (ROOT/'cassiopeia/results').exists() and not (ROOT/'results').exists()
lock=os.open(ROOT/'SUBMISSION_LOCK',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(lock)
CAS=ROOT/'cassiopeia'
cas_manifest=json.loads((CAS/'PREPARED.json').read_text())
for s_id,r in zip(cas_manifest['samples'],cas_manifest['records']):   # JZ135_n* records carry source sample_id 'JZ135'
    p=CAS/'inputs'/s_id
    for name,key in [('allele_table.csv.gz','allele_table_sha256'),('cell_metadata.csv','metadata_sha256')]:
        assert hashlib.sha256((p/name).read_bytes()).hexdigest()==r[key], (p,name)
samples=json.loads((ROOT/'PREPARED.json').read_text())['samples']
assert samples==cas_manifest['samples']==['JZ137','JZ138','JZ135_n1000','JZ135_n2000','JZ135_n3000']
active=subprocess.check_output(['squeue','-u',os.environ['USER'],'-h','-o','%j'],text=True)
assert 'rb1_fix_' not in active
(ROOT/'logs').mkdir(exist_ok=True)
env='/groups/ms3625_gp/zw2994/mamba/envs/'
receipt=dict(submitted_unix=time.time(),jobs=[])
def save():(ROOT/'SUBMISSION.json').write_text(json.dumps(receipt,indent=2))
def submit(name,tasks,mem,hours,parallel,dep=None):
    table=ROOT/(name+'_tasks.json');table.write_text(json.dumps(tasks,indent=2))
    command=f"export PYTHONHASHSEED=42 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MPLCONFIGDIR={ROOT}/logs/mpl_$SLURM_ARRAY_TASK_ID; {env}cassiopeia_env/bin/python {shlex.quote(str(ROOT/'scripts/dispatch_tree_viewer.py'))} {shlex.quote(str(ROOT))} {shlex.quote(str(table))}"
    cmd=['sbatch','--parsable','--job-name=rb1_fix_'+name,'--cpus-per-task=1',f'--mem={mem}G',f'--time={hours}:00:00',f'--array=0-{len(tasks)-1}%{parallel}',f'--output={ROOT}/logs/{name}_%A_%a.out',f'--error={ROOT}/logs/{name}_%A_%a.err']
    if dep:cmd.append(f'--dependency=afterok:{dep}')
    cmd+=['--wrap',command]
    job=subprocess.check_output(cmd,text=True).strip().split(';')[0];assert job.isdigit()
    receipt['jobs'].append(dict(name=name,id=job,tasks=tasks,command=cmd));save();print(name,job,flush=True)
    return job
save()
g=submit('greedy',[dict(sample=s,method='Greedy') for s in samples],32,6,4)
py=env+'cassiopeia_env/bin/python'
v=['sbatch','--parsable','--job-name=rb1_fix_verify',f'--dependency=afterany:{g}','--cpus-per-task=1','--mem=4G','--time=00:15:00',
   f'--output={ROOT}/logs/verify_%j.out',f'--error={ROOT}/logs/verify_%j.err','--wrap',f'{py} {ROOT}/scripts/verify_cassiopeia_eight.py {CAS}']
receipt['jobs'].append(dict(name='verify',id=subprocess.check_output(v,text=True).strip().split(';')[0],command=v));save()
full=[s for s in samples if not s.startswith('JZ135_n')]
submit('variants',[dict(sample=s,method=m) for s in full for m in ['NJ','Spectral','MaxCut','Percolation']]
                  +[dict(sample=s,method='NJ') for s in samples if s.startswith('JZ135_n')],16,6,4,dep=g)
submit('lbi',[dict(sample=s,method='LBI') for s in samples],8,2,4,dep=g)
submit('phylovelo',[dict(sample=s,method='PhyloVelo') for s in samples],32,8,2,dep=g)
print(json.dumps(receipt,indent=2))

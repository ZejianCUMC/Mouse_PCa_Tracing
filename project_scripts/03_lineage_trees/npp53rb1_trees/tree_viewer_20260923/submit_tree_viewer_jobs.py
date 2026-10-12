"""Submit the requested bounded analyses once; no automatic retries."""
from pathlib import Path
import os,json,subprocess,shlex,time
ROOT=Path(__file__).resolve().parents[1]
assert str(ROOT).startswith('/groups/')
assert not (ROOT/'SUBMISSION.json').exists()
lock=os.open(ROOT/'SUBMISSION_LOCK',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(lock)
samples=json.loads((ROOT/'PREPARED.json').read_text())['samples']
assert samples==['JZ137','JZ138','JZ256','JZ259','JZ261','JZ262','JZ263','JZ135_n1000','JZ135_n2000','JZ135_n3000']
active=subprocess.check_output(['squeue','-u',os.environ['USER'],'-h','-o','%j'],text=True)
assert 'rb1_view_' not in active
(ROOT/'logs').mkdir(exist_ok=True)
env='/groups/ms3625_gp/zw2994/mamba/envs/'
receipt=dict(submitted_unix=time.time(),jobs=[],held_runs=['MJZ026','MJZ027','MJZ028'],
    JZ135_alternatives='Full-cell run held. Subsamples Greedy+NJ supported; Spectral/MaxCut/Percolation withheld for ambiguous states',
    bounds='CPU only. Variants: 1CPU/16GB/6h max4. LBI:1CPU/8GB/2h max4. PhyloVelo:1CPU/32GB/8h max2. Subsample Greedy:1CPU/16GB/6h max3; subsequent jobs depend on completion.')
def save():(ROOT/'SUBMISSION.json').write_text(json.dumps(receipt,indent=2))
def submit(name,tasks,mem,hours,parallel,dependency=None):
    table=ROOT/(name+'_tasks.json');table.write_text(json.dumps(tasks,indent=2))
    command=f"export PYTHONHASHSEED=42 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1; {env}cassiopeia_env/bin/python {shlex.quote(str(ROOT/'scripts/dispatch_tree_viewer.py'))} {shlex.quote(str(ROOT))} {shlex.quote(str(table))}"
    cmd=['sbatch','--parsable','--job-name=rb1_view_'+name,'--cpus-per-task=1',f'--mem={mem}G',f'--time={hours}:00:00',f'--array=0-{len(tasks)-1}%{parallel}',f'--output={ROOT}/logs/{name}_%A_%a.out',f'--error={ROOT}/logs/{name}_%A_%a.err']
    if dependency:cmd+=['--dependency=afterok:'+dependency]
    cmd+=['--wrap',command]
    job=subprocess.check_output(cmd,text=True).strip().split(';')[0]
    receipt['jobs'].append(dict(name=name,id=job,tasks=tasks,command=cmd));save();print(name,job,flush=True)
    return job
save()
submit('variants',[dict(sample=s,method=m) for s in samples[:7] for m in ['NJ','Spectral','MaxCut','Percolation']],16,6,4)
submit('lbi',[dict(sample=s,method='LBI') for s in samples[:7]],8,2,4)
submit('phylovelo',[dict(sample=s,method='PhyloVelo') for s in samples[:7]],32,8,2)
greedy=submit('JZ135_sub_greedy',[dict(sample=s,method='Greedy') for s in samples[7:]],16,6,3)
submit('JZ135_sub_nj',[dict(sample=s,method='NJ') for s in samples[7:]],16,6,2,greedy)
submit('JZ135_sub_lbi',[dict(sample=s,method='LBI') for s in samples[7:]],8,2,2,greedy)
submit('JZ135_sub_pv',[dict(sample=s,method='PhyloVelo') for s in samples[7:]],32,8,2,greedy)

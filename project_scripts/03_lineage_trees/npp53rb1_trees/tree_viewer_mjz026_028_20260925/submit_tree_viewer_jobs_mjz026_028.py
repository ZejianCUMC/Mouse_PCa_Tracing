"""Submit alternative trees, LBI and PhyloVelo for the 14 MJZ026-028 mice once; no automatic retries.
Adapted from scripts/submit_tree_viewer_jobs.py: same methods/resources, new sample list, job names rb1_view2_*."""
from pathlib import Path
import os,json,subprocess,shlex,time
ROOT=Path(__file__).resolve().parents[1]
assert str(ROOT).startswith('/groups/')
assert not (ROOT/'SUBMISSION.json').exists()
lock=os.open(ROOT/'SUBMISSION_LOCK',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(lock)
samples=json.loads((ROOT/'PREPARED.json').read_text())['samples']
assert len(samples)==14 and all(s.startswith('JZ2') for s in samples)
cas=ROOT.parent/'cassiopeia_mjz026_028_20260925_seed42'
for s in samples: assert (cas/'results'/s/'COMPLETE.json').exists(), f'Greedy tree missing for {s}'
active=subprocess.check_output(['squeue','-u',os.environ['USER'],'-h','-o','%j'],text=True)
assert 'rb1_view2_' not in active
(ROOT/'logs').mkdir(exist_ok=True)
env='/groups/ms3625_gp/zw2994/mamba/envs/'
receipt=dict(submitted_unix=time.time(),jobs=[],held_runs=[],
    bounds='CPU only. Variants: 1CPU/16GB/6h max4. LBI:1CPU/8GB/2h max4. PhyloVelo:1CPU/32GB/8h max2.')
def save():(ROOT/'SUBMISSION.json').write_text(json.dumps(receipt,indent=2))
def submit(name,tasks,mem,hours,parallel):
    table=ROOT/(name+'_tasks.json');table.write_text(json.dumps(tasks,indent=2))
    command=f"export PYTHONHASHSEED=42 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1; {env}cassiopeia_env/bin/python {shlex.quote(str(ROOT/'scripts/dispatch_tree_viewer.py'))} {shlex.quote(str(ROOT))} {shlex.quote(str(table))}"
    cmd=['sbatch','--parsable','--job-name=rb1_view2_'+name,'--cpus-per-task=1',f'--mem={mem}G',f'--time={hours}:00:00',f'--array=0-{len(tasks)-1}%{parallel}',f'--output={ROOT}/logs/{name}_%A_%a.out',f'--error={ROOT}/logs/{name}_%A_%a.err','--wrap',command]
    job=subprocess.check_output(cmd,text=True).strip().split(';')[0]
    receipt['jobs'].append(dict(name=name,id=job,tasks=tasks,command=cmd));save();print(name,job,flush=True)
save()
submit('variants',[dict(sample=s,method=m) for s in samples for m in ['NJ','Spectral','MaxCut','Percolation']],16,6,4)
submit('lbi',[dict(sample=s,method='LBI') for s in samples],8,2,4)
submit('phylovelo',[dict(sample=s,method='PhyloVelo') for s in samples],32,8,2)

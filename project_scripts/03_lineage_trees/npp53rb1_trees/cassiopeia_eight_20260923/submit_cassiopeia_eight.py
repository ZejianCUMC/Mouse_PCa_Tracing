"""Guarded HPC submission with durable IDs; never retry a submitted array."""
import json,subprocess,hashlib,os,datetime
from pathlib import Path
root=Path(__file__).resolve().parents[1]
assert str(root).startswith('/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/')
assert not (root/'SUBMISSION.json').exists()
assert not (root/'results').exists(), 'Existing results require review before submission'
manifest=json.loads((root/'PREPARED.json').read_text())
for record in manifest['records']:
    p=root/'inputs'/record['sample_id']
    for name,key in [('allele_table.csv.gz','allele_table_sha256'),('cell_metadata.csv','metadata_sha256')]:
        assert hashlib.sha256((p/name).read_bytes()).hexdigest()==record[key]
running=subprocess.check_output(['squeue','-u',os.environ['USER'],'-h','-o','%j'],text=True)
assert 'rb1_cas8_greedy' not in running, 'An eight-sample tree array already exists'
fd=os.open(root/'SUBMISSION_LOCK',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
cmd=['sbatch','--parsable','--output='+str(root/'logs/cas_%A_%a.out'),
     '--error='+str(root/'logs/cas_%A_%a.err'),str(root/'scripts/cassiopeia_eight.sbatch')]
job=subprocess.check_output(cmd,text=True,cwd=root).strip().split(';')[0];assert job.isdigit()
record=dict(array_job_id=job,samples=manifest['samples'],array='0-7%4',cpu_per_task=1,memory_per_task='32G',
    time_limit='06:00:00',seed=42,command=cmd,submitted_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(root/'SUBMISSION.json').write_text(json.dumps(record,indent=2))
py='/groups/ms3625_gp/zw2994/mamba/envs/cassiopeia_env/bin/python'
qc=['sbatch','--parsable','--job-name=rb1_cas8_verify','--dependency=afterany:'+job,'--cpus-per-task=1','--mem=4G','--time=00:15:00',
    '--output='+str(root/'logs/verify_%j.out'),'--error='+str(root/'logs/verify_%j.err'),
    '--wrap',f'{py} {root}/scripts/verify_cassiopeia_eight.py {root}']
record['verification_job_id']=subprocess.check_output(qc,text=True,cwd=root).strip().split(';')[0]
(root/'SUBMISSION.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record,indent=2))

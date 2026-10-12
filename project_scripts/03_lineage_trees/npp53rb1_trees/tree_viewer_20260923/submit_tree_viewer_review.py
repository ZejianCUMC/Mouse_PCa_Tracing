import json,subprocess,os,time
from pathlib import Path
root=Path(__file__).resolve().parents[1]
assert str(root).startswith('/groups/') and not (root/'REVIEW_SUBMISSION.json').exists()
lock=os.open(root/'REVIEW_SUBMISSION_LOCK',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(lock)
jobs=json.loads((root/'SUBMISSION.json').read_text())['jobs']
cmd=['sbatch','--parsable','--job-name=rb1_view_review','--cpus-per-task=1','--mem=4G','--time=00:15:00',
     '--dependency=afterany:'+':'.join(j['id'] for j in jobs),f'--output={root}/logs/review_%j.out',f'--error={root}/logs/review_%j.err',
     '--wrap',f'/groups/ms3625_gp/zw2994/mamba/envs/cassiopeia_env/bin/python {root}/scripts/review_tree_viewer_jobs.py']
ident=subprocess.check_output(cmd,text=True).strip().split(';')[0]
(root/'REVIEW_SUBMISSION.json').write_text(json.dumps(dict(id=ident,command=cmd,submitted_unix=time.time(),purpose='Verify and assemble viewer; no retries, no full-JZ135 run, no CARTA'),indent=2))
print('afterany review and viewer build:',ident)

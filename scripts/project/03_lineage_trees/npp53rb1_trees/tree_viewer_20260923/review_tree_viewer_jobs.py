"""Durable afterany review; record failure without retry, build explicit partial viewer."""
import json,subprocess,os,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
assert str(ROOT).startswith('/groups/')
submission=json.loads((ROOT/'SUBMISSION.json').read_text())
ids=','.join(j['id'] for j in submission['jobs'])
text=subprocess.check_output(['sacct','-j',ids,'--format=JobID,State,Elapsed,ExitCode','-P','--noheader'],text=True)
(ROOT/'FINAL_SACCT.txt').write_text(text)
raw={line.split('|')[0]:line.split('|')[1:] for line in text.splitlines() if '|' in line}
status={}
for job in submission['jobs']:
    for i,task in enumerate(job['tasks']):
        row=raw.get(job['id']+'_'+str(i),['UNKNOWN','',''])
        status[task['sample']+'/'+task['method']]=dict(state=row[0].split()[0],elapsed=row[1],exit_code=row[2],array=job['id'],task=i)
(ROOT/'SCHEDULER_STATUS.json').write_text(json.dumps(status,indent=2))
os.environ['RB1_PROJECT']=str(ROOT.parents[1])
os.environ['RB1_VIEWER_TEMPLATE']=str(ROOT/'reference/viewer.html')
os.environ['RB1_VIEWER_PALETTE']=str(ROOT/'reference/palette.json')
from build_tree_viewer import build
build()
review=json.loads((ROOT/'VIEWER_BUILD.json').read_text())
rows=review['status'];expected=['Greedy','NJ','LBI','PhyloVelo']
issues=[]
for row in rows:
    methods=expected if row['analysis'].startswith('JZ135_n') else expected+['Spectral','MaxCut','Percolation']
    for m in methods:
        if not row[m].startswith('Verified'):issues.append(dict(analysis=row['analysis'],method=m,status=row[m]))
(ROOT/'FINAL_REVIEW.json').write_text(json.dumps(dict(status='VERIFIED' if not issues else 'NEEDS_REVIEW',
    reviewed_unix=time.time(),issues=issues,html_sha256=review['html_sha256'],biological_samples=8,datasets=10,
    full_JZ135='HELD',carta='NOT RUN',no_automatic_retry=True),indent=2))
print('FINAL_REVIEW',len(issues),'issues; HTML assembled from verified artifacts only')

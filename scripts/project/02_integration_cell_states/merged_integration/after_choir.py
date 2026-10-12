"""Execute dependent statistics once; stop on any failed prerequisite/boundary."""
import concurrent.futures,json,subprocess,sys,time
from pathlib import Path
project=Path(__file__).resolve().parents[1];out=project/'analysis/merged_20260923_seed42';states=out/'states_seed42';scripts=project/'scripts'
python='/home/jrchristin/miniconda3/envs/ST/bin/python';r='/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript'
lib=project.parent/'Npp53Rb1_Tracing/analysis/tumor_states_20260922/R_library'
status=out/'DEPENDENT_STATS_STATUS.json';assert not status.exists(),'Refuse duplicate dependency worker'
def record(stage,**kw):status.write_text(json.dumps({'stage':stage,'time_unix':time.time(),**kw},indent=2));print(stage,kw,flush=True)
def run(name,command,hours):
 return subprocess.run([python,str(scripts/'run_local.py'),'--log',str(out/'logs'/f'{name}.log'),'--hours',str(hours),*command],check=True)
record('WAITING_FOR_CHOIR_AND_LEIDEN')
deadline=time.monotonic()+6*3600
while not ((states/'CHOIR_COMPLETE.json').exists() and (states/'PREPARED.json').exists()):
 for name in ['choir','prepare_states']:
  receipt=out/'logs'/f'{name}.log.exit.json'
  if receipt.exists() and json.loads(receipt.read_text())['exit_code']!=0:
   record('STOPPED_FAILED_PREREQUISITE',prerequisite=name);sys.exit(1)
 if time.monotonic()>deadline:record('STOPPED_PREREQUISITE_TIMEOUT');sys.exit(1)
 time.sleep(20)
record('RANKING_AMI')
run('rank',[python,str(scripts/'rank.py'),str(states)],.5)
rank=json.loads((states/'RANKED.json').read_text())
fixed_path=out/'USER_FIXED_RESOLUTION_RANGE.json'
fixed=json.loads(fixed_path.read_text()) if fixed_path.exists() else None
if fixed:
 assert fixed['minimum']==.1 and fixed['maximum']==1.2
 assert all(fixed['minimum']<=r<=fixed['maximum'] for r in rank['candidates'])
if rank['boundary_warning'] and not fixed:
 record('NEEDS_RESOLUTION_RANGE_REVIEW',reason='Highest AMI at maximum scanned resolution; do not silently assume optimum')
 sys.exit(0)
record('SCSHC_TOP5_RUNNING',candidates=rank['candidates'],concurrent_jobs=2,user_fixed_range=fixed,boundary_warning=rank['boundary_warning'])
def candidate(i):
 resolution=rank['candidates'][i-1];name=f'scshc_r{resolution:.1f}'
 run(name,[r,str(scripts/'scshc.R'),str(states),str(i),str(lib)],3)
 return name
try:
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  for result in pool.map(candidate,range(1,6)):print('Verified process completion',result,flush=True)
except Exception as e:
 record('STOPPED_STATISTICAL_STAGE_FAILURE',error=str(e));raise
record('SCSHC_TOP5_COMPLETE_REQUIRES_BIOLOGICAL_REVIEW',candidates=rank['candidates'])

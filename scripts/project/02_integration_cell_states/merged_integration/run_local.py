"""Bounded local stage supervisor with durable exit/resource receipts."""
import argparse,json,os,resource,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--log',type=Path,required=True);p.add_argument('--hours',type=float,required=True);p.add_argument('command',nargs=argparse.REMAINDER);args=p.parse_args()
assert 0<args.hours<24 and args.command
receipt=Path(str(args.log)+'.exit.json');assert not receipt.exists() and not args.log.exists(), 'No silent retry/overwrite'
args.log.parent.mkdir(parents=True,exist_ok=True)
start=time.time();env=os.environ.copy()
env.update({'OMP_NUM_THREADS':'8','OPENBLAS_NUM_THREADS':'8','MKL_NUM_THREADS':'8','NUMBA_NUM_THREADS':'8','NUMBA_CACHE_DIR':'/tmp/npp_merged_numba','MPLCONFIGDIR':'/tmp/npp_merged_mpl'})
with args.log.open('w') as log:
 proc=subprocess.Popen(args.command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
 Path(str(args.log)+'.started.json').write_text(json.dumps({'pid':proc.pid,'command':args.command,'start_unix':start,'timeout_hours':args.hours,'supervisor_pid':os.getpid()},indent=2))
 try:code=proc.wait(timeout=args.hours*3600)
 except subprocess.TimeoutExpired:
  import signal
  os.killpg(proc.pid,signal.SIGTERM)
  try:proc.wait(timeout=10)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
  code=124
receipt.write_text(json.dumps({'exit_code':code,'elapsed_seconds':time.time()-start,'command':args.command,'max_child_rss_kb':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss},indent=2))
print(receipt.read_text(),flush=True);raise SystemExit(code)

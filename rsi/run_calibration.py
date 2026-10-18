"""Serial, resumable calibration with bounded child process groups and wall time."""
import runtime as rt
import subprocess,time,sys,os,signal,argparse
from pathlib import Path
children=[];cleaned=set();start=time.monotonic();birth=getattr((rt.ROOT/'PRIOR_ART.md').stat(),'st_birthtime',(rt.ROOT/'PRIOR_ART.md').stat().st_ctime);deadline=start+max(0,birth+5.8*3600-time.time());os.environ['VCC_ASH_WORKERS']='1'
# Controller + scoring process + at most two ASH workers = four scientific processes.
def stop(signum,frame):raise KeyboardInterrupt
signal.signal(signal.SIGTERM,stop)
def kill_group(p):
 if p.pid in cleaned:return
 try:os.killpg(p.pid,signal.SIGTERM)
 except ProcessLookupError:return
 try:p.wait(timeout=30)
 except subprocess.TimeoutExpired:
  try:os.killpg(p.pid,signal.SIGKILL)
  except ProcessLookupError:pass
  p.wait()

def owned(cmd,logpath):
 with open(logpath,'a') as log:
  p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);children.append(p);rt.save('artifacts/process_ownership.json',dict(orchestrator=os.getpid(),children=[c.pid for c in children if c.poll() is None],active_process_group=p.pid,command=cmd))
  with open('artifacts/process_launches.jsonl','a') as f:f.write(rt.json.dumps(dict(pid=p.pid,pgid=p.pid,command=cmd,time=time.time()))+'\n')
  last=time.monotonic()
  while p.poll() is None:
   if time.monotonic()>deadline:kill_group(p);break
   if time.monotonic()-last>1800:rt.status('1 calibration',f'Running {" ".join(cmd[1:])}. Calibration artifacts update after each completed recipe; HOLDOUT sealed.');last=time.monotonic()
   time.sleep(10)
  p.wait();kill_group(p);cleaned.add(p.pid) # Also reap any descendant left by an abnormal scoring-process exit.
  return p.returncode
try:
 parser=argparse.ArgumentParser();parser.add_argument('--tasks',default='H1_rich_search,CD4_matched_search');parser.add_argument('--recipes',default='v1,v2b,v3b,v4a,v4b');parser.add_argument('--skip-preparation',action='store_true');parser.add_argument('--hours',type=float,help='Explicit new invocation wall-clock allocation; otherwise retain original job deadline');args=parser.parse_args()
 if args.hours is not None:
  assert 0<args.hours<=6,'Each invocation is limited to six hours';deadline=start+args.hours*3600
 required=args.tasks.split(',');last=0
 while not all(Path('artifacts/tasks/'+t+'.h5ad').exists() and Path('artifacts/tasks/'+t+'.json').exists() for t in required):
  if time.monotonic()-start>3600:raise RuntimeError('Acquisition did not complete within one hour')
  if time.monotonic()-last>1800:rt.status('1 acquisition','Frozen SEARCH/HOLDOUT panels; waiting for bounded raw extraction. Loop gated.');last=time.monotonic()
  time.sleep(20)
 for script in ([] if args.skip_preparation else ['reproduce_atlas.py','prepare_views.py']):
  assert owned([sys.executable,script],'logs/'+script+'.log')==0,script
 tasks=args.tasks.split(',')
 schedule=[(task,recipe) for task in tasks for recipe in args.recipes.split(',')]
 for task,recipe in schedule:
  if time.monotonic()>deadline:break
  rc=owned([sys.executable,'benchmark.py','--task',task,'--recipes',recipe],f'logs/score_{task}_{recipe}.log')
  if rc:rt.save(f'artifacts/{task}_{recipe}_failure.json',dict(returncode=rc,log=f'logs/score_{task}_{recipe}.log'))

finally:
 for p in children:kill_group(p)
 rt.save('artifacts/process_ownership.json',dict(orchestrator=os.getpid(),children=[],complete=True))

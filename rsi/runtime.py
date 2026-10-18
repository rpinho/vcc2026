"""Resource limits and atomic receipts, all writes confined to this job."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','TOKIO_WORKER_THREADS']:
 os.environ[k]='1'
os.environ['NUMBA_NUM_THREADS']='2'
from pathlib import Path
import json,hashlib,shutil,time,sys,threading
ROOT=Path(__file__).resolve().parent
os.environ['TMPDIR']=str(ROOT/'artifacts/tmp');Path(os.environ['TMPDIR']).mkdir(exist_ok=True,parents=True)
os.environ['NUMBA_CACHE_DIR']=str(ROOT/'artifacts/numba_cache')
sys.path.insert(0,str(ROOT/'vendor'));sys.path.insert(0,str(ROOT/'vendor/deps'))
try:os.nice(19)
except OSError:pass

_guard_lock=threading.Lock()
def guard(extra=0):
 if not os.environ.get('VCC_RSI_STEP'):
  if shutil.disk_usage(ROOT).free-extra<40*1024**3:raise RuntimeError('40 GiB disk reserve')
  return
 # Live workers wait with their durable state intact. The retired one-shot
 # reserve must not abort readers that the supervisor's live gates permit.
 from supervisor import gates
 kind=os.environ.get('VCC_RSI_STEP_KIND','download')
 floor=gates.FREE_DOWNLOAD_FLOOR if kind=='download' else gates.FREE_HARD_FLOOR
 with _guard_lock:
  waiting=False
  while True:
   free=shutil.disk_usage(ROOT).free-extra
   size=gates.job_bytes()+extra
   if free>=floor and size<=gates.JOB_CAP:
    if waiting:print('Live disk guard resumed',flush=True)
    return
   if not waiting:
    print(f'Live disk guard waiting: free after reserve {free/1e9:.2f} GB, job with reserve {size/1e9:.2f} GB',flush=True)
    waiting=True
   time.sleep(30)
def save(path,value):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_name(p.name+'.'+str(os.getpid())+'.tmp');q.write_text(json.dumps(value,indent=2,default=str,allow_nan=False));q.replace(p)
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def status(phase,detail):
 path=ROOT/'STATUS.md';marker='<!-- supervisor:end -->'
 existing=path.read_text() if path.exists() else ''
 prefix=existing.split(marker,1)[0]+marker+'\n\n' if marker in existing else ''
 body=f'# RSI status\n\n{time.strftime("%Y-%m-%d %H:%M UTC",time.gmtime())}\n\nPhase: {phase}\n\n{detail}\n\nDisk free: {shutil.disk_usage(ROOT).free/1024**3:.2f} GiB. Worker PID: {os.getpid()}. No submissions.\n'
 temporary=path.with_name(path.name+'.'+str(os.getpid())+'.tmp')
 temporary.write_text(prefix+body);temporary.replace(path)

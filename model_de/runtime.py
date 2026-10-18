"""Private cache paths, one worker/up to two DE threads, background priority."""
import os,sys
sys.dont_write_bytecode=True
os.environ['PYTHONDONTWRITEBYTECODE']='1'
for key in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS',
            'NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','TOKIO_WORKER_THREADS'):
    os.environ[key] = '1'
_de_threads=int(os.environ.get('VCC_DE_THREADS','1'))
if _de_threads not in (1,2): raise ValueError('VCC_DE_THREADS must be 1 or 2')
os.environ['NUMBA_NUM_THREADS']=str(_de_threads)
import pathlib, json, hashlib, shutil, subprocess, time
ROOT = pathlib.Path(__file__).resolve().parent
VCC = ROOT.parent
for env, sub in [('TMPDIR','tmp'),('NUMBA_CACHE_DIR','numba_cache'),('XDG_CACHE_HOME','cache')]:
    path = ROOT/'artifacts'/sub
    path.mkdir(parents=True, exist_ok=True)
    os.environ[env] = str(path)
try: os.nice(19)
except OSError: pass
if pathlib.Path('/usr/bin/taskpolicy').exists():
    subprocess.run(['/usr/bin/taskpolicy','-b','-p',str(os.getpid())],capture_output=True)
def guard(extra=0, **kwargs):
    while True:
        busy = False
        if shutil.disk_usage(ROOT).free-extra > 24*1024**3:
            # A bounded background worker runs below the kids' scheduler
            # priority. Yield between bounded operations during their active jobs.
            if busy: time.sleep(5)
            return
        print('Waiting for disk reserve',flush=True);time.sleep(20)
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()
def save(path, data):
    path=pathlib.Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2,default=str));tmp.replace(path)
def free_gib(): return shutil.disk_usage(ROOT).free/1024**3

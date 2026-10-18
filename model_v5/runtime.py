"""Job-wide resource and provenance helpers; anonymous reads only."""
import os
for k in ['OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
for k in ['OMP_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','TOKIO_WORKER_THREADS']:os.environ[k]='2'
os.environ['TMPDIR']=os.path.abspath('artifacts/tmp');os.makedirs(os.environ['TMPDIR'],exist_ok=True)
os.environ['NUMBA_CACHE_DIR']=os.path.abspath('artifacts/numba_cache')
os.environ['XDG_CACHE_HOME']=os.path.abspath('artifacts/cache')
os.environ['HF_HOME']=os.path.abspath('artifacts/hf');os.environ['HF_HUB_DISABLE_IMPLICIT_TOKEN']='1'
import pathlib,json,shutil,time,hashlib
ROOT=pathlib.Path(__file__).resolve().parent

def guard(extra=0):
 if shutil.disk_usage('/').free-extra<40*1024**3:raise RuntimeError('40 GiB reserve')
def save(path,data):
 p=pathlib.Path(path);p.parent.mkdir(exist_ok=True,parents=True);q=p.with_suffix('.tmp');q.write_text(json.dumps(data,indent=2,default=str));q.replace(p)
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()

def other_prep_active():
 import subprocess
 listing=subprocess.check_output(['ps','-axo','pid,command'],text=True)
 return any('/.local/bin/vcc prep ' in line and str(ROOT)+'/' not in line for line in listing.splitlines())

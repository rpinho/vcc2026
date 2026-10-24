"""Single numerical worker, low priority, bounded disk; no credentials or remote writes."""
import os,pathlib,json,shutil,time,hashlib
for k in ['OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','TOKIO_WORKER_THREADS']:os.environ[k]='1'
ROOT=pathlib.Path(__file__).resolve().parent
os.environ['TMPDIR']=str(ROOT/'artifacts/tmp');os.environ['NUMBA_CACHE_DIR']=str(ROOT/'artifacts/numba_cache')
(ROOT/'artifacts/tmp').mkdir(parents=True,exist_ok=True)
os.environ['HF_HUB_DISABLE_IMPLICIT_TOKEN']='1'
def guard(extra=0):
 if shutil.disk_usage(ROOT).free-extra<30*1024**3:raise RuntimeError('30 GiB reserve')
def save(path,data):
 p=pathlib.Path(path);p.parent.mkdir(exist_ok=True,parents=True);q=p.with_suffix('.tmp');q.write_text(json.dumps(data,indent=2,default=str));q.replace(p)
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()

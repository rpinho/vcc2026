import os,pathlib,time,json,hashlib,shutil,sys,subprocess
ROOT=pathlib.Path(__file__).resolve().parent
for k in ['OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','TOKIO_WORKER_THREADS']:os.environ[k]='1'
os.environ['NUMBA_CACHE_DIR']=str(ROOT/'artifacts/numba_cache')
os.environ['TMPDIR']=str(ROOT/'artifacts/tmp')
(ROOT/'artifacts/tmp').mkdir(parents=True,exist_ok=True)
(ROOT/'artifacts/tmp').mkdir(parents=True,exist_ok=True)
def save(path,data):
 p=pathlib.Path(path);p.parent.mkdir(exist_ok=True,parents=True);q=p.with_suffix('.tmp');q.write_text(json.dumps(data,indent=2,default=str));q.replace(p)
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def scorer_hashes():
 import cell_eval2
 root=pathlib.Path(cell_eval2.__file__).parent
 return {p:sha(root/p) for p in ['run.py','de.py','de_compute.py','catalog.py','moments.py','prep.py','metrics/de.py','metrics/delta.py','metrics/direction.py','configs/vcc2026.yaml']}
def kid_jobs():
 return False
_guard_after=0.
def guard():
 global _guard_after
 if time.monotonic()<_guard_after:return
 while shutil.disk_usage(ROOT).free<25*1024**3:
  print('PAUSED: disk reserve',flush=True);time.sleep(15)
 _guard_after=time.monotonic()+2.

"""Job-local state; no writes to predecessor models or RSI."""
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS']:
    os.environ[key] = '2'
os.environ['TMPDIR'] = str(ROOT/'artifacts/tmp')
os.environ['NUMBA_CACHE_DIR'] = str(ROOT/'artifacts/numba_cache')
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
(ROOT/'artifacts/tmp').mkdir(parents=True,exist_ok=True)
import json, hashlib, time

def save(path, value):
    p = Path(path); p.parent.mkdir(parents=True,exist_ok=True)
    q = p.with_suffix('.tmp'); q.write_text(json.dumps(value,indent=2,allow_nan=False,default=str)); q.replace(p)

def sha(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(8*1024**2),b''): h.update(chunk)
    return h.hexdigest()

def status(detail):
    print(detail, flush=True)

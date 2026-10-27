"""Read-only inputs, confined outputs and CPU/disk guards for this diagnostic."""
import os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
VCC = ROOT.parents[2]
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','TOKIO_WORKER_THREADS','ARROW_NUM_THREADS']:
    os.environ[key] = '1'
for key in ['TMPDIR','NUMBA_CACHE_DIR','MPLCONFIGDIR','XDG_CACHE_HOME','HF_HOME']:
    p = ROOT/'cache'/key
    p.mkdir(parents=True, exist_ok=True)
    os.environ[key] = str(p)
os.environ['CUDA_VISIBLE_DEVICES'] = ''
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True
import json, shutil, time, hashlib, subprocess
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024**2),b''): h.update(b)
    return h.hexdigest()
def guard(launch=False):
    free=shutil.disk_usage(ROOT).free/1024**3
    if launch:
        print(subprocess.check_output(['df','-g','/'],text=True),flush=True)
    if free < (31.5 if launch else 30.0):
        raise RuntimeError(f'Disk guard: {free:.4f} GiB free')
    return free
def save(name,obj):
    (ROOT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def status(message):
    stamp=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    with open(ROOT/'STATUS.md','a') as f: f.write(f'\n{stamp}: {message}\n')
    print(stamp,message,flush=True)
def check(path,write=False):
    if isinstance(path,(int,type(None))): return
    p=Path(os.fsdecode(path)).resolve()
    if ('h2h2_1026/artifacts/de_cache' in str(p) or p.name=='final_billmA.vcc'):
        raise PermissionError('Forbidden input: '+str(p))
    if write and not p.is_relative_to(ROOT): raise PermissionError('Write outside job: '+str(p))
def audit(event,args):
    if event=='open':
        p,mode,flags=args
        check(p,bool(flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)))
    elif event in ['os.mkdir','os.remove','os.rmdir','os.chmod','os.truncate']: check(args[0],True)
    elif event in ['os.rename','os.link','os.symlink']: check(args[0],True);check(args[1],True)
sys.addaudithook(audit)

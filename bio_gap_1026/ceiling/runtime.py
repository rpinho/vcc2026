"""Job-local version of ../runtime.py: do not import its write-side effects."""
import os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','ARROW_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True
import json, shutil, subprocess, time

def guard():
    result = subprocess.run(['df','-g','/'],capture_output=True,text=True,check=True)
    free = shutil.disk_usage('/').free
    size = sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file())
    if free < 30*1024**3: raise RuntimeError(f'ABORT disk below 30 GiB: {free}')
    if size > 200_000_000: raise RuntimeError('ABORT job size exceeds 200 MB')
    return dict(df_g=result.stdout,free_bytes=free,job_bytes=size)

TRAIN_URL = 'https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/adata_Training.h5ad'
URLS = {TRAIN_URL:15482497461,
 'https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/validation/adata_Validation.h5ad':6928967541,
 'https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/test/adata_Test.h5ad':11950739168}
def enforce(network=False):
    assert os.getpriority(os.PRIO_PROCESS,0) >= 10
    def check(p,write=False):
        if isinstance(p,(int,type(None))): return
        p=Path(os.fsdecode(p)).resolve()
        if p.name=='final_billmA.vcc' or ('h2h2_1026' in p.parts and 'de_cache' in p.parts):
            raise PermissionError('Forbidden input')
        if write and not p.is_relative_to(ROOT): raise PermissionError('Write outside job: '+str(p))
    def audit(event,args):
        if event=='open': check(args[0],bool(args[2] & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)))
        elif event in ('os.mkdir','os.remove','os.rmdir','os.chmod','os.truncate'): check(args[0],True)
        elif event in ('os.rename','os.link','os.symlink'): check(args[0],True); check(args[1],True)
        elif event.startswith('socket.') and not network: raise PermissionError('No network')
    sys.addaudithook(audit)

def save(name,obj):
    path=ROOT/name
    assert path.resolve().is_relative_to(ROOT)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def status(message):
    line=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())+': '+message
    with (ROOT/'STATUS.md').open('a') as f: f.write('\n'+line+'\n')
    print(line,flush=True)

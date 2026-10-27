"""Local resource limits and write/read boundary."""
import os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
REPO = PARENT.parent
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS',
            'RAYON_NUM_THREADS','ARROW_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True
import hashlib, json, shutil, time

def guard(heavy=False):
    free = shutil.disk_usage('/').free / 2**30
    if free < (31.5 if heavy else 30):
        raise RuntimeError(f'Disk reserve: {free:.4f} GiB')
    return free

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8*1024**2), b''):
            h.update(block)
    return h.hexdigest()

def save(path, data):
    p = ROOT/path
    assert p.resolve().is_relative_to(ROOT)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2, allow_nan=False, default=str)+'\n')

def status(message):
    s = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())+': '+message
    with (ROOT/'STATUS.md').open('a') as f:
        f.write('\n'+s+'\n')
    print(s, flush=True)

def enforce():
    assert os.getpriority(os.PRIO_PROCESS, 0) >= 10
    def check(path, write=False):
        if isinstance(path, (int, type(None))):
            return
        p = Path(os.fsdecode(path)).resolve()
        if p.name == 'final_billmA.vcc' or any(x in p.parts for x in ('h2h2_1026','gatecall_1026')):
            raise PermissionError('Forbidden/unneeded held-out input: '+str(p))
        if write and not p.is_relative_to(ROOT):
            raise PermissionError('Write outside align: '+str(p))
    def audit(event, args):
        if event == 'open':
            check(args[0], bool(args[2] & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)))
        elif event in ('os.mkdir','os.remove','os.rmdir','os.chmod','os.truncate'):
            check(args[0], True)
        elif event in ('os.rename','os.link','os.symlink'):
            check(args[0], True); check(args[1], True)
        elif event.startswith('socket.'):
            raise PermissionError('Numerical job has no network')
    sys.addaudithook(audit)

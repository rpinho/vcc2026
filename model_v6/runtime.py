"""model_v6 job-wide helpers: thread caps, scratch dirs inside model_v6, disk floor, kids'-job yield.
Anonymous public reads only. Never submits."""
import os
for k in ['OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[k] = '1'
for k in ['OMP_NUM_THREADS', 'NUMBA_NUM_THREADS', 'POLARS_MAX_THREADS', 'RAYON_NUM_THREADS', 'TOKIO_WORKER_THREADS']:
    os.environ.setdefault(k, '2')
import pathlib, json, shutil, time, hashlib, subprocess, sys
ROOT = pathlib.Path(__file__).resolve().parent
VCC = ROOT.parent
os.environ['TMPDIR'] = str(ROOT / 'artifacts/tmp'); os.makedirs(os.environ['TMPDIR'], exist_ok=True)
os.environ['NUMBA_CACHE_DIR'] = str(ROOT / 'artifacts/numba_cache')
os.environ['XDG_CACHE_HOME'] = str(ROOT / 'artifacts/cache')
os.environ['HF_HOME'] = str(ROOT / 'artifacts/hf'); os.environ['HF_HUB_DISABLE_IMPLICIT_TOKEN'] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.path.insert(0, str(ROOT / 'vendor/deps'))
try:
    os.nice(19)
except OSError:
    pass

FLOOR_GIB = float(os.environ.get('VCC_FLOOR_GIB', 24))  # rsi pauses heavy steps under 20 GB and downloads under 30 GB; stay clear of both where we can
KIDS = ()


def free_gib():
    return shutil.disk_usage('/').free / 1024**3


_KIDS_CACHE = [0.0, False]
def kids_busy():
    return False


def guard(extra=0, floor=FLOOR_GIB, yield_kids=True):
    """Block (never abort) while disk is under the floor or a kids' job is running."""
    t = 0
    while True:
        low = shutil.disk_usage('/').free - extra < floor * 1024**3
        kids = yield_kids and kids_busy()
        if not low and not kids:
            return
        if time.time() - t > 120:
            print('guard waiting:', 'disk %.1f GiB' % free_gib() if low else '', 'kids job' if kids else '', flush=True)
            t = time.time()
        time.sleep(20)


def save(path, data):
    p = pathlib.Path(path); p.parent.mkdir(exist_ok=True, parents=True)
    q = p.with_name(p.name + '.tmp'); q.write_text(json.dumps(data, indent=2, default=str)); q.replace(p)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(b)
    return h.hexdigest()

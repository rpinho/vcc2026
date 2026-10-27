import os, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parent
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','ARROW_NUM_THREADS']:
    os.environ[k]='1'
for k in ['TMPDIR','NUMBA_CACHE_DIR','MPLCONFIGDIR','XDG_CACHE_HOME','HF_HOME']:
    p=ROOT/'scratch'/k;p.mkdir(parents=True,exist_ok=True);os.environ[k]=str(p)
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['PYTHONDONTWRITEBYTECODE']='1';sys.dont_write_bytecode=True
import ast, json, hashlib, shutil, time
def guard(heavy=False):
    free=shutil.disk_usage('/').free/1024**3
    if free < (31.5 if heavy else 30):raise RuntimeError(f'Disk reserve: {free:.3f} GiB')
    return free
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
def save(name,obj):
    p=ROOT/name;assert p.resolve().is_relative_to(ROOT);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,allow_nan=False,default=str)+'\n')
def status(s):
    t=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    with (ROOT/'STATUS.md').open('a') as f:f.write(f'\n{t}: {s}\n')
    print(t,s,flush=True)
def definitions(path,names,ns):
    tree=ast.parse(Path(path).read_text());nodes=[n for n in tree.body if isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name in names]
    assert len(nodes)==len(names);exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),ns)
def enforce():
    assert os.getpriority(os.PRIO_PROCESS,0)>=10
    def check(p,write=False):
        if isinstance(p,(int,type(None))):return
        p=Path(os.fsdecode(p)).resolve()
        if p.name=='final_billmA.vcc' or ('h2h2_1026' in p.parts and 'de_cache' in p.parts):raise PermissionError('Forbidden input')
        if write and not p.is_relative_to(ROOT):raise PermissionError('Write outside job: '+str(p))
    def audit(event,args):
        if event=='open':check(args[0],bool(args[2] & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)))
        elif event in ['os.mkdir','os.remove','os.rmdir','os.chmod','os.truncate']:check(args[0],True)
        elif event in ['os.rename','os.link','os.symlink']:check(args[0],True);check(args[1],True)
        elif event.startswith('socket.'):raise PermissionError('No network in numerical jobs')
    sys.addaudithook(audit)

"""Confined runtime for read-only frozen numerical code."""
import os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parent
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','TOKIO_WORKER_THREADS','ARROW_NUM_THREADS']:
 os.environ[k]='1'
for k in ['TMPDIR','NUMBA_CACHE_DIR','MPLCONFIGDIR','XDG_CACHE_HOME','HF_HOME']:
 p=ROOT/'scratch'/k;p.mkdir(parents=True,exist_ok=True);os.environ[k]=str(p)
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['PYTHONDONTWRITEBYTECODE']='1';sys.dont_write_bytecode=True
import ast,json,hashlib,shutil,time,types,importlib.util
def guard(extra=0):
 if shutil.disk_usage(ROOT).free-extra<30_000_000_000:raise RuntimeError('30 GB decimal disk reserve')
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(path,obj):
 guard();p=(ROOT/path).resolve();assert p.is_relative_to(ROOT);p.parent.mkdir(parents=True,exist_ok=True)
 q=p.with_name(p.name+'.'+str(os.getpid())+'.tmp')
 q.write_text(json.dumps(obj,indent=2,allow_nan=False,default=str)+'\n');q.replace(p)
def utc():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def definitions(filename,names,ns):
 p=ROOT/'inputs'/filename;tree=ast.parse(p.read_text());nodes=[n for n in tree.body if isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name in names]
 assert len(nodes)==len(names);exec(compile(ast.Module(body=nodes,type_ignores=[]),str(p),'exec'),ns);return ns
def frozen():
 import numpy as np,torch
 from torch import nn
 from collections import Counter
 torch.set_num_threads(1)
 ns=dict(np=np,torch=torch,nn=nn,Counter=Counter,rt=types.SimpleNamespace(Path=Path,ROOT=REPO/'model_billmB'))
 definitions('dataset.py',['Sources'],ns);definitions('network.py',['ResponseMLP','blend'],ns);return ns
def modules():
 rt=types.ModuleType('runtime');rt.ROOT=ROOT;rt.VCC=REPO;rt.sha=sha;rt.save=save;rt.guard=guard;rt.json=json
 sys.modules['runtime']=rt
 if 'bounded_score' not in sys.modules:load('bounded_score',ROOT/'inputs/vendor_runtime/bounded_score_one.py')
 return rt
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
_installed=False
def enforce_io():
 global _installed
 if _installed:return
 _installed=True
 def check(p,write=False):
  if isinstance(p,(int,type(None))):return
  p=Path(os.fsdecode(p)).resolve()
  if 'preserve_1025' in p.parts or p.name=='final_billmA.vcc' or '5NeeoWn96lfm23hgExcm' in str(p):raise PermissionError('Forbidden input')
  if write and not p.is_relative_to(ROOT):raise PermissionError('Write outside job: '+str(p))
 def audit(event,args):
  if event=='open':
   p,mode,flags=args;check(p,bool(flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)))
  elif event in ['os.mkdir','os.remove','os.rmdir','os.chmod','os.truncate']:check(args[0],True)
  elif event in ['os.rename','os.link','os.symlink']:check(args[0],True);check(args[1],True)
  elif event.startswith('socket.'):raise PermissionError('Network disabled')
 sys.addaudithook(audit)

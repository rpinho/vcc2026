"""Lossless memoization of the official one-group moment reducer.
No statistic is reimplemented. Exact CSR bytes, shape, genes and normalization
identify the input; the single group's name cannot affect its numeric moments.
"""
import runtime as rt
from pathlib import Path
import hashlib,json,os,collections
import numpy as np
import bounded_score
ORIGINAL=bounded_score.stat;memory=collections.OrderedDict();hits=0;misses=0;namespace=None

def cached(a,cfg):
 global hits,misses
 assert a.obs.target_gene.astype(str).nunique()==1
 x=a.X;h=hashlib.sha256();h.update(namespace);h.update(str((x.shape,cfg.bulk_target_sum)).encode());h.update('\0'.join(a.var_names).encode())
 from scipy import sparse
 if sparse.issparse(x):
  x=x.tocsr();h.update(b'csr');arrays=[x.data,x.indices,x.indptr]
 else:h.update(b'dense');arrays=[np.asarray(x)]
 for v in arrays:h.update(v.dtype.str.encode());h.update(v.tobytes())
 key=h.hexdigest();path=Path('artifacts/stat_cache')/(key+'.npz')
 if key in memory:result=memory.pop(key);hits+=1
 elif path.exists():
  with np.load(path,allow_pickle=False) as d:result=(d['mean'],d['count'][()],d['sumsq'],d['jk'])
  hits+=1
 else:
  misses+=1;result=ORIGINAL(a,cfg);path.parent.mkdir(exist_ok=True);rt.guard(2*1024**2);temp=path.with_name(key+'.'+str(os.getpid())+'.tmp.npz');np.savez_compressed(temp,mean=result[0],count=result[1],sumsq=result[2],jk=result[3]);temp.replace(path)
 memory[key]=result
 while len(memory)>4:memory.popitem(last=False)
 return result

def install_unchecked():
 global namespace
 import official_scoring
 namespace=json.dumps(dict(wrapper_sha256=rt.sha(__file__),official=official_scoring.provenance()),sort_keys=True).encode();bounded_score.stat=cached

def install():
 receipt=rt.json.loads(Path('artifacts/stat_cache_fixture.json').read_text());assert receipt['status']=='PASS bitwise metrics, anchors and relabeled single-group moments' and receipt['implementation_sha256']==rt.sha(__file__)
 import official_scoring
 assert receipt['official_files']==official_scoring.provenance()['official_files'];install_unchecked()

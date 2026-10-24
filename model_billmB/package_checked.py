"""Installed prep only after full learned-profile gates; independent full decoded readback."""
import sys,os,json,hashlib,tarfile,tempfile,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parent.parent/'model_billmB';REPO=R.parent
sys.path.insert(0,str(R));import common as job
import io_guard;io_guard.install()
import numpy as np,pandas as pd,h5py,zstandard as zstd
from scipy import sparse
from fractions import Fraction as F
sys.path.insert(0,str(REPO/'model_v7'));from verify_counts import audit_h5

def group_hashes(path,records):
 assert len(records)==900
 with h5py.File(path) as h:
  X=h['X']
  for k,r in enumerate(records):
   ptr=X['indptr'][k*400:(k+1)*400+1];a,b=map(int,[ptr[0],ptr[-1]])
   x=sparse.csr_matrix((X['data'][a:b],X['indices'][a:b],ptr-a),shape=(400,18533)).toarray()
   assert hashlib.sha256(x.astype('<i4').tobytes()).hexdigest()==r['counts_sha256']
   assert r['context']=='ABC'[k//300] and r['index']==k%300
   if k%100==0:print('GROUP HASH',path.name,k,flush=True)

def main():
 gatepath=REPO/'model_billmB/full_gates.json';g=json.loads(gatepath.read_text());assert g['status']=='PASS all900' and g['full900'] and len(g['records'])==900
 # Independently recheck exact inequalities after rounding interval endpoints
 # outward; the raw integer-count auditor and full decoded audit remain required.
 subprocess.run(['node',str(REPO/'model_billmB/verify_renderer_receipts.mjs'),'--full'],check=True)
 for r in g['records']:
  assert r['all_gates_pass'] and all(r[k] for k in ['totals_exact','depths_exact','nonDE_exact','own_target_exact'])
  assert F(r['jk'])>=F(r['v5b_jk'])
  assert F(r['mean_certificate']['mean_L1_upper'])<F('.015') and F(r['mean_certificate']['bulk_L1_upper'])<F('.0015')
 genes=pd.read_csv(REPO/'data/gene_names.csv').gene_name.to_numpy(str);targets=pd.read_csv(REPO/'data/pert_counts.csv').target_gene.to_numpy(str)
 path=R/'out/billmB.h5ad';group_hashes(path,g['records']);structure=audit_h5(path,genes,targets)
 job.save(R/'artifacts/billmB_structure.json',structure)
 from cache_cleanup import remove
 cache=R/g['cache_directory'];assert cache.resolve().is_relative_to(R/'artifacts/validation_groups')
 caches=list(cache.glob('*.npz'))
 for c in 'ABC':caches+=list((R/f'artifacts/ref_{c}').glob('*.npy'))
 remove(caches,R/'artifacts/validation_preprep_cleanup.json','Retained full900 H5 passed independent structure and all900 dense count hashes; renderer checkpoints no longer required by prep/decode')
 # Keep the validated storage adapter byte-identical, with this job's fail-closed write guard.
 support=REPO/'model_billmB/prep_support';support.mkdir(exist_ok=True)
 adapter=REPO/'model_v7/prep_support/sitecustomize.py'
 code='import sys\nsys.path.insert(0,'+repr(str(REPO/'model_billmB'))+')\nimport io_guard\nio_guard.install()\nexec(compile(open('+repr(str(adapter))+').read(),'+repr(str(adapter))+',"exec"))\n'
 (support/'sitecustomize.py').write_text(code)
 env=os.environ.copy();env.update(PYTHONPATH=str(support),PYTHONDONTWRITEBYTECODE='1',VCC_RESOURCE_LIMITS='1',VCC_MMAP_DATA='1',VCC_DATA_U16='1' if structure['max_count']<=65535 else '0',VCC_DEFER_ENCODING='1',VCC_DELETE_INPUT='0',VCC_KEEP_U16='0',VCC_RESERVE_GIB='25',TMPDIR=str(R/'artifacts/tmp'))
 package=R/'billmB.vcc';assert not list(R.glob('*.vcc')), 'No unverified or duplicate package may be overwritten'
 command=['vcc','prep',str(path),'-g',str(REPO/'data/gene_names.csv'),'--perts',str(REPO/'data/pert_counts.csv'),'--contexts','A,B,C','-o',str(package)]
 p=subprocess.Popen(command,env=env);job.save(R/'artifacts/billmB_prep_running.json',dict(pid=p.pid,command=command));print('PREP PID',p.pid,flush=True)
 # A kids-first SIGSTOP is expected waiting. Do not force-kill a prep child
 # on an interrupted wait; its PID is recorded above for guard-aware recovery.
 code=p.wait()
 job.save(R/'artifacts/billmB_prep.json',dict(exit_code=code,command=command,storage_adapter_sha256=job.sha(adapter)))
 assert code==0
 try:
  with tarfile.open(package) as t:
   assert set(t.getnames())=={'meta.json','pred.h5ad.zst'};meta=json.load(t.extractfile('meta.json'));assert (meta['n_obs'],meta['n_vars'],meta['nnz'])==(360000,18533,structure['nnz'])
   with tempfile.TemporaryDirectory(dir=R/'artifacts/tmp') as tmp:
    decoded=Path(tmp)/'decoded.h5ad'
    with decoded.open('wb') as dst,zstd.ZstdDecompressor().stream_reader(t.extractfile('pred.h5ad.zst')) as reader:
     for block in iter(lambda:reader.read(8*1024**2),b''):dst.write(block)
    exact=audit_h5(decoded,genes,targets);assert exact==structure;group_hashes(decoded,g['records'])
  job.save(R/'artifacts/billmB_package.json',dict(status='PASS installed prep, independent full canonical CSR/axes equality, all900 decoded group hashes and same-profile renderer gates',bytes=package.stat().st_size,sha256=job.sha(package),structure=structure,gates_sha256=job.sha(gatepath),v5b_MSE_preservation='INAPPLICABLE: learned pooled direction changed'))
 except BaseException:
  package.unlink(missing_ok=True);raise
 print('PACKAGE VERIFIED',flush=True)
if __name__=='__main__':main()

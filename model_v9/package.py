"""Local installed vcc prep and independent exact archive decode. Never submit."""
import runtime,sys,os,pathlib,subprocess,signal,time,tarfile,tempfile
import numpy as np,pandas as pd,zstandard as zstd
R=runtime.ROOT
sys.path.append(str(R.parent/'model_v7'))
from verify_counts import audit_h5
def main(name):
 assert name in ['v9a','v9b'];runtime.guard();(R/'artifacts/tmp').mkdir(exist_ok=True)
 genes=pd.read_csv(R.parent/'data/gene_names.csv').gene_name.to_numpy(str);targets=pd.read_csv(R.parent/'data/pert_counts.csv').target_gene.to_numpy(str)
 path=R/f'out/{name}.h5ad';structure=audit_h5(path,genes,targets)
 runtime.save(R/f'artifacts/{name}_structure.json',structure)
 env=os.environ.copy();support=R.parent/'model_v7/prep_support'
 env.update(PYTHONPATH=str(support),PYTHONDONTWRITEBYTECODE='1',VCC_RESOURCE_LIMITS='1',VCC_MMAP_DATA='1',VCC_DATA_U16='1' if structure['max_count']<=65535 else '0',VCC_DEFER_ENCODING='1',VCC_DELETE_INPUT='0',VCC_KEEP_U16='0',VCC_RESERVE_GIB='25',TMPDIR=str(R/'artifacts/tmp'))
 output=R/f'{name}.vcc'
 present=list(R.glob('*.vcc'))
 assert {p.name for p in present}<={'v9a.vcc','v9b.vcc'}
 assert output.exists() or len(present)<2
 cmd=['vcc','prep',str(path),'-g',str(R.parent/'data/gene_names.csv'),'--perts',str(R.parent/'data/pert_counts.csv'),'--contexts','A,B,C','-o',str(output)]
 p=subprocess.Popen(cmd,env=env);paused=False;print('PREP PID',p.pid,flush=True)
 try:
  while p.poll() is None:
   want=bool(runtime.kid_jobs()) or runtime.shutil.disk_usage(R).free<25*1024**3
   if want and not paused:p.send_signal(signal.SIGSTOP);paused=True
   elif paused and not want:p.send_signal(signal.SIGCONT);paused=False
   time.sleep(1)
  code=p.returncode
 finally:
  if p.poll() is None:
   if paused:p.send_signal(signal.SIGCONT)
   p.terminate();p.wait()
 public_command=['vcc','prep',f'model_v9/out/{name}.h5ad','-g','data/gene_names.csv','--perts','data/pert_counts.csv','--contexts','A,B,C','-o',f'model_v9/{name}.vcc']
 runtime.save(R/f'artifacts/{name}_prep.json',dict(cli_exit=code,command=public_command,adapter_sha256=runtime.sha(support/'sitecustomize.py'),structure=structure))
 assert code==0
 with tarfile.open(output,'r') as tar:
  assert set(tar.getnames())=={'meta.json','pred.h5ad.zst'};meta=runtime.json.load(tar.extractfile('meta.json'))
  assert (meta['n_obs'],meta['n_vars'],meta['nnz'])==(360000,18533,structure['nnz'])
  with tempfile.TemporaryDirectory(dir=R/'artifacts/tmp') as tmp:
   decoded=pathlib.Path(tmp)/'decoded.h5ad'
   with decoded.open('wb') as dst,zstd.ZstdDecompressor().stream_reader(tar.extractfile('pred.h5ad.zst')) as reader:
    for block in iter(lambda:reader.read(8*1024**2),b''):runtime.guard();dst.write(block)
   decoded_structure=audit_h5(decoded,genes,targets)
 assert decoded_structure==structure
 receipt=dict(status='PASS installed vcc prep and independent decoded CSR/axes identity',package_bytes=output.stat().st_size,package_sha256=runtime.sha(output),metadata=meta,structure=structure)
 runtime.save(R/f'artifacts/{name}_package.json',receipt);print('PACKAGE VERIFIED',name,receipt['package_sha256'],flush=True)
if __name__=='__main__':main(sys.argv[1])

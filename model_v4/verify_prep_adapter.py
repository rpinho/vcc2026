"""Exercise real installed prep, forced mixed CSR offsets, mapped encoding and exact decoding."""
import runtime,os,subprocess,pathlib,tarfile,tempfile
import numpy as np,pandas as pd,anndata as ad,zstandard as zstd,h5py
from scipy import sparse
from verify_counts import column
rt=runtime.ROOT;folder=rt/'artifacts/prep_fixture';folder.mkdir(exist_ok=True);genes=np.array([f'G{i}' for i in range(64)]);targets=genes[:3];rng=np.random.default_rng(883);x=sparse.csr_matrix(rng.poisson(3,(3600,64)).astype(np.int32));obs=pd.DataFrame(dict(context=np.repeat(list('ABC'),1200),target_gene=np.tile(np.repeat(targets,400),3)),index=[f'fixture_{i}' for i in range(3600)]);a=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=genes));path=folder/'input.h5ad';a.write_h5ad(path,compression='gzip');pd.DataFrame({'gene_name':genes}).to_csv(folder/'genes.csv',index=False);pd.DataFrame({'target_gene':targets,'n_cells':400}).to_csv(folder/'perts.csv',index=False)
env=os.environ.copy();env.update(TMPDIR=str(rt/'artifacts/tmp'),PYTHONPATH=str(rt/'prep_support'),VCC_RESOURCE_LIMITS='1',VCC_MMAP_DATA='1',VCC_DATA_U16='1',VCC_DEFER_ENCODING='1',VCC_MMAP_THRESHOLD='1',VCC_MIXED_INDICES='1',VCC_DELETE_INPUT='1',VCC_RESERVE_GIB='45')
output=folder/'fixture.vcc';cmd=[(__import__('shutil').which('vcc') or 'vcc'),'prep',str(path),'-g',str(folder/'genes.csv'),'--perts',str(folder/'perts.csv'),'--contexts','A,B,C','--expected-gene-dim','64','-o',str(output),'--force'];r=subprocess.run(cmd,env=env,capture_output=True,text=True);print(r.stdout,flush=True);print(r.stderr,flush=True);assert r.returncode==0;assert 'defer exact float32 encoding to bounded HDF5 writer' in r.stderr+r.stdout;assert not path.exists()
with tarfile.open(output,'r') as tar:
 assert set(tar.getnames())=={'meta.json','pred.h5ad.zst'}
 with tempfile.TemporaryDirectory(dir=rt/'artifacts/tmp') as temp:
  decoded=pathlib.Path(temp)/'decoded.h5ad'
  with decoded.open('wb') as f,zstd.ZstdDecompressor().stream_reader(tar.extractfile('pred.h5ad.zst')) as reader:
   for block in iter(lambda:reader.read(2**20),b''):f.write(block)
  with h5py.File(decoded,'r') as h:
   assert h['X/data'].dtype==np.float32;assert h['X/indptr'].dtype==np.int64;assert h['X/indices'].dtype==np.int32
   assert np.array_equal(h['X/data'][:],x.data);assert np.array_equal(h['X/indices'][:],x.indices);assert np.array_equal(h['X/indptr'][:],x.indptr)
   assert np.array_equal(column(h['obs/context']),obs.context.to_numpy());assert np.array_equal(column(h['obs/target_gene']),obs.target_gene.to_numpy());assert np.array_equal(column(h['var'][h['var'].attrs['_index']]),genes)
runtime.save(rt/'artifacts/prep_adapter_verification.json',dict(status='VERIFIED real installed prep and exact decoded CSR/axis equality under forced int64-offset/int32-column mode',shape=[3600,64],nnz=x.nnz,input_deleted_after_mapping=True,encoding='uint16 disk staging with exact bounded row sums and deferred float32 HDF5 encoding; no validator changed'))
output.unlink()
# Ensure the bounded sum preserves the installed validator's count-cap refusal.
bad=a.copy();bad.X[0,:]=20000;bad.write_h5ad(path,compression='gzip');env['VCC_DELETE_INPUT']='0'
r=subprocess.run(cmd,env=env,capture_output=True,text=True);assert r.returncode!=0 and 'per-cell count total above' in r.stderr+r.stdout,(r.returncode,r.stderr);assert not output.exists()
print('Prep adapter VERIFIED, including installed over-depth refusal',flush=True)

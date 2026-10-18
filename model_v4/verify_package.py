"""Independently decode the final package and match exact canonical arrays."""
import runtime,tarfile,tempfile,pathlib,sys
import pandas as pd,zstandard as zstd
from verify_counts import audit_h5
rt=runtime.ROOT;arm=sys.argv[1];assert arm in ['a','b'];name='prediction_v4'+arm;path=rt/'out'/f'{name}.vcc';before=runtime.json.loads((rt/'artifacts'/f'{name}_structure.json').read_text());genes=pd.read_csv(rt.parent/'data/gene_names.csv').gene_name.to_numpy(str);targets=pd.read_csv(rt.parent/'data/pert_counts.csv').target_gene.to_numpy(str)
with tarfile.open(path,'r') as tar:
 assert set(tar.getnames())=={'meta.json','pred.h5ad.zst'};meta=runtime.json.load(tar.extractfile('meta.json'));assert meta['n_obs']==360000 and meta['n_vars']==18533 and meta['nnz']==before['nnz'] and meta['schema']==1
 with tempfile.TemporaryDirectory(dir=rt/'artifacts/tmp') as tmp:
  decoded=pathlib.Path(tmp)/'decoded.h5ad'
  with decoded.open('wb') as f,zstd.ZstdDecompressor().stream_reader(tar.extractfile('pred.h5ad.zst')) as reader:
   for block in iter(lambda:reader.read(8*1024**2),b''):
    runtime.guard(len(block));f.write(block)
  after=audit_h5(decoded,genes,targets);runtime.save(rt/'artifacts'/f'{name}_decoded_structure.json',after);assert before==after,{k:(before[k],after[k]) for k in before if before[k]!=after[k]}
manifest=runtime.json.loads((rt/'out'/f'{name}.manifest.json').read_text())
for c in manifest['controls']:assert runtime.sha(rt.parent/'data'/c['file'])==c['sha256'],'Input controls changed'
cert=dict(status='VERIFIED installed prep success plus exact decoded CSR/axes equality',input=before,package_metadata=meta,package_bytes=path.stat().st_size,package_sha256=runtime.sha(path))
runtime.save(rt/'artifacts'/f'{name}_package.json',cert);print(runtime.json.dumps(cert,indent=2),flush=True)

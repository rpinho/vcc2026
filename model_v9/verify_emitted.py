"""Bind every independent DE/JK gate record to the actually written cell block."""
import runtime,sys,hashlib
import numpy as np,h5py
from scipy import sparse
R=runtime.ROOT
name=sys.argv[1];assert name in ['v9a','v9b']
records=runtime.json.loads((R/f'artifacts/{name}_gates.json').read_text())['records'];assert len(records)==900
with h5py.File(R/f'out/{name}.h5ad') as f:
 X=f['X'];ptr=X['indptr'][:]
 for k,rec in enumerate(records):
  runtime.guard();lo=k*400;a,b=int(ptr[lo]),int(ptr[lo+400]);x=sparse.csr_matrix((X['data'][a:b],X['indices'][a:b],ptr[lo:lo+401]-a),shape=(400,18533)).toarray()
  actual=hashlib.sha256(x.astype('<i4').tobytes(order='C')).hexdigest();assert actual==rec['counts_sha256'],(name,k,'gate counts differ from written counts')
  assert rec['context']=='ABC'[k//300] and rec['index']==k%300
  if k%100==0:print('EXACT WRITTEN GROUP HASH',name,k,flush=True)
runtime.save(R/f'artifacts/{name}_emitted_identity.json',dict(status='PASS all 900 written cell blocks exactly match independent DE/JK/moment gate input hashes',groups=900,hash='SHA-256 canonical little-endian int32 dense rows',gates_sha256=runtime.sha(R/f'artifacts/{name}_gates.json')))

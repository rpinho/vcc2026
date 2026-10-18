"""Compare all index addresses and selected data with ordinary HDF5, no network."""
import runtime as rt
import os,io,collections,threading,h5py,numpy as np
from pathlib import Path
from tree_prefetch import PrefetchRemote,read_uncompressed_vector
path=Path('artifacts/tmp/tree_fixture.h5')
with h5py.File(path,'w',libver='earliest') as h:h.create_dataset('x',data=np.arange(1_000_000,dtype=np.int32),chunks=(31,))
def index(h):
 out=[];h['x'].id.chunk_iter(lambda z:out.append((tuple(z.chunk_offset),z.byte_offset,z.size,z.filter_mask)));return out
with h5py.File(path,'r') as h:expected=index(h)
class Local(PrefetchRemote):
 def __init__(self,path):
  io.RawIOBase.__init__(self);self.fd=os.open(path,os.O_RDONLY);self._byte_lock=threading.Lock();self._local=threading.local();self._pool=None;self.size=path.stat().st_size;self.pos=0;self.bytes=0;self.cap=20_000_000;self.block=65536;self.cache=collections.OrderedDict()
 def _read_range(self,start,n):return os.pread(self.fd,n,start)
r=Local(path)
with h5py.File(r,'r') as h:np.testing.assert_array_equal(read_uncompressed_vector(h['x'],r),np.arange(1_000_000,dtype=np.int32))
r.enable_tree_prefetch()
with h5py.File(r,'r') as h:
 actual=index(h);assert actual==expected;ix=np.array([0,31,9999,999999]);np.testing.assert_array_equal(h['x'][ix],ix.astype(np.int32))
internal=0;children=0
for address in r._seen:
 node=os.pread(r.fd,2096,address)
 if node[4]==1 and node[5]>0:
  internal+=1;n=int.from_bytes(node[6:8],'little')
  for i in range(n):
   child=int.from_bytes(node[48+32*i:56+32*i],'little');assert os.pread(r.fd,4,child)==b'TREE';children+=1
assert internal>0 and children>0,'Read-ahead hints were not exercised'
r.disable_tree_prefetch();os.close(r.fd);path.unlink();rt.save('artifacts/tree_prefetch_fixture.json',dict(status='PASS bitwise index addresses and data',chunks=len(actual),bytes=r.bytes,internal_nodes_checked=internal,child_hints_verified=children,implementation_sha256=rt.sha('tree_prefetch.py'),scope='Read-ahead only; HDF5 retains all interpretation and chunk-address decisions'));print('PASS HDF5 tree read-ahead',len(actual),'chunks',r.bytes,'bytes',internal,'internal nodes',children,'verified children')

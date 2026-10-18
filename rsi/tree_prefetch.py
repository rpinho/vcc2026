"""Read-ahead for HDF5 v1 chunk-tree metadata. HDF5 still interprets every node.
Only speculative HTTP GETs are scheduled; no chunk address returned to the caller
is computed here. Layout: HDF Group's File Format Specification, v1 B-link trees.
"""
import runtime as rt
from remote_h5 import Remote
import threading,requests,time
import numpy as np
from concurrent.futures import ThreadPoolExecutor

def read_uncompressed_vector(dataset,remote):
 """Official HDF5 chunk locations, exact payload ranges, bounded parallel reads."""
 assert dataset.ndim==1 and dataset.chunks and dataset.compression is None and not dataset.shuffle
 info=[];dataset.id.chunk_iter(lambda z:info.append((int(z.chunk_offset[0]),z.byte_offset)));out=np.empty(dataset.shape,dataset.dtype);chunk=dataset.chunks[0]
 def fetch(item):
  start,address=item;n=min(chunk,len(out)-start);return start,n,remote.fetch(address,n*dataset.dtype.itemsize)
 with ThreadPoolExecutor(max_workers=32) as pool:
  for base in range(0,len(info),64):
   for start,n,value in pool.map(fetch,info[base:base+64]):out[start:start+n]=np.frombuffer(value,dtype=dataset.dtype,count=n)
 return out

class PrefetchRemote(Remote):
 def __init__(self,*args,**kwargs):
  self._byte_lock=threading.Lock();self._local=threading.local();self._pool=None;super().__init__(*args,**kwargs)
 def _read_range(self,start,n):
  if not hasattr(self._local,'http'):self._local.http=requests.Session();self._local.http.trust_env=False
  headers={'Range':f'bytes={start}-{start+n-1}','Accept-Encoding':'identity'}
  if self.etag:headers['If-Match']=self.etag
  with self._local.http.get(self.url,headers=headers,stream=True,timeout=(10,60)) as response:
   assert response.status_code==206;value=response.raw.read(n+1);assert len(value)==n;return value
 def fetch(self,start,n):
  for attempt in range(3):
   rt.guard()
   with self._byte_lock:assert self.bytes+n<=self.cap,'Metadata transfer cap';self.bytes+=n
   try:return self._read_range(start,n)
   except Exception:
    if attempt==2:raise
    time.sleep(1)
 def enable_tree_prefetch(self):
  assert self._pool is None;self.block=4096;self.cache.clear();self._pending={};self._nodes={};self._active_node=None;self._seen=set();self._pool=ThreadPoolExecutor(max_workers=32)
 def _schedule(self,block):
  if block not in self._pending:
   start=block*self.block
   if start<self.size:self._pending[block]=self._pool.submit(self.fetch,start,min(self.block,self.size-start))
 def _get(self,block):self._schedule(block);return self._pending[block].result()
 def _peek(self,start,n):
  if start in self._nodes and n<=2096:
   value=self._nodes[start].result();self._active_node=(start,value);return value[:n]
  if self._active_node is not None:
   address,value=self._active_node
   if address<=start and start+n<=address+len(value):return value[start-address:start-address+n]
  end=min(self.size,start+n);out=[]
  while start<end:
   b=start//self.block;offset=start%self.block;v=self._get(b);take=min(end-start,len(v)-offset);out.append(v[offset:offset+take]);start+=take
  return b''.join(out)
 def read(self,n=-1):
  if self._pool is None:return super().read(n)
  if n<0:n=self.size-self.pos
  start=self.pos;value=self._peek(start,n);self.pos+=len(value)
  if value[:4]==b'TREE' and start not in self._seen:
   self._seen.add(start);header=self._peek(start,8);kind,level=header[4],header[5];count=int.from_bytes(header[6:8],'little')
   if kind==1 and 0<level<8 and 0<count<=64:
    # Source arrays are rank one, with the datatype-element placeholder dimension;
    # 8-byte offsets, 24-byte keys. A wrong hint can only waste a capped read.
    node=self._peek(start,48+32*count)
    for i in range(count):
     child=int.from_bytes(node[48+32*i:56+32*i],'little')
     if 0<child<self.size-2096 and child not in self._nodes:self._nodes[child]=self._pool.submit(self.fetch,child,2096)
  return value
 def disable_tree_prefetch(self):
  if self._pool is not None:self._pool.shutdown(wait=True);self._pool=None;self._pending.clear();self._nodes.clear();self._active_node=None;self.cache.clear()

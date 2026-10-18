"""Resource-only VCC CLI adapter; all installed validation rules run unchanged.
Keep count values in RAM and column indices on disk. After the installed writer
serializes the CSR, release the index mapping before compression/tar duplication.
The installed writer subsequently reads only nnz/shape (indptr retained in RAM).
"""
import os
if os.environ.get('VCC_RESOURCE_LIMITS')=='1':
 import atexit,pathlib,tempfile,shutil,sys
 import anndata as ad,numpy as np,h5py
 from scipy.sparse import csr_matrix
 import vcc.prep as prep
 prep._ZSTD_MAX_WORKERS=3
 _maps=[]
 _mapped_matrix=None
 original_astype=csr_matrix.astype
 def bounded_encoding(self,dtype,casting='unsafe',copy=True):
  if os.environ.get('VCC_DEFER_ENCODING')=='1' and self is _mapped_matrix and self.dtype==np.uint16 and np.dtype(dtype)==np.float32 and not copy:
   print('VCC resource adapter: defer exact float32 encoding to bounded HDF5 writer',file=sys.stderr,flush=True)
   return self
  # Preserve int32 columns plus int64 row pointers beyond 2^31 NNZ.
  # Encoding changes only exactly representable count values, in the existing mapped buffer.
  if self is _mapped_matrix and np.dtype(dtype)==np.dtype('float32') and self.data.dtype==np.dtype('int32') and not copy and casting=='unsafe':
   counts=self.data
   if counts.size and (counts.min()<0 or counts.max()>1000000):return original_astype(self,dtype,casting=casting,copy=copy)
   encoded=counts.view(np.float32)
   for lo in range(0,len(counts),2**20):encoded[lo:lo+2**20]=counts[lo:lo+2**20].astype(np.float32)
   self.data=encoded
   print('VCC resource adapter: exact bounded float32 encoding without index duplication',file=sys.stderr,flush=True)
   return self
  return original_astype(self,dtype,casting=casting,copy=copy)
 csr_matrix.astype=bounded_encoding
 original_sum=csr_matrix.sum
 def bounded_row_sum(self,axis=None,dtype=None,out=None):
  # Integral sums are exact, with SciPy's same default signed/unsigned accumulator.
  # Avoid SciPy promoting the ENTIRE mapped data array to int64 for one row reduction.
  if self is _mapped_matrix and axis in (1,-1) and dtype is None and self.data.dtype.kind in 'iu':
   dt=np.uint64 if self.data.dtype.kind=='u' else np.int64
   totals=np.empty(self.shape[0],dtype=dt);ptr=self.indptr
   for lo in range(0,self.shape[0],128):
    hi=min(lo+128,self.shape[0]);a,b=int(ptr[lo]),int(ptr[hi]);cs=np.empty(b-a+1,dtype=dt);cs[0]=0;np.cumsum(self.data[a:b],dtype=dt,out=cs[1:]);totals[lo:hi]=cs[ptr[lo+1:hi+1]-a]-cs[ptr[lo:hi]-a]
   result=np.asmatrix(totals[:,None])
   if out is not None:out[...] = result;return out
   print('VCC resource adapter: exact bounded integral row sums',file=sys.stderr,flush=True)
   return result
  return original_sum(self,axis=axis,dtype=dtype,out=out)
 csr_matrix.sum=bounded_row_sum
 def release_indices():
  while _maps:
   target,directory=_maps.pop();target.flush();target._mmap.close()
   try:directory.cleanup()
   except Exception:pass
 atexit.register(release_indices)
 original_write=ad.AnnData.write_h5ad
 def compressed_write(self,*args,**kwargs):
  kwargs.setdefault('compression','gzip');kwargs.setdefault('compression_opts',4)
  if args and pathlib.Path(args[0]).name=='pred.h5ad':print(f'VCC serialized index dtypes: {self.X.indices.dtype}, {self.X.indptr.dtype}',file=sys.stderr,flush=True)
  result=original_write(self,*args,**kwargs)
  # Installed _write_vcc's only full serialization, after every validation passed.
  if args and pathlib.Path(args[0]).name=='pred.h5ad':
   release_indices();print('VCC resource adapter: serialized CSR; released index maps before compression',file=sys.stderr,flush=True)
   if os.environ.get('VCC_KEEP_U16')!='1':_restore_float32(pathlib.Path(args[0]))
   else:print('VCC resource adapter: keeping uint16 X/data in the serialized pred.h5ad (values are the same integers)',file=sys.stderr,flush=True)
  return result
 def _restore_float32(path):
  with h5py.File(path,'r') as h:
   if h['X/data'].dtype==np.float32:return
  tmp=path.with_name('pred.f32.h5ad')
  with h5py.File(path,'r') as src,h5py.File(tmp,'w') as dst:
   for k,v in src.attrs.items():dst.attrs[k]=v
   def copy(name,obj):
    if name=='X/data':
     d=dst.create_dataset(name,shape=obj.shape,dtype=np.float32,chunks=obj.chunks,compression=obj.compression,compression_opts=obj.compression_opts,shuffle=obj.shuffle)
     for lo in range(0,obj.shape[0],2**22):d[lo:lo+2**22]=obj[lo:lo+2**22].astype(np.float32)
     for k,v in obj.attrs.items():d.attrs[k]=v
    elif isinstance(obj,h5py.Group):
     g=dst.require_group(name)
     for k,v in obj.attrs.items():g.attrs[k]=v
    else:
     src.copy(obj,dst,name=name)
   src.visititems(copy)
  with h5py.File(tmp,'r') as h:assert h['X/data'].dtype==np.float32
  os.replace(tmp,path);print('VCC resource adapter: X/data restored to float32 in serialized pred.h5ad (values unchanged)',file=sys.stderr,flush=True)
 ad.AnnData.write_h5ad=compressed_write
 original_create=h5py.Group.create_dataset
 def fixed_csr_chunks(self,name,*args,**kwargs):
  if self.name.endswith('/X') and name in ('data','indices','indptr'):
   value=kwargs.get('data');length=value.shape[0] if value is not None else int(kwargs.get('shape',(0,))[0])
   if length:
    kwargs['chunks']=(min(length,4096 if name=='indptr' else 262144),);kwargs['shuffle']=True
  return original_create(self,name,*args,**kwargs)
 h5py.Group.create_dataset=fixed_csr_chunks
 original_read=prep.read_h5ad
 def mapped_read(path):
  global _mapped_matrix
  with h5py.File(path,'r') as h:
   if 'X' not in h or not isinstance(h['X'],h5py.Group) or h['X'].attrs.get('encoding-type')!='csr_matrix':return original_read(path)
   g=h['X'];nnz=len(g['data'])
   if nnz<int(os.environ.get('VCC_MMAP_THRESHOLD','100000000')):return original_read(path)
   shape=tuple(map(int,g.attrs['shape']));index_dtype=np.int32
   mixed=nnz>np.iinfo(np.int32).max or os.environ.get('VCC_MIXED_INDICES')=='1'
   pointer_dtype=np.int64 if mixed else np.int32
   # At most one additional HDF5 plus mapped indices before serialization completes;
   # the index map is then closed before the compressed copy and tar are created.
   maps_bytes=nnz*np.dtype(index_dtype).itemsize+(nnz*(2 if os.environ.get('VCC_DATA_U16')=='1' else g['data'].dtype.itemsize) if os.environ.get('VCC_MMAP_DATA')=='1' else 0)
   # index/data maps are released before the compressed copy is made, so the peak is the larger of the two stages
   needed=max(maps_bytes,1.3*os.path.getsize(path))+512*1024**2
   if shutil.disk_usage('/').free-needed<float(os.environ.get('VCC_RESERVE_GIB','60'))*1024**3:raise prep.PrepError(f"Hybrid mapped prep would cross the {os.environ.get('VCC_RESERVE_GIB','60')} GiB disk reserve")
   directory=tempfile.TemporaryDirectory(prefix='vcc-indices-',dir=os.environ['TMPDIR']);arrays={}
   for name,dtype in [('data',g['data'].dtype),('indices',index_dtype),('indptr',pointer_dtype)]:
    source=g[name]
    if name=='data' and os.environ.get('VCC_DATA_U16')=='1':dtype=np.uint16
    if name=='indices' or (name=='data' and os.environ.get('VCC_MMAP_DATA')=='1'):
     target=np.memmap(pathlib.Path(directory.name)/f'{name}.bin',mode='w+',dtype=dtype,shape=source.shape);_maps.append((target,directory))
    else:target=np.empty(source.shape,dtype=dtype)
    for lo in range(0,len(source),2**20):
     if shutil.disk_usage('/').free<float(os.environ.get('VCC_RESERVE_GIB','60'))*1024**3:raise prep.PrepError(f"{os.environ.get('VCC_RESERVE_GIB','60')} GiB disk reserve reached")
     block=source[lo:lo+2**20]
     if name=='data' and dtype==np.uint16:
      if block.size and (block.max()>65535 or block.min()<0 or np.any(block!=np.floor(block))):raise prep.PrepError('uint16 data map requires integer counts in [0,65535]')
     target[lo:lo+2**20]=block
    if isinstance(target,np.memmap):target.flush()
    arrays[name]=target
   if mixed:
    ptr=arrays['indptr'];ix=arrays['indices']
    if ptr.shape!=(shape[0]+1,) or ptr[0]!=0 or ptr[-1]!=nnz or np.any(np.diff(ptr)<0):raise prep.PrepError('Invalid CSR row pointers')
    # HDF5 CSR supports int32 columns and int64 offsets independently. Verify canonical
    # rows before marking their true status, avoiding SciPy's unnecessary int64 column
    # expansion. Every installed structural/value/count validator still runs afterwards.
    for row in range(shape[0]):
     segment=ix[ptr[row]:ptr[row+1]]
     if np.any(np.diff(segment)<=0):raise prep.PrepError('Mixed-index resource mode requires sorted unique CSR columns')
    matrix=csr_matrix(shape,dtype=arrays['data'].dtype);matrix.data=arrays['data'];matrix.indices=ix;matrix.indptr=ptr
    matrix.has_sorted_indices=True;matrix.has_canonical_format=True
   else:matrix=csr_matrix((arrays['data'],arrays['indices'],arrays['indptr']),shape=shape,copy=False)
   _mapped_matrix=matrix
   result=ad.AnnData(matrix,obs=ad.io.read_elem(h['obs']),var=ad.io.read_elem(h['var']))
   for name in ('uns','obsm','varm','obsp','varp','layers'):
    if name in h:
     value=ad.io.read_elem(h[name])
     if name=='layers':
      if value:result.layers.update(value)
     else:setattr(result,name,value)
   print(f'VCC resource adapter: {nnz:,} entries; data/indices mapped (data dtype {arrays["data"].dtype}); installed validation unchanged',file=sys.stderr,flush=True)
   if os.environ.get('VCC_DELETE_INPUT')=='1':
    os.remove(path);print(f'VCC resource adapter: removed regenerable input {path} after mapping (seeded rebuild)',file=sys.stderr,flush=True)
   return result
 prep.read_h5ad=mapped_read

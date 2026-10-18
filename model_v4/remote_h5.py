"""Read-only, bounded HTTP range file with explicit anonymous sessions."""
import runtime
import io,collections,requests,time
class Remote(io.RawIOBase):
 def __init__(self,url,cap=3_000_000_000,block=262144):
  self.url=url;self.http=requests.Session();self.http.trust_env=False
  h=self.http.head(url,timeout=40);h.raise_for_status();self.size=int(h.headers.get('Content-Length',0));self.etag=h.headers.get('ETag')
  if not self.size:
   with self.http.get(url,headers={'Range':'bytes=0-0','Accept-Encoding':'identity'},stream=True,timeout=40) as probe:
    assert probe.status_code==206;self.size=int(probe.headers['Content-Range'].split('/')[-1]);self.etag=probe.headers.get('ETag');assert len(probe.raw.read(2))==1
  self.pos=0;self.bytes=0;self.cap=cap;self.block=block;self.cache=collections.OrderedDict()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):
  self.pos=n if whence==0 else self.pos+n if whence==1 else self.size+n
  return self.pos
 def fetch(self,start,n):
  runtime.guard();assert self.bytes+n<=self.cap,'HTTP transfer cap'
  headers={'Range':f'bytes={start}-{start+n-1}','Accept-Encoding':'identity'}
  if self.etag:headers['If-Match']=self.etag
  for attempt in range(5):
   try:
    with self.http.get(self.url,headers=headers,timeout=(10,60),stream=True) as r:
     assert r.status_code==206,(r.status_code,self.url);assert int(r.headers.get('Content-Length',n))==n
     value=r.raw.read(n+1);assert len(value)==n
    self.bytes+=n;return value
   except Exception:
    self.bytes+=n
    if attempt==4:raise
    time.sleep(2)
 def readinto(self,b):
  v=self.read(len(b));b[:len(v)]=v;return len(v)
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  end=min(self.size,self.pos+n);parts=[]
  while self.pos<end:
   block=self.pos//self.block;start=block*self.block
   if block not in self.cache:
    self.cache[block]=self.fetch(start,min(self.block,self.size-start))
    if len(self.cache)>512:self.cache.popitem(last=False)
   value=self.cache[block];offset=self.pos-start;take=min(end-self.pos,len(value)-offset)
   parts.append(value[offset:offset+take]);self.pos+=take
  return b''.join(parts)

"""Bounded in-memory range adapter, adapted from model_v6/remote_h5.py."""
import runtime as rt
import io, collections, time, requests

class Remote(io.RawIOBase):
    def __init__(self,url=rt.TRAIN_URL):
        assert url in rt.URLS
        self.http=requests.Session(); self.http.trust_env=False
        self.url=url; self.pos=0; self.bytes=0; self.requests=0
        self.block=4*2**20; self.cache=collections.OrderedDict()
        r=self.http.head(self.url,timeout=40,allow_redirects=False)
        assert r.status_code==200
        self.size=int(r.headers['Content-Length']); self.etag=r.headers['ETag']
        assert self.size==rt.URLS[url]
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self,n,whence=0):
        self.pos=n if whence==0 else self.pos+n if whence==1 else self.size+n
        assert 0<=self.pos<=self.size
        return self.pos
    def fetch(self,start,n):
        assert self.bytes+n<35_000_000_000
        for attempt in range(10):
            try:
                with self.http.get(self.url,headers={'Range':f'bytes={start}-{start+n-1}',
                    'Accept-Encoding':'identity','If-Match':self.etag},stream=True,
                    timeout=(15,90),allow_redirects=False) as r:
                    assert r.status_code==206
                    assert r.headers['Content-Range']==f'bytes {start}-{start+n-1}/{self.size}'
                    value=r.raw.read(n+1); assert len(value)==n
                self.bytes+=n; self.requests+=1
                return value
            except Exception as e:
                print('range retry',attempt+1,start,n,str(e)[:180],flush=True)
                if attempt==9: raise
                time.sleep(min(2**attempt,30))
    def readinto(self,b):
        v=self.read(len(b)); b[:len(v)]=v; return len(v)
    def read(self,n=-1):
        if n<0:n=self.size-self.pos
        end=min(self.size,self.pos+n); parts=[]
        while self.pos<end:
            block=self.pos//self.block; start=block*self.block
            if block not in self.cache:
                self.cache[block]=self.fetch(start,min(self.block,self.size-start))
                if len(self.cache)>16:self.cache.popitem(last=False)
            value=self.cache[block]; offset=self.pos-start
            take=min(end-self.pos,len(value)-offset)
            parts.append(value[offset:offset+take]); self.pos+=take
        return b''.join(parts)
    def close(self):
        self.http.close(); super().close()

"""Ephemeral anonymous HTTP range bridge for public Hugging Face Lance files.
Caches public redirect URLs in memory, never reads account tokens or config.
"""
import runtime
from http.server import ThreadingHTTPServer,BaseHTTPRequestHandler
import threading,requests,time,urllib.parse
class PublicStore:
 def __init__(self):
  self.http=requests.Session();self.http.trust_env=False;self.http.mount('https://',requests.adapters.HTTPAdapter(pool_connections=64,pool_maxsize=64));self.links={};self.lock=threading.Lock();self.bytes=0;self.requests=0;self.errors=0
  parent=self
  class Handler(BaseHTTPRequestHandler):
   protocol_version='HTTP/1.1'
   def log_message(self,*args):pass
   def do_HEAD(self):self.serve(True)
   def do_GET(self):self.serve(False)
   def serve(self,head):
    try:
     path=urllib.parse.urlsplit(self.path).path
     assert path.startswith('/public/');key=path[len('/public/'):]
     assert '..' not in key
     with parent.lock:
      if key not in parent.links:
       s=parent.http
       url='https://huggingface.co/datasets/slaf-project/X-Atlas-Orion/resolve/main/'+key
       for attempt in range(6):
        rr=s.head(url,allow_redirects=True,timeout=60)
        if rr.status_code==429:time.sleep(40);continue
        rr.raise_for_status();parent.links[key]=(rr.url,requests.structures.CaseInsensitiveDict(rr.headers));break
       else:raise RuntimeError('Public resolver rate limit')
      url,headers=parent.links[key]
     if head:
      self.send_response(200)
      for k in ['Content-Length','ETag','Last-Modified','Content-Type']:
       if k in headers:self.send_header(k,headers[k])
      self.end_headers();return
     reqheaders={'Accept-Encoding':'identity'}
     if self.headers.get('Range'):reqheaders['Range']=self.headers['Range']
     s=parent.http
     with s.get(url,headers=reqheaders,stream=True,timeout=(10,90)) as r:
      r.raise_for_status();expected=int(r.headers['Content-Length'])
      with parent.lock:
       # Conservative total across all restarted processes is accounted separately.
       assert parent.bytes+expected<8_000_000_000,'Per-source network cap'
       parent.bytes+=expected;parent.requests+=1
      self.send_response(r.status_code)
      for k in ['Content-Length','Content-Range','ETag','Last-Modified','Content-Type']:
       if k in r.headers:self.send_header(k,r.headers[k])
      self.end_headers()
      for b in r.iter_content(1024*1024):self.wfile.write(b)
    except Exception as e:
     parent.errors+=1
     try:self.send_error(503,type(e).__name__)
     except Exception:pass
  self.server=ThreadingHTTPServer(('127.0.0.1',0),Handler);self.thread=threading.Thread(target=self.server.serve_forever,daemon=True);self.thread.start()
 def options(self):return {'endpoint':f'http://127.0.0.1:{self.server.server_port}','allow_http':'true','skip_signature':'true','region':'us-east-1','virtual_hosted_style_request':'false'}
 def close(self):self.server.shutdown();self.server.server_close();self.thread.join();self.http.close()

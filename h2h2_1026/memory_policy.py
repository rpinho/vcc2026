"""Process-local heap release and proactive disk/swap reserve watchdog."""
import local_runtime as j
import ctypes,gc,os,sys,threading,time
_lib=ctypes.CDLL(None)
_relief=getattr(_lib,'malloc_zone_pressure_relief',None)
if _relief is not None:
 _relief.argtypes=[ctypes.c_void_p,ctypes.c_size_t];_relief.restype=ctypes.c_size_t
def relieve():
 if _relief is not None:_relief(None,0)
def install():
 original=gc.collect
 def collect(*a,**k):
  result=original(*a,**k);relieve();return result
 gc.collect=collect
 stop=threading.Event();minimum=[j.shutil.disk_usage(j.ROOT).free]
 def watch():
  while not stop.wait(.25):
   free=j.shutil.disk_usage(j.ROOT).free;minimum[0]=min(minimum[0],free)
   if free<30_500_000_000:
    message=f'RESOURCE STOP: free disk {free} bytes below proactive 30.5 GB decimal threshold; no continued compute.\n'
    os.write(2,message.encode())
    record=dict(utc=j.utc(),free_bytes=free,exit_code=75,argv=sys.argv,threshold_bytes=30_500_000_000)
    try:j.save(f'results/resource_stop_{os.getpid()}.json',record)
    except Exception:
     # A tiny failure receipt must survive an abrupt allocation across the
     # hard guard too. No numeric work or large storage allocation follows.
     try:(j.ROOT/f'results/resource_stop_{os.getpid()}.json').write_text(j.json.dumps(record,indent=2)+'\n')
     except Exception:pass
    os._exit(75)
 thread=threading.Thread(target=watch,name='disk-reserve-watchdog',daemon=True);thread.start()
 def finish():
  stop.set();thread.join(1);j.save(f'results/memory_policy_{os.getpid()}.json',dict(utc=j.utc(),minimum_sampled_free_bytes=minimum[0],threshold_bytes=30_500_000_000,poll_seconds=.25,process_pid=os.getpid(),argv=sys.argv,policy='process-local heap relief and 30.5 GB decimal proactive reserve under ORCHESTRATOR_AUTHORIZATION.md; existing warm cache retained; no numerical formula changes'))
 return finish

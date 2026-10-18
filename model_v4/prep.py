"""Invoke installed vcc prep with storage-only resource adapter; never submit."""
import runtime,os,sys,subprocess,pathlib
rt=runtime.ROOT;name='prediction_v4'+sys.argv[1];assert sys.argv[1] in ['a','b']
import h5py,time,signal
def stop_wrapper(signum,frame):raise KeyboardInterrupt
signal.signal(signal.SIGTERM,stop_wrapper)
first=rt/'artifacts'/f'{name}_structure_first.json'
if first.exists():assert runtime.json.loads(first.read_text())==runtime.json.loads((rt/'artifacts'/f'{name}_structure.json').read_text()),'Regenerated count matrix differs from first exact certificate'
input_path=rt/'out'/f'{name}.h5ad'
with h5py.File(input_path,'r') as h:nnz=len(h['X/data'])
needed=max(6*nnz,1.3*input_path.stat().st_size)+512*1024**2
start=time.time();last=0
while runtime.other_prep_active() or runtime.shutil.disk_usage('/').free-needed<43*1024**3:
 runtime.guard()
 if time.time()-start>1800:raise RuntimeError('Packaging capacity unavailable after 30 minutes; input retained')
 if time.time()-last>60:
  print('Waiting for packaging capacity: free GiB',round(runtime.shutil.disk_usage('/').free/1024**3,2),'required GiB',round((needed+43*1024**3)/1024**3,2),flush=True);last=time.time()
 time.sleep(5)
env=os.environ.copy();env.update(TMPDIR=str(rt/'artifacts/tmp'),PYTHONPATH=str(rt/'prep_support'),VCC_RESOURCE_LIMITS='1',VCC_MMAP_DATA='1',VCC_DATA_U16='1',VCC_DEFER_ENCODING='1',VCC_DELETE_INPUT='1',VCC_KEEP_U16='0',VCC_RESERVE_GIB='43',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OMP_NUM_THREADS='4',NUMBA_NUM_THREADS='4',POLARS_MAX_THREADS='4')
cmd=[(__import__('shutil').which('vcc') or 'vcc'),'prep',str(rt/'out'/f'{name}.h5ad'),'-g',str(rt.parent/'data/gene_names.csv'),'--perts',str(rt.parent/'data/pert_counts.csv'),'--contexts','A,B,C','-o',str(rt/'out'/f'{name}.vcc'),'--force']
p=subprocess.Popen(cmd,env=env);print('prep worker PID',p.pid,flush=True)
minimum=runtime.shutil.disk_usage('/').free
try:
 while p.poll() is None:
  free=runtime.shutil.disk_usage('/').free;minimum=min(minimum,free)
  if free<43*1024**3:
   import signal
   p.send_signal(signal.SIGINT)
   try:p.wait(timeout=10)
   except subprocess.TimeoutExpired:p.terminate();p.wait()
   raise RuntimeError('43 GiB packaging reserve reached; stopped own CLI')
  time.sleep(.1)
 code=p.returncode
 runtime.save(rt/'artifacts'/f'{name}_prep_resource.json',dict(minimum_sampled_free_bytes=minimum,sampling_seconds=.1,packaging_reserve_gib=43,cli_exit=code))
finally:
 if p.poll() is None:
  p.send_signal(signal.SIGINT)
  try:p.wait(timeout=15)
  except subprocess.TimeoutExpired:p.terminate();p.wait()
print('prep exit',code,flush=True);sys.exit(code)

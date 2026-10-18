"""Sequential resumable queue; foreground-owned children, no detached workers."""
import fcntl,json,os,signal,subprocess,sys,time
from pathlib import Path
P=Path(__file__).resolve().parent
lock=(P/'artifacts/queue_v2.lock').open('a')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
child=None

def stop(signum,frame):
    if child is not None and child.poll() is None:
        os.killpg(child.pid,signal.SIGTERM)
    raise SystemExit(128+signum)
for s in (signal.SIGTERM,signal.SIGINT):signal.signal(s,stop)
env=os.environ.copy()
for k in ('PROXY_DE_THREADS','NUMBA_NUM_THREADS'):env[k]='2'
for k in ('OMP_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','TOKIO_WORKER_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[k]='1'
jobs=[]
for dest in ('H1','CD4'):
    missing=[r for r in ('v1','v2b','v3b','v4a','v4b','v5a','v5b') if not (P/f'artifacts/scores_v2/{dest}/{r}.json').exists()]
    if missing:jobs.append((dest,','.join(missing)))
jobs.append(('HCT116',None))
invalid=set()
for dest,recipe in jobs:
    if dest in invalid:continue
    args=[sys.executable,'-u',str(P/'run_v2.py'),'--dest',dest]
    args+=['--recipes',recipe] if recipe else ['--anchors-only']
    log=P/f'logs/v2_{dest}_{"arms" if recipe else "anchors"}.log'
    start=time.time()
    with log.open('a') as out:
        child=subprocess.Popen(args,env=env,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
        print('START',dest,recipe,'pid',child.pid,flush=True)
        (P/'artifacts/queue_v2_current.json').write_text(json.dumps(dict(dest=dest,recipe=recipe,pid=child.pid,log=str(log),started=start)))
        rc=child.wait()
    print('END',dest,recipe,'rc',rc,'seconds',round(time.time()-start),flush=True)
    if rc==0 and recipe:
        with (P/f'logs/v2_{dest}_fast_prepare.log').open('a') as out:
            child=subprocess.Popen([sys.executable,'-u',str(P/'fast_components.py'),'--dest',dest,'--prepare'],env=env,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
            print('FAST_PREPARE_START',dest,'pid',child.pid,flush=True)
            fast_rc=child.wait()
        print('FAST_PREPARE_END',dest,'rc',fast_rc,flush=True)
    if rc:
        (P/'artifacts/queue_v2_failures.jsonl').open('a').write(json.dumps(dict(dest=dest,recipe=recipe,rc=rc,log=str(log)))+'\n')
        # Invalid anchors affect every arm; continue only to the next destination.
        anchor=P/f'artifacts/scores_v2/{dest}/anchors.json'
        if anchor.exists() and json.loads(anchor.read_text()).get('status')=='INVALID':
            invalid.add(dest)
child=None
print('QUEUE_COMPLETE',flush=True)

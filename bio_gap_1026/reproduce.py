"""One sequential local reproduction command; at most two numerical workers."""
import runtime as rt
rt.enforce()
import subprocess,sys,os,time,signal
def stop(signum,frame):raise KeyboardInterrupt
signal.signal(signal.SIGTERM,stop)
begin=time.monotonic();p=None
try:
    for name in ['reconstruct','render_bulk','analyze','inventory','cpu_baselines','benchmark_cpu','verify','write_report']:
        rt.guard(True)
        if time.monotonic()-begin>5*3600:raise RuntimeError('5h stage boundary; total cap 6h')
        with (rt.ROOT/f'{name}.log').open('w') as log:
            p=subprocess.Popen([sys.executable,'-B',str(rt.ROOT/f'{name}.py')],cwd=rt.ROOT,stdout=log,stderr=subprocess.STDOUT)
            rt.save('results/reproduction_running.json',dict(stage=name,pid=p.pid,controller=os.getpid()))
            code=p.wait();assert code==0,(name,code)
        p=None
    rt.status(f'Reproduction complete in {time.monotonic()-begin:.1f}s; all child processes reaped.')
finally:
    if p is not None and p.poll() is None:p.terminate();p.wait()

"""Isolate unmodified model_v6 bank readers from all predecessor writes."""
import common as job
import sys, types, shutil, time, os
from pathlib import Path

def guard(extra=0, floor=30, **kwargs):
    while shutil.disk_usage(job.ROOT).free-extra < floor*1024**3:
        print('WAIT disk reserve',floor,'GiB',flush=True); time.sleep(30)

def install():
    rt=types.ModuleType('runtime')
    rt.ROOT=job.ROOT; rt.VCC=job.REPO; rt.json=job.json; rt.sha=job.sha; rt.save=job.save
    rt.guard=guard; rt.shutil=shutil
    rt.free_gib=lambda:shutil.disk_usage(job.ROOT).free/1024**3
    rt.status=lambda *a,**kw:None
    sys.modules['runtime']=rt
    os.environ.update(XDG_CACHE_HOME=str(job.ROOT/'artifacts/cache'),HF_HOME=str(job.ROOT/'artifacts/hf'),
                      HF_HUB_DISABLE_IMPLICIT_TOKEN='1',RAYON_NUM_THREADS='2',TOKIO_WORKER_THREADS='2',ARROW_NUM_THREADS='2')
    import pyarrow
    pyarrow.set_cpu_count(2);pyarrow.set_io_thread_count(2)
    for p in ['model_v6/vendor/deps','model_v4','model_v6','model_v8','validation']:
        sys.path.append(str(job.REPO/p))
    return rt

def readers():
    install()
    import importlib.util
    spec=importlib.util.spec_from_file_location('unchanged_panel_banks',job.REPO/'model_v6/banks.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

"""Fourth numerical worker: reverse tail shards, unchanged reader, no collision with the main fetch."""
import common as job
import io_guard;io_guard.install()
import numpy as np,pandas as pd,ctypes,os,errno,time,gc,argparse
from source_runtime import readers,guard

def publish_no_replace(source,destination):
    # macOS SDK sys/stdio.h: RENAME_EXCL=4. An existing main-worker directory is NEVER replaced.
    assert source.resolve().is_relative_to(job.ROOT) and destination.resolve().is_relative_to(job.ROOT)
    libc=ctypes.CDLL('/usr/lib/libSystem.B.dylib',use_errno=True)
    rename=libc.renamex_np;rename.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_uint];rename.restype=ctypes.c_int
    if rename(os.fsencode(source),os.fsencode(destination),4)==0:return True
    code=ctypes.get_errno()
    if code in [errno.EEXIST,errno.ENOTEMPTY]:return False
    raise OSError(code,os.strerror(code))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--indices',default='7,6,5');ap.add_argument('--label',default='tail');a=ap.parse_args()
    assert a.label.replace('_','').isalnum()
    frozen=job.json.loads((job.ROOT/'artifacts/expansion_targets.json').read_text())
    for f,h in frozen['files'].items():assert job.sha(job.ROOT/f)==h
    targets=sorted(set(pd.read_csv(job.ROOT/'targets_remote.csv').target_gene.astype(str)))
    shards=[targets[i:i+200] for i in range(0,len(targets),200)];m=readers()
    state=job.ROOT/f'artifacts/cd4_{a.label}_state.json';records=[]
    for i in map(int,a.indices.split(',')):
        assert 0<=i<len(shards)
        final=job.ROOT/f'banks/shards/CD4_{i:03d}'
        if final.exists():
            records.append(dict(shard=i,skipped='main fetch already owns or completed this directory'));continue
        stage=job.ROOT/f'banks/supplement/CD4_{i:03d}';stage.mkdir(parents=True,exist_ok=True);done=stage/'DONE.json'
        if not done.exists():
            guard();start=time.time();out,info=m.cd4_banks(shards[i])
            for name,bank in out.items():m.savez(stage/f'{name}.npz',**bank)
            job.save(done,dict(source='CD4',shard=i,targets=shards[i],info=info,seconds=time.time()-start,
                upstream_banks_sha256=job.sha(job.REPO/'model_v6/banks.py'),
                hashes={p.name:job.sha(p) for p in stage.glob('*.npz')},worker='reverse-tail fourth numerical worker'))
            del out;gc.collect()
        rec=job.json.loads(done.read_text());assert rec['targets']==shards[i]
        for f,h in rec['hashes'].items():assert job.sha(stage/f)==h
        published=publish_no_replace(stage,final)
        records.append(dict(shard=i,published=published,info=rec['info'],
            note='Main sequential fetch reuses DONE.json if published; otherwise private result is retained and counted as extra transfer'))
        job.save(state,dict(status='PROGRESS',pid=os.getpid(),records=records))
        print('CD4 TAIL',i,'published',published,rec['info'],flush=True)
        total=sum(job.json.loads(p.read_text())['info'].get('bytes',0) for folder in ['banks/shards','banks/supplement'] for p in (job.ROOT/folder).glob('*/DONE.json'))
        assert total<=38_000_000_000,'Reserve 2 GB of the 40 GB budget for mock GWPS reads'
    job.save(state,dict(status='COMPLETE',pid=os.getpid(),records=records))
if __name__=='__main__':main()

"""Public reads only; resumable <=200-target shards, exact upstream range/integrity kernels."""
import common as job
import io_guard;io_guard.install()
import argparse,time,gc,numpy as np,pandas as pd
from source_runtime import readers,guard

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',required=True,choices=['HEK293T','HCT116','CD4'])
    ap.add_argument('--targets',default=str(job.ROOT/'targets_remote.csv'));ap.add_argument('--limit',type=int)
    a=ap.parse_args();m=readers()
    frozen=job.json.loads((job.ROOT/'artifacts/expansion_targets.json').read_text())
    for f,h in frozen['files'].items():assert job.sha(job.ROOT/f)==h
    targets=sorted(set(pd.read_csv(a.targets).target_gene.astype(str)))
    # Fixed response-blind shard order. Keep each under upstream 4/8 GB caps.
    shards=[targets[i:i+200] for i in range(0,len(targets),200)]
    for i,tt in enumerate(shards):
        if a.limit is not None and i>=a.limit:break
        folder=job.ROOT/f'banks/shards/{a.source}_{i:03d}';folder.mkdir(parents=True,exist_ok=True)
        done=folder/'DONE.json'
        if done.exists():
            rec=job.json.loads(done.read_text());assert rec['targets']==tt
            for f,h in rec['hashes'].items():assert job.sha(folder/f)==h
            print('REUSE',a.source,i,flush=True);continue
        guard();start=time.time()
        if a.source=='CD4':out,info=m.cd4_banks(tt)
        else:
            bank,info=m.orion_bank(a.source,tt);out={a.source:bank}
        for name,bank in out.items():m.savez(folder/f'{name}.npz',**bank)
        job.save(done,dict(source=a.source,shard=i,targets=tt,info=info,seconds=time.time()-start,
            upstream_banks_sha256=job.sha(job.REPO/'model_v6/banks.py'),
            hashes={p.name:job.sha(p) for p in folder.glob('*.npz')}))
        total=sum(job.json.loads(p.read_text())['info'].get('bytes',0) for p in (job.ROOT/'banks/shards').glob('*/DONE.json'))
        print('SHARD DONE',a.source,i,info,'TOTAL BYTES',total,flush=True)
        assert total<=40_000_000_000,'40GB public-read budget exhausted'
        # Upstream writes only an aggregate moment accumulator, never raw cells.
        accumulator=job.ROOT/f'artifacts/tmp/acc_{a.source}.npz'
        if accumulator.exists():accumulator.unlink()
        del out;gc.collect()
if __name__=='__main__':main()

"""Unchanged v6 public bank API, isolated runtime; any future final panel."""
import common as job
import io_guard;io_guard.install()
import argparse,pandas as pd,numpy as np,time
from pathlib import Path
from source_runtime import readers

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--targets',required=True);ap.add_argument('--banks',required=True)
    ap.add_argument('--sources',default='K562_GWPS_raw,HCT116,HEK293T,CD4,H1,promoters');a=ap.parse_args()
    assert Path(a.banks).resolve().is_relative_to(job.ROOT)
    frame=pd.read_csv(a.targets);targets=sorted(set(frame['target_gene' if 'target_gene' in frame else frame.columns[0]].astype(str))-{'non-targeting'})
    assert targets,'empty target panel'
    sources=a.sources.split(',');assert len(sources)==len(set(sources))
    dest=Path(a.banks).resolve();dest.mkdir(parents=True,exist_ok=True)
    request=dict(targets=targets,sources=sources);request_path=dest/'panel_request.json'
    if request_path.exists():assert job.json.loads(request_path.read_text())==request,'bank directory belongs to another panel/source request; use a fresh --banks path'
    else:job.save(request_path,request)
    m=readers();chunks=[targets[i:i+200] for i in range(0,len(targets),200)]
    # Each source/chunk is independent. Upstream APIs and range integrity checks are unchanged.
    # A completed chunk binds its target list, upstream code and every output hash.
    outputs=[];records=[];start=time.time()
    for source in sources:
        for index,tt in enumerate(chunks):
            folder=dest/f'_shards/{source}_{index:04d}';folder.mkdir(parents=True,exist_ok=True)
            done=folder/'DONE.json';code_sha=job.sha(job.REPO/'model_v6/banks.py')
            if done.exists():
                r=job.json.loads(done.read_text())
                assert r['targets']==tt and r['source']==source and r['upstream_banks_sha256']==code_sha
                for name,h in r['hashes'].items():assert job.sha(folder/name)==h
                print('REUSE FINAL BANK CHUNK',source,index,flush=True)
            else:
                info=m.build_panel_banks(tt,folder,[source])
                files=list(folder.glob('*.npz'))+list(folder.glob('promoter_pairs.csv'))
                r=dict(source=source,index=index,targets=tt,upstream_banks_sha256=code_sha,info=info,
                    hashes={p.name:job.sha(p) for p in files})
                job.save(done,r)
            outputs.append(folder);records.append(r)
    from finalize_sources import merge_bank
    names=sorted({p.name for folder in outputs for p in folder.glob('*.npz') if not p.stem.endswith('.partial')})
    for name in names:merge_bank([folder/name for folder in outputs if (folder/name).exists()],dest/name)
    promoters=[folder/'promoter_pairs.csv' for folder in outputs if (folder/'promoter_pairs.csv').exists()]
    if promoters:pd.concat([pd.read_csv(p) for p in promoters],ignore_index=True).drop_duplicates().to_csv(dest/'promoter_pairs.csv',index=False)
    coverage={};seen=set()
    for name in names:
        with np.load(dest/name) as z:tt=set(z['targets'].astype(str))
        coverage[Path(name).stem]=len(set(targets)&tt);seen|=tt
    job.save(dest/f'banks_receipt_{"_".join(sources)}.json',dict(targets=targets,source_chunks=records,
        panel_coverage=coverage,targets_without_any_source=sorted(set(targets)-seen),seconds=time.time()-start,
        hashes={name:job.sha(dest/name) for name in names},upstream_banks_sha256=job.sha(job.REPO/'model_v6/banks.py')))
    print('FINAL PANEL BANKS COMPLETE',coverage,flush=True)
if __name__=='__main__':main()

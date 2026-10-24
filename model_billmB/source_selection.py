"""Freeze completed source summaries; a completed helper may supply an in-flight duplicate shard."""
import common as job
import csv

def select():
    path=job.ROOT/'artifacts/source_selection.json'
    with (job.ROOT/'targets_remote.csv').open() as stream:targets=sorted({r['target_gene'] for r in csv.DictReader(stream)})
    signature=dict(targets_sha256=job.sha(job.ROOT/'targets_remote.csv'),upstream_banks_sha256=job.sha(job.REPO/'model_v6/banks.py'))
    frozen=job.json.loads((job.ROOT/'artifacts/expansion_targets.json').read_text())
    assert signature['targets_sha256']==frozen['files']['targets_remote.csv']
    if path.exists():
        plan=job.json.loads(path.read_text());assert plan['inputs']==signature
    else:
        records=[]
        for source in ['HEK293T','HCT116','CD4']:
            for index in range(8):
                candidates=[job.ROOT/f'banks/{scope}/{source}_{index:03d}/DONE.json' for scope in ['shards','supplement']]
                done=next((p for p in candidates if p.exists()),None);assert done is not None,(source,index)
                records.append(dict(source=source,index=index,receipt=str(done.relative_to(job.ROOT)),receipt_sha256=job.sha(done)))
        plan=dict(inputs=signature,records=records,rule='first completed primary summary, otherwise completed private helper; response-blind, fixed before training')
    result={s:[] for s in ['HEK293T','HCT116','CD4']}
    for record in plan['records']:
        done=job.ROOT/record['receipt'];assert done.resolve().is_relative_to(job.ROOT)
        assert job.sha(done)==record['receipt_sha256']
        r=job.json.loads(done.read_text());source=record['source'];index=record['index']
        assert r['source']==source and r['targets']==targets[index*200:(index+1)*200]
        assert r['upstream_banks_sha256']==signature['upstream_banks_sha256']
        required={source+'.npz'} if source!='CD4' else {f'CD4_{condition}.npz' for condition in ['Rest','Stim8hr','Stim48hr','DE_statistics']}
        assert required<=set(r['hashes'])
        for name,h in r['hashes'].items():assert job.sha(done.parent/name)==h
        result[source].append(done.parent)
    assert all(len(v)==8 for v in result.values())
    if not path.exists():job.save(path,plan)
    return result

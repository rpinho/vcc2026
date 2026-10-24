"""Release owned reproducible caches only after the caller verifies its retained result."""
import common as job

def remove(paths,receipt,reason):
    previous=job.json.loads(receipt.read_text()) if receipt.exists() else {'files':[],'reclaimed_bytes':0}
    removed=[]
    for p in dict.fromkeys(paths):
        if not p.exists():continue
        q=p.resolve();assert q.is_relative_to(job.ROOT) and q.is_file()
        assert (q.is_relative_to(job.ROOT/'artifacts') and q.suffix in ['.npz','.npy']) or q==job.ROOT/'billmB_gpu_input.tar'
        removed.append(dict(path=str(p.relative_to(job.ROOT)),bytes=p.stat().st_size));p.unlink()
    job.save(receipt,dict(status='PASS own regenerable caches only',reason=reason,
        files=previous['files']+removed,reclaimed_bytes=previous['reclaimed_bytes']+sum(r['bytes'] for r in removed)))
    print('CACHE CLEANUP',len(removed),sum(r['bytes'] for r in removed),flush=True)

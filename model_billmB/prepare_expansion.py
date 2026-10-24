"""Exact billmA commit kernel, bounded native moments; baseline rows win duplicates."""
import common as job
import io_guard;io_guard.install()
import numpy as np,pandas as pd,gc,argparse,shutil
from prepare import commit,G

def merge_response(name,newpath,oldpath):
    with np.load(newpath) as z:new={k:z[k] for k in z.files}
    if oldpath and oldpath.exists():
        with np.load(oldpath) as z:old={k:z[k] for k in z.files}
        keep=~np.isin(new['targets'],old['targets'])
        # Distinct measurement masks cannot share a bank. Store separate disjoint targets.
        if not np.array_equal(old['measured'],new['measured']):
            for k in ['response','control','targets','n_cells']:new[k]=new[k][keep]
            np.savez_compressed(newpath,**new)
            shutil.copyfile(oldpath,job.ROOT/f'data/{name}_cached.npz')
        else:
            for k in ['response','control','targets','n_cells']:new[k]=np.concatenate([old[k],new[k][keep]])
            np.savez_compressed(newpath,**new)

def native(name,family):
    folder=job.REPO/f'model_sidechain/approx/sources/{name}'
    assert (folder/'DONE.json').exists()
    meta=job.json.loads((folder/'meta.json').read_text());labels=np.array(meta['labels']);genes=np.array(meta['genes'])
    n=np.load(folder/'n_cells.npy');s=np.load(folder/'count_sum.npy',mmap_mode='r')
    lib=np.load(folder/'libsize_sum.npy');np.testing.assert_array_equal(s.sum(1),lib)
    ci=np.flatnonzero(np.isin(labels,['non-targeting','Non-Targeting','control']))
    assert len(ci)==1 and n[ci[0]]>0
    control=s[ci[0]].copy();keep=(n>=30)&~np.isin(np.arange(len(n)),ci)
    # 128-target files bound peak RAM and preserve masks without compressing a 10k x 18k block.
    rows=np.flatnonzero(keep)
    existing=set()
    if name=='K562_GWPS':
        with np.load(job.REPO/'model_billm/data/K562_GWPS_raw.npz') as z:existing=set(z['targets'])
        shutil.copyfile(job.REPO/'model_billm/data/K562_GWPS_raw.npz',job.ROOT/'data/K562_GWPS_raw.npz')
    # Prefer baseline H1 rows, then only exact newly available targets.
    if name=='H1':
        with np.load(job.REPO/'model_billm/data/H1.npz') as z:existing=set(z['targets'])
        shutil.copyfile(job.REPO/'model_billm/data/H1.npz',job.ROOT/'data/H1_cached.npz')
    rows=np.array([i for i in rows if labels[i] not in existing])
    for j in range(0,len(rows),128):
        ix=rows[j:j+128];bankname=f'{name}_full_{j//128:03d}'
        p=job.ROOT/f'data/{bankname}.npz'
        if p.exists():continue
        commit(bankname,family,labels[ix],genes,s[ix].copy(),control,n[ix],dict(
            source=str(folder.relative_to(job.REPO)),done_sha256=job.sha(folder/'DONE.json'),
            count_sum_sha256=job.sha(folder/'count_sum.npy') if j==0 else 'bound by first shard receipt',
            controls='global non-targeting pooled raw counts; NOT batch matched',
            rule='baseline cached rows retained; no duplicate target within this source screen'))
        gc.collect()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--local',action='store_true');a=ap.parse_args()
    if a.local:
        native('K562_GWPS','K562');native('H1','H1')
    manifest={}
    for p in (job.ROOT/'data').glob('*.npz'):
        with np.load(p) as z:rec=dict(rows=len(z['targets']),family=str(z['family']),measured_genes=int(z['measured'].sum()))
        manifest[p.stem]=dict(rec,sha256=job.sha(p),bytes=p.stat().st_size)
    job.save(job.ROOT/'data/manifest.json',manifest)
    print('MANIFEST rows',sum(r['rows'] for r in manifest.values()),flush=True)
if __name__=='__main__':main()

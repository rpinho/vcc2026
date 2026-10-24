"""Canonical response-bank union. Existing matched rows win; fresh raw GWPS supersedes global/essential duplicates."""
import common as job
import numpy as np
from pathlib import Path

ROWS=['response','control','targets','n_cells']

def read(path):
    with np.load(path) as z:return {k:z[k] for k in z.files}

def write(path,d,provenance):
    path=Path(path);temp=path.with_suffix('.partial.npz');np.savez_compressed(temp,**d);temp.replace(path)
    job.save(path.with_suffix('.json'),dict(rows=len(d['targets']),family=str(d['family']),sha256=job.sha(path),
        bytes=path.stat().st_size,measured_genes=int(d['measured'].sum()),provenance=provenance))

def append(dest,name,fresh):
    """One target per family, except CD4 states which remain distinct banks in one family."""
    dest=Path(dest);new=read(fresh);family=str(new['family']);seen=set()
    for p in dest.glob('*.npz'):
        with np.load(p) as z:
            if str(z['family'])!=family:continue
            if family=='CD4' and not (p.stem==name or p.stem.startswith(name+'_')):continue
            if name=='K562_GWPS_raw' and not p.stem.startswith('K562_GWPS_raw'):continue
            seen.update(z['targets'].astype(str))
    keep=~np.isin(new['targets'],list(seen))
    if not keep.any():return 0
    for k in ROWS:new[k]=new[k][keep]
    path=dest/(name+'.npz');added=len(new['targets']);provenance=dict(fresh_sha256=job.sha(fresh),rule='existing eligible response rows retained, new targets appended')
    matching=None
    for candidate in [path]+sorted(dest.glob(name+'_*.npz')):
        if not candidate.exists():continue
        with np.load(candidate) as z:
            if np.array_equal(z['genes'],new['genes']) and np.array_equal(z['measured'],new['measured']):matching=candidate;break
    if matching is not None:
        path=matching;old=read(path);provenance['previous_sha256']=job.sha(path)
        for k in ROWS:new[k]=np.concatenate([old[k],new[k]])
    elif path.exists():
        suffix=__import__('hashlib').sha256(new['measured'].tobytes()).hexdigest()[:12]
        path=dest/(name+'_panel_'+suffix+'.npz');assert not path.exists()
    write(path,new,provenance);return added

def raw_targets(folder):
    seen=set()
    for p in Path(folder).glob('K562_GWPS_raw*.npz'):
        with np.load(p) as z:seen.update(z['targets'].astype(str))
    return seen

def raw_priority(folder):
    folder=Path(folder);matched=raw_targets(folder)
    for p in [folder/'K562.npz',folder/'K562_GWPS_full.npz']:
        if not p.exists():continue
        d=read(p);keep=~np.isin(d['targets'],list(matched))
        if keep.all():continue
        previous=job.sha(p)
        for k in ROWS:d[k]=d[k][keep]
        write(p,d,dict(previous_sha256=previous,rule='batch-matched raw GWPS rows take precedence; no duplicate K562-family target'))

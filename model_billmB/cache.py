"""Fold-local features, computed once rather than on every optimizer step."""
import common as rt
import numpy as np
from pathlib import Path
from dataset import Sources

class FeatureCache:
    def __init__(self,fold='production'):
        self.folder=rt.ROOT/f'artifacts/features/{fold}';self.folder.mkdir(parents=True,exist_ok=True)
        signature={p.stem:rt.sha(p) for p in (rt.ROOT/'data').glob('*.npz')}
        complete=self.folder/'complete.json'
        if complete.exists():
            rec=rt.json.loads(complete.read_text());assert rec['source_sha256']==signature
            assert rec['dataset_sha256']==rt.sha(rt.ROOT/'dataset.py')
        else:
            src=Sources([] if fold=='production' else [fold]);rows=src.rows();n=len(rows);g=src.g
            arrays={k:np.lib.format.open_memmap(self.folder/f'{k}.npy',mode='w+',dtype=dt,shape=shape) for k,dt,shape in [
                ('x','float32',(n,g*5)),('y','float32',(n,g)),('mask','uint8',(n,g)),('ids','int64',(n,))]}
            for lo in range(0,n,32):
                batch=src.batch(rows[lo:lo+32]);hi=lo+len(batch[0])
                for key,value in zip(['x','y','mask','ids'],batch):arrays[key][lo:hi]=value
                if lo%512==0:print('FEATURE CACHE',fold,lo,'/',n,flush=True)
            for a in arrays.values():a.flush()
            family_rows={}
            for j,(b,i) in enumerate(rows):family_rows.setdefault(src.banks[b]['family'],[]).append(j)
            rec=dict(fold=fold,source_sha256=signature,dataset_sha256=rt.sha(rt.ROOT/'dataset.py'),rows=n,genes=g,family_rows=family_rows,
                excluded_families=[] if fold=='production' else [fold],row_family_counts={f:sum(len(d['targets']) for d in src.banks if d['family']==f) for f in {d['family'] for d in src.banks}})
            rt.save(complete,rec)
        self.rows=rec['rows'];self.g=rec['genes'];self.record=rec
        for k in ['x','y','mask','ids']:setattr(self,k,np.load(self.folder/f'{k}.npy',mmap_mode='r'))
    def batch(self,indices):
        return self.x[indices],self.y[indices],self.mask[indices].astype(np.float32),self.ids[indices]
    def sample_rows(self,rng,n=64):
        families=sorted(self.record['family_rows'])
        draw=rng.integers(len(families),size=n)
        selected=np.empty(n,np.int64)
        for j,f in enumerate(families):
            ix=np.flatnonzero(draw==j)
            selected[ix]=rng.choice(self.record['family_rows'][f],len(ix),replace=False)
        return selected

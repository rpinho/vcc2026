"""Check spawned count workers preserve the released kernel's seeded output."""
import runtime,multiprocessing
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from build import init_worker,emit_worker
if __name__=='__main__':
 d=np.load('artifacts/atlas_fixture.npz');q=d['template'];depth=d['depths'];items=[(0,'fixture',d['p'],d['bulk'],20260910),(1,'repeat',d['p'],d['bulk'],20260910)]
 with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn'),initializer=init_worker,initargs=(q,depth,'new_dual_dual')) as executor:r=list(executor.map(emit_worker,items))
 assert np.array_equal(r[0][0],d['counts']) and np.array_equal(r[1][0],d['counts'])
 assert all(v['row_depths_exact'] and v['mean_L1_upper']<.001 and v['bulk_L1_upper']<.001 for x,v in r)
 runtime.save('artifacts/parallel_verification.json',dict(status='VERIFIED spawned workers exactly reproduce released-kernel integer fixture; independent rational moment certificates pass',workers=2,seed=20260910));print('Parallel emission VERIFIED',flush=True)

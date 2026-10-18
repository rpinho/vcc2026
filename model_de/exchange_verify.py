"""Independent integer verifier for the equal-depth exchange contract.
Checking bucket margins proves sum_i x_ig and sum_i x_ig/depth_i are unchanged.
Uses no exchange/search implementation and accepts complete emitted matrices.
"""
import hashlib
import numpy as np

def certify(before,after):
    before=np.asarray(before);after=np.asarray(after)
    assert before.shape==after.shape and before.ndim==2
    assert np.issubdtype(before.dtype,np.integer) and np.issubdtype(after.dtype,np.integer)
    assert before.min()>=0 and after.min()>=0
    depth=before.sum(1,dtype=np.int64)
    assert np.array_equal(depth,after.sum(1,dtype=np.int64))
    assert np.all((depth>0)&(depth<=1000000))
    h=hashlib.sha256();groups=0;movable_cells=0
    for d in np.unique(depth):
        rows=np.flatnonzero(depth==d)
        a=before[rows].sum(0,dtype=np.int64);b=after[rows].sum(0,dtype=np.int64)
        assert np.array_equal(a,b),'depth-bucket gene margin changed'
        h.update(np.int64(d).tobytes());h.update(a.tobytes())
        if len(rows)>1:groups+=1;movable_cells+=len(rows)
    assert np.array_equal(before.sum(0,dtype=np.int64),after.sum(0,dtype=np.int64))
    return dict(status='PASS exact integer row and depth-bucket column margins',
        bucket_margin_sha256=h.hexdigest(),row_depth_sha256=hashlib.sha256(depth.tobytes()).hexdigest(),
        movable_groups=groups,movable_cells=movable_cells,changed_entries=int(np.count_nonzero(before!=after)),
        consequence='identical raw bulk and mean cell CPM in exact rational arithmetic')

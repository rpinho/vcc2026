"""Verify file-backed task counts/digests against the original in-memory transform."""
import os
for k in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):os.environ.setdefault(k,'1')
import hashlib,json,tempfile
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import mapped_task as m


def main():
    root=Path(__file__).resolve().parent;cases=[]
    for high,interleaved in [(True,False),(True,True),(False,True)]:
        with tempfile.TemporaryDirectory(dir=root/'v2_cache',prefix='mapped_check_') as scratch:
            p=Path(scratch);x=sparse.csr_matrix(np.random.default_rng(8).poisson(2000 if high else .3,(13,40)).astype(np.int32))
            mask=np.arange(13)%3==0 if interleaved else np.isin(np.arange(13),[5,6,7])
            obs=pd.DataFrame({'target_gene':np.where(mask,'non-targeting','G0')},index=[str(i) for i in range(13)])
            var=pd.DataFrame(index=[f'G{i}' for i in range(40)]);a=ad.AnnData(x,obs=obs,var=var)
            path=p/'source.h5ad';a.write_h5ad(path);cache=p/'cache';cache.mkdir()
            m.build(path,cache,dict(fixture=True),obs.target_gene.to_numpy(),var.index.to_numpy(str),mask)
            r=json.loads((cache/'receipt.json').read_text());arrays={k:np.load(cache/(k+'.npy'),mmap_mode='r') for k in ('data','indices','indptr')}
            n=r['nnz'];mapped=sparse.csr_matrix((arrays['data'][:n],arrays['indices'][:n],arrays['indptr']),shape=x.shape)
            expected=x.copy()
            if r['thinning_probability']<1:
                expected.data=np.random.default_rng(20261018).binomial(expected.data,r['thinning_probability']).astype(np.int32);expected.eliminate_zeros()
            assert np.array_equal(expected.data,mapped.data) and np.array_equal(expected.indices,mapped.indices) and np.array_equal(expected.indptr,mapped.indptr)
            assert r['native_control_median_umi']==float(np.median(np.asarray(x[mask].sum(1)).ravel()))
            h=hashlib.sha256()
            for v in (x.data,x.indices,x.indptr):h.update(v.tobytes())
            h.update('\0'.join(obs.target_gene).encode());h.update('\0'.join(var.index).encode())
            assert h.hexdigest()==r['base_matrix_sha256']
            cases.append(dict(thinning=high,interleaved_controls=interleaved))
    out=dict(status='PASS bitwise mapped counts and base digest',cases=cases,adapter_sha256=hashlib.sha256((root/'mapped_task.py').read_bytes()).hexdigest())
    (root/'artifacts/mapped_task_fixture.json').write_text(json.dumps(out,indent=2));print(out,flush=True)


if __name__=='__main__':main()

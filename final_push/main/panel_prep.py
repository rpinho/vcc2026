"""Installed prep and independent decoded canonical-CSR/axes/group equality, arbitrary labels."""
import common as job
import io_guard;io_guard.install()
import sys,os,argparse,tarfile,tempfile,subprocess,time,hashlib
from pathlib import Path
from source_runtime import install,guard
install();sys.path.insert(0,str(job.REPO/'model_v6'))
import numpy as np,h5py,zstandard as zstd
from scipy import sparse
from verify_counts import audit_h5,column
from pipeline import read_genes,read_targets
from fractions import Fraction as F

def groups(path,records,genes):
    with h5py.File(path) as h:
        cc=column(h['obs/context']);tt=column(h['obs/target_gene']);X=h['X']
        for i,r in enumerate(records):
            lo=i*400;assert np.all(cc[lo:lo+400]==r['context']) and np.all(tt[lo:lo+400]==r['target'])
            ptr=X['indptr'][lo:lo+401];a,b=map(int,[ptr[0],ptr[-1]])
            x=sparse.csr_matrix((X['data'][a:b],X['indices'][a:b],ptr-a),shape=(400,len(genes))).toarray()
            assert hashlib.sha256(x.astype('<i4').tobytes()).hexdigest()==r['dense_counts_sha256']

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--name',required=True);ap.add_argument('--targets',required=True)
    ap.add_argument('--contexts',required=True);ap.add_argument('--genes',required=True);ap.add_argument('--delete-mock',action='store_true');a=ap.parse_args()
    assert '/' not in a.name and '..' not in a.name
    start=time.time();out=job.ROOT/f'artifacts/panels/{a.name}';gpath=out/'full_gates.json'
    g=job.json.loads(gpath.read_text());build=job.json.loads((out/'build.json').read_text())
    labels=a.contexts.split(',');genes=read_genes(a.genes);panels=read_targets(a.targets,labels)
    assert g['status']=='PASS' and len(g['records'])==sum(map(len,panels.values()))
    assert g['inputs_sha256']==build['input_sha256']==job.sha(out/'inputs.json')
    bound=job.json.loads((out/'inputs.json').read_text());assert bound['targets']==job.sha(a.targets) and bound['genes']==job.sha(a.genes)
    assert set(bound['contexts'])==set(labels)
    for r in g['records']:
        assert r['all_gates_pass'] and all(r[k] for k in ['totals_exact','depths_exact','nonDE_exact','own_target_exact'])
        assert F(r['jk'])>=F(r['v5b_jk'])
        assert F(float(np.nextafter(r['mean_certificate']['mean_L1_upper'],np.inf)))<F(15,1000)
        assert F(float(np.nextafter(r['mean_certificate']['bulk_L1_upper'],np.inf)))<F(15,10000)
    inp=job.ROOT/f'out/{a.name}.h5ad';package=job.ROOT/f'out/{a.name}.vcc'
    assert not package.exists() and len(list(job.ROOT.rglob('*.vcc')))<2
    structure=audit_h5(inp,genes,panels,labels);assert structure==build['structure'];groups(inp,g['records'],genes)
    # The closed H5 now holds every independently checked count and remains available for prep retries.
    from cache_cleanup import remove
    caches=list(out.glob('group_*.npz'))+list(out.glob('profile_*.npz'))
    for folder in (job.ROOT/'artifacts').glob(f'ref_{a.name}_*'):caches+=list(folder.glob('*.npy'))
    if a.name=='mock_billmA':caches+=list((job.ROOT/'artifacts/mock_response_banks').glob('*.npz'))
    remove(caches,job.ROOT/f'artifacts/{a.name}_preprep_cleanup.json','Full retained H5 passed independent structure, canonical hashes and every dense-group hash; inference caches are no longer required by prep/decode')
    support=job.ROOT/'prep_support';support.mkdir(parents=True,exist_ok=True);adapter=job.REPO/'model_v7/prep_support/sitecustomize.py'
    code='import sys\nsys.path.insert(0,'+repr(str(job.ROOT))+')\nimport io_guard\nio_guard.install()\nexec(compile(open('+repr(str(adapter))+').read(),'+repr(str(adapter))+',"exec"))\n'
    (support/'sitecustomize.py').write_text(code)
    guard(max(2.5*structure['nnz'],1.3*inp.stat().st_size)+512*1024**2,floor=25)
    env=os.environ.copy();env.update(PYTHONPATH=str(support),PYTHONDONTWRITEBYTECODE='1',VCC_RESOURCE_LIMITS='1',VCC_MMAP_DATA='1',
        VCC_DATA_U16='1' if structure['max_count']<=65535 else '0',VCC_DEFER_ENCODING='1',VCC_DELETE_INPUT='0',VCC_KEEP_U16='0',VCC_RESERVE_GIB='25',
        TMPDIR=str(job.ROOT/'artifacts/tmp'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',POLARS_MAX_THREADS='2',NUMBA_NUM_THREADS='2',RAYON_NUM_THREADS='2',TOKIO_WORKER_THREADS='2')
    cmd=['vcc','prep',str(inp),'-g',str(Path(a.genes).resolve()),'--perts',str(Path(a.targets).resolve()),'--contexts',a.contexts,'-o',str(package)]
    p=subprocess.Popen(cmd,env=env);job.save(out/'prep_running.json',dict(pid=p.pid,command=cmd));print('PREP PID',p.pid,flush=True)
    code=p.wait();assert code==0
    prep_seconds=time.time()-start;guard(floor=25)
    with tarfile.open(package) as tar:
        assert set(tar.getnames())=={'meta.json','pred.h5ad.zst'}
        meta=job.json.load(tar.extractfile('meta.json'))
        assert (meta['n_obs'],meta['n_vars'],meta['nnz'])==(*structure['shape'],structure['nnz'])
        with tempfile.TemporaryDirectory(dir=job.ROOT/'artifacts/tmp') as tmp:
            decoded=Path(tmp)/'decoded.h5ad'
            with decoded.open('wb') as f,zstd.ZstdDecompressor().stream_reader(tar.extractfile('pred.h5ad.zst')) as reader:
                for block in iter(lambda:reader.read(8*1024**2),b''):guard(len(block),floor=20);f.write(block)
            after=audit_h5(decoded,genes,panels,labels);assert after==structure;groups(decoded,g['records'],genes)
    receipt=dict(status='PASS installed prep and independently decoded CSR/axes/all-group hashes',command=cmd,exit_code=code,
        structure=structure,bytes=package.stat().st_size,sha256=job.sha(package),gates_sha256=job.sha(gpath),
        prep_seconds=prep_seconds,total_seconds=time.time()-start,storage_adapter_sha256=job.sha(adapter),mock_deleted=False)
    job.save(out/'decode_receipt.json',receipt)
    if a.delete_mock:
        assert a.name.startswith('mock_'),'Deletion only permitted for explicitly named mock artifacts'
        size=package.stat().st_size;package.unlink();receipt.update(mock_deleted=True,reclaimed_bytes=size)
        job.save(out/'decode_receipt.json',receipt)
    print(receipt['status'],flush=True)
if __name__=='__main__':main()

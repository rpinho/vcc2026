"""Read-only imported model_v6/run_final.sh v6 arm, with outputs confined to this job."""
import common as job
import io_guard;io_guard.install()
import sys,os,argparse,time,subprocess,tempfile,tarfile,types
from pathlib import Path
from source_runtime import install,guard
runtime=install();sys.path.insert(0,str(job.REPO/'model_v6'))
from pipeline import build,read_genes,read_targets
from verify_counts import audit_h5
import zstandard as zstd

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['build','prep'],required=True)
    for key in ['name','targets','contexts','genes','banks']:ap.add_argument('--'+key,required=True)
    a=ap.parse_args();assert '/' not in a.name and '..' not in a.name
    labels=list(dict(kv.split('=',1) for kv in a.contexts.split(',')))
    genes=read_genes(a.genes);panels=read_targets(a.targets,labels)
    if a.stage=='build':
        # The imported builder hashes ROOT/model.py and ROOT/v5.py for provenance.
        # Resolve only those provenance names to the actual unchanged imported sources.
        def source_sha(path):
            for name in ['model.py','v5.py']:
                if Path(path)==job.ROOT/name:return job.sha(job.REPO/'model_v6'/name)
            return job.sha(path)
        runtime.sha=source_sha
        args=types.SimpleNamespace(arm='v6',out=a.name,targets=a.targets,contexts=a.contexts,genes=a.genes,
            banks=a.banks,promoters=str(Path(a.banks)/'promoter_pairs.csv'),exclude='',fill=False,k=None,
            rho=1.,workers=2,pie_dir=None,pie_prefix='ens4',pie_add=300)
        build(args);return
    start=time.time();inp=job.ROOT/f'out/{a.name}.h5ad';package=job.ROOT/f'out/{a.name}.vcc'
    before=job.json.loads((job.ROOT/f'artifacts/{a.name}_build.json').read_text())['structure']
    assert audit_h5(inp,genes,panels,labels)==before and not package.exists()
    assert len(list(job.ROOT.rglob('*.vcc')))<2
    guard(max(2.5*before['nnz'],1.3*inp.stat().st_size)+512*1024**2,floor=25)
    support=job.ROOT/'prep_support';support.mkdir(exist_ok=True);adapter=job.REPO/'model_v7/prep_support/sitecustomize.py'
    (support/'sitecustomize.py').write_text('import sys\nsys.path.insert(0,'+repr(str(job.ROOT))+')\nimport io_guard\nio_guard.install()\nexec(compile(open('+repr(str(adapter))+').read(),'+repr(str(adapter))+',"exec"))\n')
    env=os.environ.copy();env.update(PYTHONPATH=str(support),PYTHONDONTWRITEBYTECODE='1',VCC_RESOURCE_LIMITS='1',VCC_MMAP_DATA='1',
        VCC_DATA_U16='1' if before['max_count']<=65535 else '0',VCC_DEFER_ENCODING='1',VCC_DELETE_INPUT='0',VCC_KEEP_U16='0',VCC_RESERVE_GIB='25',
        TMPDIR=str(job.ROOT/'artifacts/tmp'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',POLARS_MAX_THREADS='2',NUMBA_NUM_THREADS='2')
    cmd=['vcc','prep',str(inp),'-g',str(Path(a.genes).resolve()),'--perts',str(Path(a.targets).resolve()),'--contexts',','.join(labels),'-o',str(package)]
    p=subprocess.Popen(cmd,env=env);job.save(job.ROOT/f'artifacts/{a.name}_prep_running.json',dict(pid=p.pid,command=cmd));code=p.wait();assert code==0
    with tarfile.open(package) as archive:
        assert set(archive.getnames())=={'meta.json','pred.h5ad.zst'}
        meta=job.json.load(archive.extractfile('meta.json'))
        assert (meta['n_obs'],meta['n_vars'],meta['nnz'])==(*before['shape'],before['nnz'])
        with tempfile.TemporaryDirectory(dir=job.ROOT/'artifacts/tmp') as tmp:
            decoded=Path(tmp)/'decoded.h5ad'
            with decoded.open('wb') as f,zstd.ZstdDecompressor().stream_reader(archive.extractfile('pred.h5ad.zst')) as reader:
                for block in iter(lambda:reader.read(8*1024**2),b''):guard(len(block),floor=20);f.write(block)
            after=audit_h5(decoded,genes,panels,labels)
    assert before==after
    rec=dict(status='PASS v6 arm installed prep and independent exact decoded CSR/axes equality',
        command=cmd,exit_code=code,bytes=package.stat().st_size,sha256=job.sha(package),structure=after,seconds=time.time()-start,
        method='unchanged model_v6 pipeline v6 arm; storage-only resource adapter; no Bill layout gate claim',mock_deleted=False)
    job.save(job.ROOT/f'artifacts/{a.name}_prep.json',rec)
    if a.name=='mock_smoke_v6':
        rec.update(mock_deleted=True,reclaimed_bytes=package.stat().st_size+inp.stat().st_size)
        package.unlink();inp.unlink();job.save(job.ROOT/f'artifacts/{a.name}_prep.json',rec)
    print('FALLBACK v6 VERIFIED',flush=True)
if __name__=='__main__':main()

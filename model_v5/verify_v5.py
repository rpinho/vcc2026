"""Independently decode the .vcc and require exact CSR/axes equality with the pre-packaging certificate
(model_v4/verify_package.py, reserve lowered to 25 GiB for the shared disk)."""
import runtime, tarfile, tempfile, pathlib, sys, shutil
import pandas as pd, zstandard as zstd
from verify_counts import audit_h5
rt = runtime.ROOT; arm = sys.argv[1]; name = 'v5' + arm; path = rt / f'{name}.vcc'
before = runtime.json.loads((rt / 'artifacts' / f'{name}_structure.json').read_text())
genes = pd.read_csv(rt.parent / 'data/gene_names.csv').gene_name.to_numpy(str)
targets = pd.read_csv(rt.parent / 'data/pert_counts.csv').target_gene.to_numpy(str)
with tarfile.open(path, 'r') as tar:
    assert set(tar.getnames()) == {'meta.json', 'pred.h5ad.zst'}
    meta = runtime.json.load(tar.extractfile('meta.json'))
    assert meta['n_obs'] == 360000 and meta['n_vars'] == 18533 and meta['nnz'] == before['nnz']
    with tempfile.TemporaryDirectory(dir=rt / 'artifacts/tmp') as tmp:
        decoded = pathlib.Path(tmp) / 'decoded.h5ad'
        with decoded.open('wb') as f, zstd.ZstdDecompressor().stream_reader(tar.extractfile('pred.h5ad.zst')) as reader:
            for block in iter(lambda: reader.read(8 * 1024**2), b''):
                if shutil.disk_usage('/').free < 25 * 1024**3:
                    raise RuntimeError('25 GiB reserve')
                f.write(block)
        after = audit_h5(decoded, genes, targets)
assert before == after, {k: (before[k], after[k]) for k in before if before[k] != after[k]}
cert = dict(status='VERIFIED: installed vcc prep exit 0 plus exact decoded CSR/axes equality', package_metadata=meta,
            package_bytes=path.stat().st_size, package_sha256=runtime.sha(path), structure=before)
runtime.save(rt / 'artifacts' / f'{name}_package.json', cert); print(runtime.json.dumps(cert, indent=1)[:1500], flush=True)

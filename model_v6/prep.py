"""Package with the installed official CLI (`vcc prep`), generalized: any perts file, gene list and context labels.
Same storage-only resource adapter as model_v4/v5 (prep_support/sitecustomize.py). `vcc prep` runs its full local
validation (genes, contexts, per-context official targets, 400 cells/target, integer counts, caps). Never submits.
Then independently decodes the .vcc and requires exact CSR/axes equality with the pre-packaging audit.

python prep.py NAME --perts panel.csv --contexts D,E,F [--genes gene_names.csv]
"""
import runtime, os, sys, subprocess, pathlib, time, signal, shutil, argparse, tarfile, tempfile, json
import h5py, pandas as pd, zstandard as zstd
from verify_counts import audit_h5
from pipeline import read_genes, read_targets

rt = runtime.ROOT


def prep(name, perts, contexts, genes_csv=None, reserve_gib=22, delete_input=False):
    inp = rt / 'out' / f'{name}.h5ad'; out = rt / 'out' / f'{name}.vcc'
    genes_csv = genes_csv or str(runtime.VCC / 'data/gene_names.csv')
    with h5py.File(inp, 'r') as h:
        nnz = len(h['X/data'])
    needed = max(2.5 * nnz, 1.3 * inp.stat().st_size) + 512 * 1024**2  # measured: v5 prep peaks ~2 bytes/nnz; the 20 GiB abort floor still applies
    runtime.guard(needed, floor=reserve_gib, yield_kids=False)
    (rt / 'artifacts/tmp').mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(TMPDIR=str(rt / 'artifacts/tmp'), PYTHONPATH=str(rt / 'prep_support'), VCC_RESOURCE_LIMITS='1', VCC_MMAP_DATA='1',
               VCC_DATA_U16='1', VCC_DEFER_ENCODING='1', VCC_DELETE_INPUT='1' if delete_input else '0', VCC_KEEP_U16='0', VCC_RESERVE_GIB=str(reserve_gib),
               OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', OMP_NUM_THREADS='3', NUMBA_NUM_THREADS='3', POLARS_MAX_THREADS='3')
    cmd = ['nice', '-n', '19', (shutil.which('vcc') or 'vcc'), 'prep', str(inp), '-g', genes_csv, '--perts', str(perts),
           '--contexts', contexts, '-o', str(out), '--force']
    t0 = time.time(); p = subprocess.Popen(cmd, env=env); print('prep PID', p.pid, flush=True); minimum = shutil.disk_usage('/').free
    try:
        while p.poll() is None:
            free = shutil.disk_usage('/').free; minimum = min(minimum, free)
            if free < 20 * 1024**3:
                p.send_signal(signal.SIGINT); p.wait(timeout=30); raise RuntimeError('20 GiB floor reached; stopped own CLI')
            time.sleep(0.5)
        code = p.returncode
    finally:
        if p.poll() is None:
            p.send_signal(signal.SIGINT)
            try:
                p.wait(timeout=15)
            except subprocess.TimeoutExpired:
                p.terminate(); p.wait()
    rec = dict(cli_exit=code, prep_seconds=round(time.time() - t0, 1), minimum_free_bytes=minimum, command=' '.join(cmd[2:]), output=out.name,
               bytes=out.stat().st_size if out.exists() else None)
    if code == 0:
        rec.update(verify(name, perts, contexts.split(','), genes_csv))
    runtime.save(rt / 'artifacts' / f'{name}_prep.json', rec)
    print('prep exit', code, json.dumps({k: v for k, v in rec.items() if k != 'structure'}), flush=True)
    return rec


def verify(name, perts, labels, genes_csv):
    path = rt / 'out' / f'{name}.vcc'
    before = json.loads((rt / 'artifacts' / f'{name}_build.json').read_text())['structure']
    genes = read_genes(genes_csv); panels = read_targets(perts, labels)
    t0 = time.time()
    with tarfile.open(path, 'r') as tar:
        assert set(tar.getnames()) == {'meta.json', 'pred.h5ad.zst'}
        meta = json.load(tar.extractfile('meta.json'))
        assert meta['n_obs'] == before['shape'][0] and meta['n_vars'] == len(genes) and meta['nnz'] == before['nnz']
        with tempfile.TemporaryDirectory(dir=rt / 'artifacts/tmp') as tmp:
            decoded = pathlib.Path(tmp) / 'decoded.h5ad'
            with decoded.open('wb') as f, zstd.ZstdDecompressor().stream_reader(tar.extractfile('pred.h5ad.zst')) as reader:
                for block in iter(lambda: reader.read(8 * 1024**2), b''):
                    runtime.guard(len(block), floor=20, yield_kids=False); f.write(block)
            after = audit_h5(decoded, genes, panels, labels)
    assert before == after, {k: (before[k], after[k]) for k in before if before[k] != after[k]}
    return dict(verify='VERIFIED: vcc prep exit 0 plus exact decoded CSR/axes equality', package_metadata=meta, package_sha256=runtime.sha(path),
                verify_seconds=round(time.time() - t0, 1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('name'); ap.add_argument('--perts', required=True); ap.add_argument('--contexts', required=True)
    ap.add_argument('--genes', default=None); ap.add_argument('--delete-input', action='store_true')
    a = ap.parse_args(); r = prep(a.name, a.perts, a.contexts, a.genes, delete_input=a.delete_input); sys.exit(r['cli_exit'])

"""Package with the installed official CLI (`vcc prep`, same command and storage-only adapter as model_v4/prep.py).
`vcc prep` runs the CLI's full local validation (genes, contexts, official targets, 400 cells/target, raw
integer counts, per-cell cap, nnz cap) before writing. Never submits."""
import runtime, os, sys, subprocess, pathlib, time, signal, shutil
import h5py
rt = runtime.ROOT; arm = sys.argv[1]; assert arm in ['a', 'b']; name = 'v7' + arm
RESERVE = 22 * 1024**3  # rsi loop pauses heavy steps below 20 GB free
inp = rt / 'out' / f'{name}.h5ad'
with h5py.File(inp, 'r') as h:
    nnz = len(h['X/data'])
needed = max(6 * nnz, 1.3 * inp.stat().st_size) + 512 * 1024**2
t0 = time.time()
while shutil.disk_usage('/').free - needed < RESERVE:
    if time.time() - t0 > 3600:
        raise RuntimeError('no packaging capacity after 1 h; input retained')
    print('waiting for disk: free GiB', round(shutil.disk_usage('/').free / 1024**3, 1), 'need GiB', round((needed + RESERVE) / 1024**3, 1), flush=True)
    time.sleep(60)
(rt / 'artifacts/tmp').mkdir(parents=True, exist_ok=True)
env = os.environ.copy()
env.update(TMPDIR=str(rt / 'artifacts/tmp'), PYTHONPATH=str(rt / 'prep_support'), VCC_RESOURCE_LIMITS='1', VCC_MMAP_DATA='1',
           VCC_DATA_U16='1', VCC_DEFER_ENCODING='1', VCC_DELETE_INPUT='0', VCC_KEEP_U16='0', VCC_RESERVE_GIB='22',
           OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', OMP_NUM_THREADS='1', NUMBA_NUM_THREADS='1', POLARS_MAX_THREADS='1')
out = rt / f'{name}.vcc'
cmd = ['nice', '-n', '0', str(pathlib.Path.home() / '.local/bin/vcc'), 'prep', str(inp), '-g', os.environ.get('V7_GENES',str(rt.parent / 'data/gene_names.csv')),
       '--perts', os.environ.get('V7_TARGETS',str(rt.parent / 'data/pert_counts.csv')), '--contexts', os.environ.get('V7_CONTEXTS','A,B,C'), '-o', str(out), '--force']
p = subprocess.Popen(cmd, env=env); print('prep PID', p.pid, flush=True); minimum = shutil.disk_usage('/').free
try:
    while p.poll() is None:
        free = shutil.disk_usage('/').free; minimum = min(minimum, free)
        if free < 20 * 1024**3:
            p.send_signal(signal.SIGINT); p.wait(timeout=30)
            raise RuntimeError('20 GiB floor reached; stopped own CLI')
        time.sleep(0.5)
    code = p.returncode
    runtime.save(rt / 'artifacts' / f'{name}_prep.json', dict(cli_exit=code, minimum_free_bytes=minimum, command=' '.join(cmd[3:]),
                 output=str(out.name), bytes=out.stat().st_size if out.exists() else None,
                 sha256=runtime.sha(out) if out.exists() and code == 0 else None))
finally:
    if p.poll() is None:
        p.send_signal(signal.SIGINT)
        try:
            p.wait(timeout=15)
        except subprocess.TimeoutExpired:
            p.terminate(); p.wait()
print('prep exit', code, flush=True); sys.exit(code)

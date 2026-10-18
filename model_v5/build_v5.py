"""Build v5a / v5b prediction H5AD for official contexts A/B/C (same writer and checks as model_v4/build.py).
Label-free: uses only public source banks, the official control files and our own prediction.
"""
import runtime, argparse, time, gc, json
import numpy as np, pandas as pd, anndata as ad, h5py
from scipy import sparse
from model import template, heads, emit, atlas, load_banks, aligned_bank, cd4_family
from verify_counts import moment_certificate, audit_h5
import v5

rt = runtime.ROOT
RESERVE_GIB = 25  # rsi loop pauses heavy steps below 20 GB free; stay above it


def free_gib():
    return runtime.shutil.disk_usage('/').free / 1024**3


class Writer:
    def __init__(self, path, genes):
        self.n = 0; self.nnz = 0; self.obs = []; self.f = h5py.File(path, 'w'); self.genes = genes
        self.g = self.f.create_group('X')
        self.f.attrs.update({'encoding-type': 'anndata', 'encoding-version': '0.1.0'})
        self.g.attrs.update({'encoding-type': 'csr_matrix', 'encoding-version': '0.1.0'})
        for key, dt in [('data', 'int32'), ('indices', 'int32'), ('indptr', 'int64')]:
            self.g.create_dataset(key, shape=(1 if key == 'indptr' else 0,), maxshape=(None,), dtype=dt,
                                  chunks=(4096 if key == 'indptr' else 262144,), compression='gzip',
                                  compression_opts=1, shuffle=True)
        self.g['indptr'][0] = 0

    def append(self, x, c, t):
        x = sparse.csr_matrix(x)
        while free_gib() < RESERVE_GIB:
            print('waiting for disk, free GiB', round(free_gib(), 1), flush=True); time.sleep(30)
        for key, value in [('data', x.data.astype(np.int32)), ('indices', x.indices.astype(np.int32)),
                           ('indptr', x.indptr[1:].astype(np.int64) + self.nnz)]:
            ds = self.g[key]; start = len(ds); ds.resize((start + len(value),)); ds[start:] = value
        self.nnz += x.nnz; self.n += x.shape[0]; self.obs.extend([(c, t)] * x.shape[0])

    def close(self):
        self.g.attrs['shape'] = (self.n, len(self.genes))
        obs = pd.DataFrame(self.obs, columns=['context', 'target_gene'], index=[f'pred_{i}' for i in range(self.n)])
        for k in obs:
            obs[k] = pd.Categorical(obs[k])
        ad.io.write_elem(self.f, 'obs', obs); ad.io.write_elem(self.f, 'var', pd.DataFrame(index=pd.Index(self.genes)))
        for k in ['uns', 'obsm', 'varm', 'obsp', 'varp', 'layers']:
            ad.io.write_elem(self.f, k, {})
        self.f.close()


_W = None
def init_worker(q, d):
    global _W
    _W = (q, d)


def emit_worker(item):
    """v4a's emitter unchanged: pooled control template (pool_k=4), pooled-moment fit, exact integer
    per-gene totals (AtlasShift repair_bulk_totals)."""
    i, b, seed, want_jk, sharp = item; q, d = _W
    x = emit(q, None, b, d, seed, 'anchor_single_mean').astype(np.int32)
    d = x.sum(1)
    assert d.min() > 0 and d.max() <= 1_000_000
    cert = moment_certificate(x, d, b, b)
    x = sparse.csr_matrix(x)
    jk = None
    if want_jk:
        jk = float(v5.jackknife_correction(x, np.zeros(x.shape[0], np.intp), 1, v5.TS, chunk=32)[0])
    return i, x, cert, jk


def profile(arm, targets, genes, a, c, ci, log):
    """Return the final per-target bulk profile b [T, G] for one context, plus the template."""
    q, d, mean, bulk = template(a, genes, pool_k=4, seed=20260910 + ci)  # identical to v4a
    (p, b), _ = heads(targets, genes, mean, bulk, 'anchor_single_mean')  # v4a heads, unchanged
    del p
    Lc = v5.logbulk(bulk)[None, :]; L = v5.logbulk(b)
    info = dict(context=c, B_v4a_profile=v5.budget(b, targets, genes)); mask = None
    if arm == 'b':
        V, s = v5.control_axes(a.X, k=50, seed=ci)
        L = Lc + v5.deflate(L - Lc, V)
        names = load_banks(None, True)
        src = [aligned_bank(n, targets, genes, 'log2fc') for n in names]
        cd = cd4_family(targets, genes, 'log2fc')
        if cd is not None:
            src.append(cd)
        mask, net = v5.de_call_mask(src, L - Lc, k=ARGS.k)
        del src
        L = Lc + np.where(mask, L - Lc, 0.0)
        info['calls_k'] = int(ARGS.k)
        info.update(calls_per_target_mean=float(mask.sum(1).mean()), calls_zero_targets=int((mask.sum(1) == 0).sum()),
                    deflation_k=50, control_sv_top5=s[:5].tolist())
    return q, d, L, Lc, mask, info


ARGS = None
def main(args):
    global ARGS
    ARGS = args
    genes = pd.read_csv(rt.parent / 'data/gene_names.csv').gene_name.to_numpy(str)
    targets = pd.read_csv(rt.parent / 'data/pert_counts.csv').target_gene.to_numpy(str)
    assert (len(genes), len(targets)) == (18533, 300)
    name = 'v5' + args.arm; path = rt / 'out' / f'{name}.h5ad'; assert not path.exists()
    writer = Writer(path, genes); mom = []; records = []; start = time.time()
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor
    try:
        for ci, c in enumerate('ABC'):
            a = ad.read_h5ad(rt.parent / 'data' / f'context_{c}.h5ad'); assert a.n_obs == 18400
            q, d, L, Lc, mask, info = profile(args.arm, targets, genes, a, c, ci, print)
            sh = (lambda i: None) if mask is None else (lambda i: mask[i])
            del a; gc.collect()
            rng = np.random.default_rng(ci)
            probe = np.sort(rng.choice(len(targets), 40, replace=False))
            b0 = v5.from_logbulk(L)
            with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context('spawn'),
                                     initializer=init_worker, initargs=(q, d)) as ex:
                # 1) claim of our own pinned cells (40-target probe, mean jk x sum of row weights)
                jks = [r[3] for r in ex.map(emit_worker, [(i, b0[i], atlas.seed_for(c + ':' + targets[i], 20260910), True, sh(i)) for i in probe])]
                w = v5._row_weights(len(targets), len(genes), v5.excl_cols(targets, genes))
                claim_est = float(np.mean(jks) * w.sum())
                alpha, B = v5.solve_alpha(L, targets, genes, args.rho * claim_est, center=Lc)
                if args.cap1:
                    alpha = min(alpha, 1.0); B = v5.budget(v5.from_logbulk(v5.rescale(L, alpha, Lc)), targets, genes)
                bf = v5.from_logbulk(v5.rescale(L, alpha, Lc)); del b0
                info.update(claim_est_probe=claim_est, alpha=alpha, rho=args.rho, B_target_profile=B)
                print(name, c, json.dumps(info), flush=True)
                # 2) final emission, streamed; exact jk recorded for every target
                jk_all = np.zeros(len(targets)); P = np.zeros((len(targets), len(genes)))
                items = [(i, bf[i], atlas.seed_for(c + ':' + targets[i], 20260910), True, sh(i)) for i in range(len(targets))]
                pending = {}; nxt = 0
                for j in range(len(items)):
                    while nxt < min(len(items), j + 6):
                        pending[nxt] = ex.submit(emit_worker, items[nxt]); nxt += 1
                    i, x, cert, jk = pending.pop(j).result()
                    assert x.shape[0] == 400
                    writer.append(x, c, targets[i]); jk_all[i] = jk; P[i] = np.asarray(x.sum(0)).ravel()
                    mom.append(dict(context=c, target=targets[i], **cert))
                    if i % 25 == 0:
                        print(name, c, i, 'nnz', writer.nnz, 's', round(time.time() - start), 'free', round(free_gib(), 1), flush=True)
            Pn = P / P.sum(1, keepdims=True)
            info.update(B_emitted=v5.budget(Pn, targets, genes), claim_pred_emitted=float((jk_all * w).sum()))
            info['ratio_B_over_claim_pred'] = info['B_emitted'] / info['claim_pred_emitted']
            records.append(info); print(name, c, 'done', json.dumps(info), flush=True)
            del q, L, bf, P; gc.collect()
    finally:
        writer.close()
    pd.DataFrame(mom).to_csv(rt / 'artifacts' / f'{name}_moments.csv', index=False)
    runtime.save(rt / 'artifacts' / f'{name}_build.json', dict(arm=args.arm, contexts=records, nnz=writer.nnz,
                 shape=[writer.n, len(genes)], v5_sha256=runtime.sha(rt / 'v5.py'), model_sha256=runtime.sha(rt / 'model.py'),
                 builder_sha256=runtime.sha(__file__), seconds=time.time() - start))
    verified = audit_h5(path, genes, targets); runtime.save(rt / 'artifacts' / f'{name}_structure.json', verified)
    print('VERIFIED', name, verified['nnz'], flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('arm', choices=['a', 'b']); ap.add_argument('--rho', type=float, default=1.0)
    ap.add_argument('--workers', type=int, default=3); ap.add_argument('--k', type=int, default=275)
    ap.add_argument('--cap1', action='store_true'); main(ap.parse_args())

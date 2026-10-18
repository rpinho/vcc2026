"""model_v6 panel- and context-agnostic builder for v4a / v4b / v5a / v5b (and v6 = v5b + coverage fill).

Inputs are arguments, never data/pert_counts.csv or contexts A/B/C:
  --targets   CSV with a target_gene column (one panel for all contexts) or a context,target_gene table
  --contexts  LABEL=controls.h5ad[,LABEL=...]   (each file: that context's control cells, var = gene axis)
  --genes     headerless/`gene_name` CSV of the output gene axis (default: official 18,533)
  --banks     source-bank directory (default: public_data/model_v4_signature_banks); see banks.py
  --exclude   source families to hold out (e.g. K562,CD4) for leave-one-line-out tests
Arms (all label-free; each uses only public source banks + that context's control cells):
  v4a = AtlasShift anchor heads, pooled 4-cell template, one-moment pinned emission (model_v4 a)
  v4b = IV + ash + dual summary + guide-aware promoters, dual-moment emission (model_v4 b)
  v5a = v4a + B = rho x claim pooled amplitude about the control bulk (model_v5 a)
  v5b = v5a + top-50 control-PC deflation + cross-source DE call gate (model_v5 b, uncapped gate)
  v6  = v5b + coverage_fill for targets without any source row (model.py) and the gate left open
        on those rows (they have no sources to vote)
`--fill` adds coverage_fill to any arm. Emission, seeds and all numerics follow model_v5/build_v5.py.
"""
import runtime, argparse, time, gc, json, sys
import numpy as np, pandas as pd, anndata as ad, h5py
from pathlib import Path
from scipy import sparse
import model
from model import template, heads, emit, atlas, load_banks, aligned_bank, cd4_family
from verify_counts import moment_certificate, audit_h5
import v5

rt = runtime.ROOT
OFFICIAL_GENES = runtime.VCC / 'data/gene_names.csv'
ARMS = ('v4a', 'v4b', 'v5a', 'v5b', 'v6')
VARIANT = {'v4a': 'anchor_single_mean', 'v4b': 'iv_ash_dual_guide', 'v5a': 'anchor_single_mean',
           'v5b': 'anchor_single_mean', 'v6': 'anchor_single_mean'}


def read_genes(path=None):
    p = Path(path) if path else OFFICIAL_GENES
    g = pd.read_csv(p, header=None).iloc[:, 0].astype(str).to_numpy()
    if g[0] in ('gene_name', 'gene', 'genes'):
        g = g[1:]
    return g


def read_targets(path, contexts):
    t = pd.read_csv(path)
    col = 'target_gene' if 'target_gene' in t else t.columns[0]
    t = t[t[col].astype(str) != 'non-targeting']
    if 'context' in t:
        return {c: t.loc[t.context.astype(str) == c, col].astype(str).to_numpy() for c in contexts}
    panel = t[col].astype(str).to_numpy()
    return {c: panel for c in contexts}


def load_controls(path, genes):
    a = ad.read_h5ad(path)
    if 'target_gene' in a.obs:
        lab = a.obs.target_gene.astype(str).to_numpy()
        if (lab != 'non-targeting').any():
            a = a[lab == 'non-targeting'].copy()
    names = a.var_names.to_numpy(str)
    if not np.array_equal(names, genes):
        # reorder/intersect onto the requested axis; refuse silent zero-filling of the OUTPUT axis
        ix = pd.Index(names).get_indexer(genes)
        assert (ix >= 0).all(), f'{path}: {int((ix < 0).sum())} output genes missing from controls'
        a = a[:, ix].copy()
    a.X = sparse.csr_matrix(a.X)
    return a


class Writer:
    def __init__(self, path, genes):
        self.n = 0; self.nnz = 0; self.obs = []; self.f = h5py.File(path, 'w'); self.genes = genes
        self.g = self.f.create_group('X')
        self.f.attrs.update({'encoding-type': 'anndata', 'encoding-version': '0.1.0'})
        self.g.attrs.update({'encoding-type': 'csr_matrix', 'encoding-version': '0.1.0'})
        for key, dt in [('data', 'int32'), ('indices', 'int32'), ('indptr', 'int64')]:
            self.g.create_dataset(key, shape=(1 if key == 'indptr' else 0,), maxshape=(None,), dtype=dt,
                                  chunks=(4096 if key == 'indptr' else 262144,), compression='gzip', compression_opts=1, shuffle=True)
        self.g['indptr'][0] = 0

    def append(self, x, c, t):
        x = sparse.csr_matrix(x); runtime.guard(x.nnz * 8, yield_kids=False)
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
def init_worker(q, d, variant, banks, exclude, promoters):
    global _W
    model.configure(banks, exclude, promoters)
    _W = (q, d, variant)


def emit_worker(item):
    """v4 emitters unchanged: one-moment pinned (`_mean`) or AtlasShift dual-moment; exact integer totals."""
    i, p, b, seed, want_jk = item; q, d, variant = _W
    x = emit(q, p, b, d, seed, variant).astype(np.int32)
    dd = x.sum(1); assert dd.min() > 0 and dd.max() <= 1_000_000
    cert = moment_certificate(x, dd, b if p is None else p, b)
    x = sparse.csr_matrix(x); jk = None
    if want_jk:
        jk = float(v5.jackknife_correction(x, np.zeros(x.shape[0], np.intp), 1, v5.TS, chunk=32)[0])
    return i, x, cert, jk


def load_pie(path):
    """PIE output npz (model_v6/pie): targets, genes, p_de [T, G] (float16)."""
    d = np.load(path)
    return dict(targets=d['targets'].astype(str), genes=d['genes'].astype(str), p_de=d['p_de'].astype(np.float32), file=str(path))


def pie_padding(mask, delta, targets, genes, pie, k_add=300):
    """Chen Zhang #82 pattern ('centered PIE p_de rank on padding candidates'): keep the cross-source gate, then add
    the k_add genes per target with the highest gene-centred PIE p_de among genes the gate left out and the fused
    response is nonzero (direction = fused consensus). Targets/genes PIE lacks keep the cross-source call."""
    ri = pd.Index(pie['targets']).get_indexer(targets); ci = pd.Index(pie['genes']).get_indexer(genes)
    P = np.full(mask.shape, np.nan, np.float32); rr = np.flatnonzero(ri >= 0); cc = np.flatnonzero(ci >= 0)
    P[np.ix_(rr, cc)] = pie['p_de'][np.ix_(ri[rr], ci[cc])]
    C = P - np.nanmean(P, 0, keepdims=True)
    out = mask.copy(); added = np.zeros(len(targets), int)
    for i in rr:
        cand = np.flatnonzero(~mask[i] & np.isfinite(C[i]) & (np.abs(delta[i]) > 1e-9))
        if cand.size and k_add > 0:
            top = cand[np.argsort(-C[i, cand], kind='stable')[:k_add]]; out[i, top] = True; added[i] = len(top)
    return out, dict(file=pie.get('file'), targets_with_pie=int(len(rr)), genes_on_pie_axis=int(len(cc)), k_add=int(k_add),
                     added_per_target_mean=float(added.mean()), calls_after_mean=float(out.sum(1).mean()))


class Context:
    """Everything one context needs: profile heads, gate, alpha; then a streamed emission."""

    def __init__(self, arm, targets, genes, ctrl, label, ci, fill=False, k=None, rho=1.0, workers=3, log=print,
                 embed_k=32, pie=None, pie_add=300):
        assert arm in ARMS
        self.arm = arm; self.targets = np.asarray(targets); self.genes = np.asarray(genes); self.label = label
        self.ci = ci; self.workers = workers; self.rho = rho; self.log = log
        self.fill = fill or arm == 'v6'; self.k = k or len(genes); self.variant = VARIANT[arm]
        t0 = time.time()
        q, d, mean, bulk = template(ctrl, self.genes, pool_k=4, seed=20260910 + ci)
        self.q, self.d = q, d
        info = dict(context=label, arm=arm, n_targets=len(self.targets), n_genes=len(self.genes),
                    banks=str(model.CFG['banks']), exclude=sorted(model.CFG['exclude']), sources=load_banks(None, True),
                    cd4_de=model.cd4_de_path() is not None, n_controls=int(ctrl.n_obs))
        fillcfg = None
        if self.fill:
            fillcfg = dict(embedding=model.control_embedding(ctrl.X, self.genes, k=embed_k, seed=ci))
        (p, b), cov = heads(self.targets, self.genes, mean, bulk, self.variant, fill=fillcfg, info=info)
        covered = (np.asarray(cov[0]) > .01).any(1)  # filled rows carry 1e-3 synthetic coverage
        info['targets_with_source_row'] = int(covered.sum())
        info['targets_without_source_row'] = [str(t) for t in self.targets[~covered]]
        self.p = p if arm == 'v4b' else None
        self.mask = None
        if arm in ('v4a', 'v4b'):
            self.b = b; self.alpha = None
        else:
            Lc = v5.logbulk(bulk)[None, :]; L = v5.logbulk(b)
            info['B_v4a_profile'] = v5.budget(b, self.targets, self.genes)
            if arm in ('v5b', 'v6'):
                V, s = v5.control_axes(ctrl.X, k=50, seed=ci)
                L = Lc + v5.deflate(L - Lc, V)
                src = [aligned_bank(n, self.targets, self.genes, 'log2fc') for n in load_banks(None, True)]
                cd = cd4_family(self.targets, self.genes, 'log2fc')
                if cd is not None:
                    src.append(cd)
                if src:
                    mask, net = v5.de_call_mask(src, L - Lc, k=self.k)
                else:
                    mask = np.zeros(L.shape, bool)
                del src
                if self.fill:
                    mask[~covered] = True  # no sources to vote: keep the filled profile
                if pie is not None:
                    mask, pinfo = pie_padding(mask, L - Lc, self.targets, self.genes, pie, pie_add)
                    info['pie'] = pinfo
                L = Lc + np.where(mask, L - Lc, 0.0)
                self.mask = mask
                info.update(calls_k=int(self.k), calls_per_target_mean=float(mask.sum(1).mean()),
                            calls_zero_targets=int((mask.sum(1) == 0).sum()), deflation_k=50, control_sv_top5=s[:5].tolist())
            self.L, self.Lc = L, Lc
            self.b = v5.from_logbulk(L)
        self.info = info
        info['profile_seconds'] = round(time.time() - t0, 1)
        log(label, arm, 'profile', json.dumps({k: v for k, v in info.items() if k not in ('targets_without_source_row', 'coverage_fill')}), flush=True)

    def seed(self, t):
        return atlas.seed_for(self.label + ':' + t, 20260910)

    def stream(self):
        """Yield (target, csr 400 x G int32) in panel order; calibrates alpha first for v5*/v6."""
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor
        T = len(self.targets); w = v5._row_weights(T, len(self.genes), v5.excl_cols(self.targets, self.genes))
        jk_all = np.zeros(T); P = np.zeros((T, len(self.genes))); self.moments = []
        with ProcessPoolExecutor(max_workers=self.workers, mp_context=multiprocessing.get_context('spawn'), initializer=init_worker,
                                 initargs=(self.q, self.d, self.variant, model.CFG['banks'], model.CFG['exclude'], model.CFG['promoters'])) as ex:
            bf = self.b
            if self.arm not in ('v4a', 'v4b'):
                rng = np.random.default_rng(self.ci)
                probe = np.sort(rng.choice(T, min(40, T), replace=False))
                jks = [r[3] for r in ex.map(emit_worker, [(i, None, self.b[i], self.seed(self.targets[i]), True) for i in probe])]
                claim_est = float(np.mean(jks) * w.sum())
                alpha, B = v5.solve_alpha(self.L, self.targets, self.genes, self.rho * claim_est, center=self.Lc)
                bf = v5.from_logbulk(v5.rescale(self.L, alpha, self.Lc)); self.alpha = alpha
                self.info.update(claim_est_probe=claim_est, alpha=alpha, rho=self.rho, B_target_profile=B)
                self.log(self.label, self.arm, 'alpha', alpha, 'B', B, 'claim', claim_est, flush=True)
            items = [(i, None if self.p is None else self.p[i], bf[i], self.seed(self.targets[i]), True) for i in range(T)]
            pending = {}; nxt = 0
            for j in range(T):
                while nxt < min(T, j + 6):
                    pending[nxt] = ex.submit(emit_worker, items[nxt]); nxt += 1
                i, x, cert, jk = pending.pop(j).result()
                assert x.shape[0] == 400
                jk_all[i] = jk; P[i] = np.asarray(x.sum(0)).ravel()
                self.moments.append(dict(context=self.label, target=self.targets[i], **cert))
                yield self.targets[i], x
        Pn = P / P.sum(1, keepdims=True)
        self.info.update(B_emitted=v5.budget(Pn, self.targets, self.genes), claim_pred_emitted=float((jk_all * w).sum()))
        self.info['ratio_B_over_claim_pred'] = self.info['B_emitted'] / max(self.info['claim_pred_emitted'], 1e-300)


def build(args):
    t_start = time.time()
    contexts = dict(kv.split('=', 1) for kv in args.contexts.split(','))
    labels = list(contexts)
    genes = read_genes(args.genes)
    panels = read_targets(args.targets, labels)
    model.configure(args.banks, [e for e in (args.exclude or '').split(',') if e], args.promoters)
    name = args.out; path = rt / 'out' / f'{name}.h5ad'
    if path.exists():
        path.unlink()
    writer = Writer(path, genes); records = []; moments = []
    try:
        for ci, c in enumerate(labels):
            ctrl = load_controls(contexts[c], genes)
            pie = load_pie(Path(args.pie_dir) / f'{args.pie_prefix}_{c}.npz') if args.pie_dir else None
            ctx = Context(args.arm, panels[c], genes, ctrl, c, ci, fill=args.fill, k=args.k, rho=args.rho, workers=args.workers, pie=pie, pie_add=args.pie_add)
            del ctrl; gc.collect()
            for n, (t, x) in enumerate(ctx.stream()):
                writer.append(x, c, t)
                if n % 50 == 0:
                    print(name, c, n, 'nnz', writer.nnz, 's', round(time.time() - t_start), 'free', round(runtime.free_gib(), 1), flush=True)
            ctx.info['seconds_cumulative'] = round(time.time() - t_start, 1)
            records.append(ctx.info); moments += ctx.moments
            print(name, c, 'done', json.dumps({k: v for k, v in ctx.info.items() if k not in ('targets_without_source_row', 'coverage_fill')}), flush=True)
            del ctx; gc.collect()
    finally:
        writer.close()
    pd.DataFrame(moments).to_csv(rt / 'artifacts' / f'{name}_moments.csv', index=False)
    verified = audit_h5(path, genes, panels, labels)
    rec = dict(arm=args.arm, out=str(path), contexts=records, nnz=writer.nnz, shape=[writer.n, len(genes)],
               targets_file=str(args.targets), contexts_files=contexts, banks=str(model.CFG['banks']),
               exclude=sorted(model.CFG['exclude']), fill=bool(args.fill or args.arm == 'v6'),
               model_sha256=runtime.sha(rt / 'model.py'), v5_sha256=runtime.sha(rt / 'v5.py'), builder_sha256=runtime.sha(__file__),
               build_seconds=round(time.time() - t_start, 1), structure=verified)
    runtime.save(rt / 'artifacts' / f'{name}_build.json', rec)
    print('VERIFIED', name, verified['nnz'], 'build seconds', rec['build_seconds'], flush=True)
    return rec


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--arm', choices=ARMS, required=True)
    ap.add_argument('--targets', required=True)
    ap.add_argument('--contexts', required=True)
    ap.add_argument('--genes', default=None)
    ap.add_argument('--banks', default=None)
    ap.add_argument('--promoters', default=None)
    ap.add_argument('--exclude', default='')
    ap.add_argument('--out', required=True)
    ap.add_argument('--fill', action='store_true')
    ap.add_argument('--k', type=int, default=None)
    ap.add_argument('--rho', type=float, default=1.0)
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--pie-dir', default=None, help='dir with <prefix>_<context>.npz PIE outputs (v5b/v6 arms)')
    ap.add_argument('--pie-prefix', default='ens4')
    ap.add_argument('--pie-add', type=int, default=300)
    build(ap.parse_args())

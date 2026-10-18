"""Local check on the only real-cell held-out split on disk: public H1 native fold (11 targets,
400 real cells each, 2,000 input / 2,000 reference controls; H1 source bank excluded).
Official cell-eval2 0.18.0 compute_metrics; scaled with model_v4's cached baseline/replicate anchors
(same file, sha-checked), exactly as model_v4/evaluate.py scored v4a.
"""
import runtime, argparse, json, time, sys
from dataclasses import replace
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
from cell_eval2 import EvalConfig, compute_metrics, aggregate_metrics
from cell_eval2.config import DEParams
from model import heads, template, emit, atlas, load_banks, aligned_bank, cd4_family
from verify_counts import scale_record
import v5

rt = runtime.ROOT; NT = 'non-targeting'
SIX = ['pds_cosine', 'expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae',
       'de_wilcoxon_direction_fidelity_yield_raw', 'de_wilcoxon_direction_reach_raw', 'de_wilcoxon_sig_jaccard']
FOLD = 'H1'
EXACT = False
DUALP = None
DUALFAIL = []
REFBULK = None
OUT = rt / 'artifacts/evaluation_native'


def thin(a, seed):
    """model_v4/evaluate.py thin(): binomial thinning to <=20k UMI/cell (official contexts ~20k)."""
    a = a.copy(); x = a.X.tocsr(); d = np.asarray(x.sum(1)).ravel(); rng = np.random.default_rng(seed)
    p = np.minimum(1, 20000 / d)
    x.data = rng.binomial(x.data.astype(np.int64), np.repeat(p, np.diff(x.indptr))).astype(np.int32)
    x.eliminate_zeros(); a.X = x; return a


def identity(a):
    import hashlib
    h = hashlib.sha256()
    for v in [a.X.data, a.X.indices, a.X.indptr]:
        h.update(v.tobytes())
    h.update('\0'.join(a.var_names).encode()); h.update('\0'.join(a.obs_names).encode()); return h.hexdigest()


def load_fold(view='native'):
    global OUT
    a = ad.read_h5ad(rt.parent / 'public_data/v2_diagnostic/H1_balanced.h5ad')
    role = a.obs.evaluation_role.astype(str)
    ctrl = a[role == 'input_control'].copy(); real = a[role != 'input_control'].copy(); del a
    if view == 'native':
        anchors = json.loads((OUT / 'anchors.json').read_text())
        assert anchors['input_arrays_sha256'] == identity(ctrl) and anchors['reference_arrays_sha256'] == identity(real)
        return ctrl, real, anchors
    ctrl = thin(ctrl, 816); real = thin(real, 815)
    prev = rt.parent / 'model_v2/artifacts/balanced_h1'
    anchors = json.loads((prev / 'anchors.json').read_text())
    assert anchors['input_arrays_sha256'] == identity(ctrl) and anchors['reference_arrays_sha256'] == identity(real)
    anchors['tiled_baseline'] = json.loads((prev / 'baseline_audit.json').read_text())['tile_oracle']
    OUT = rt / 'artifacts/evaluation_depth20k'; OUT.mkdir(exist_ok=True)
    return ctrl, real, anchors


def build(variant, targets, genes, ctrl, base_b, src, V, rho=1.0, k_calls=275, min_net=1, base_p=None):
    """Return (pred cells, codes, jk, info) for one variant.
    '_shape' variants keep the FULL (alpha-calibrated) profile for the pinned totals and choose each
    gene's cell-level shape: called genes rank-maximising, all others rank-neutral (v5.shaped_emit)."""
    pool_k = 4 if variant == 'v4a' else 1
    if '_ipf' in variant:
        pool_k = int(variant.split('_p')[-1].split('_')[0]) if '_ipf_p' in variant else 4
    q, d, mean, bulk = template(ctrl, genes, pool_k=pool_k, seed=20260910)
    shape = '_shape' in variant
    global EXACT, DUALP
    EXACT = '_ipfx' in variant
    DUALP = base_p if '_dual' in variant else None
    mode = 'ipf' if (variant == 'v4a' or '_ipf' in variant) else ('shape' if shape else 'pin')
    global _CTRL
    _CTRL = sparse.csr_matrix(ctrl.X)
    Lc = v5.logbulk(bulk)[None, :]
    L = v5.logbulk(base_b)
    info = dict(variant=variant, pool_k=pool_k, emission=mode); sharp_mask = None
    if 'deflate' in variant or variant.startswith('v5b'):
        L = Lc + v5.deflate(L - Lc, V)
    if 'calls' in variant or variant.startswith('v5b') or '_mag' in variant:
        delta = L - Lc
        if '_mag' in variant:
            mask = np.zeros(delta.shape, bool)
            for i in range(len(targets)):
                mask[i, np.argsort(-np.abs(delta[i]))[:k_calls]] = True
        else:
            mask, net = v5.de_call_mask(src, delta, k=k_calls, min_net=min_net)
        if '_dual' in variant and DUALP is not None:
            # call-set control on the DE side only: CPM moment keeps v4a's effect on called genes and
            # the control CPM elsewhere; the pooled (MSE) moment keeps the full profile.
            pm = np.where(mask, DUALP, mean[None, :]); DUALP = pm / pm.sum(1, keepdims=True)
        elif not shape:
            L = Lc + np.where(mask, delta, 0.0)
        info['sharp'] = 'nosharp' not in variant
        sharp_mask = mask if info['sharp'] else None
        info['calls_per_target'] = mask.sum(1).tolist()
    center = None if 'mcenter' in variant else Lc
    if '_refc' in variant:
        # Officially the released controls ARE the scorer's reference controls (manifest: 18,400 control
        # cells in ground truth), so our centre equals the scorer's control. Locally input and reference
        # controls are disjoint; centring on the reference bulk reproduces the official situation.
        center = v5.logbulk(REFBULK)[None, :]
    b = v5.from_logbulk(L)
    info['B_before'] = v5.budget(b, targets, genes)
    if variant.startswith('v5') and 'noalpha' not in variant:
        # claim of realistic pinned cells (minimal-change emission ~ reference-level scatter)
        X, codes = emit_all(q, b, d, targets, 'pin' if shape else mode, None if shape else sharp_mask)
        claim, jk = v5.claim_from_cells(X, codes, len(targets), targets, genes)
        alpha, B = v5.solve_alpha(L, targets, genes, rho * claim, center=center)
        if '_cap1' in variant:
            alpha = min(alpha, 1.0)
        b = v5.from_logbulk(v5.rescale(L, alpha, center))
        info.update(claim_pred_initial=claim, alpha=alpha, rho=rho)
    X, codes = emit_all(q, b, d, targets, mode, sharp_mask)
    claim, jk = v5.claim_from_cells(X, codes, len(targets), targets, genes)
    info.update(B_final=v5.budget(pseudobulk(X, codes, len(targets)), targets, genes), claim_pred=claim)
    return X, codes, jk, info


def emit_all(q, b, d, targets, mode='ipf', sharp=None):
    chunks = []
    for i, t in enumerate(targets):
        seed = atlas.seed_for(FOLD + ':' + t, 20260910)
        if mode == 'ipf':
            if DUALP is not None:
                try:   # AtlasShift dual-summary emission: CPM moment = v4a's unscaled CPM head (DE/LFC),
                       # pooled moment = alpha-scaled bulk (MSE). Fallback to pooled-only fit on failure.
                    x = emit(q, DUALP[i], b[i], d, seed, 'anchor_single_dual')
                except ValueError:
                    DUALFAIL.append(t); x = emit(q, None, b[i], d, seed, 'anchor_single_mean')
            else:
                x = emit(q, None, b[i], d, seed, 'anchor_single_mean')
            if EXACT:
                x = v5.pin_emit(x, b[i], np.random.default_rng(seed + 1))
        else:
            rng = np.random.default_rng(seed)
            rows = np.sort(rng.choice(_CTRL.shape[0], 400, replace=False))
            if mode == 'shape':
                x = v5.shaped_emit(_CTRL[rows].toarray(), b[i], rng, called=None if sharp is None else sharp[i])
            else:
                x = v5.pin_emit(_CTRL[rows].toarray(), b[i], rng, sharp=None if sharp is None else sharp[i])
        chunks.append(sparse.csr_matrix(x))
    X = sparse.vstack(chunks, format='csr')
    codes = np.repeat(np.arange(len(targets)), [c.shape[0] for c in chunks])
    return X, codes


def pseudobulk(X, codes, n):
    S = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(n, len(codes)))
    P = np.asarray((S @ X).todense(), np.float64)
    return P / P.sum(1, keepdims=True)


def main(args):
    ctrl, real, anchors = load_fold(args.view)
    targets = np.array(anchors['targets']); genes = real.var_names.to_numpy()
    labels = real.obs.target_gene.astype(str).to_numpy()
    base, rep = anchors['tiled_baseline'], anchors['replicate']
    # v4a heads, H1 family excluded (identical call to model_v4/evaluate.py)
    _, _, mean4, bulk4 = template(ctrl, genes, pool_k=4, seed=20260910)
    (p, b), cov = heads(targets, genes, mean4, bulk4, 'anchor_single_mean', exclude=FOLD)
    names = load_banks(FOLD, True)
    src = [aligned_bank(n, targets, genes, 'log2fc') for n in names]
    cd = cd4_family(targets, genes, 'log2fc')
    if cd is not None:
        src.append(cd)
    V, s = v5.control_axes(ctrl.X, k=50)
    # real-side jackknife (diagnostic only; never used to build a prediction)
    rl = np.array([np.flatnonzero(labels == t) for t in targets], dtype=object)
    rcodes = np.concatenate([np.full(len(r), i) for i, r in enumerate(rl)])
    Xr = real.X[np.concatenate(rl)]
    claim_real_max, jk_real = v5.claim_from_cells(Xr, rcodes, len(targets), targets, genes)
    global REFBULK
    rb = np.asarray(real.X[labels == NT].sum(0)).ravel().astype(np.float64); REFBULK = rb / rb.sum()
    from evaluate_v5 import pseudobulk as _pb
    print('REAL B/claim_max', v5.budget(_pb(Xr, rcodes, len(targets)), targets, genes) / claim_real_max, flush=True)
    cfg = replace(EvalConfig.from_preset('vcc2026'), pert_col='target_gene', device='cpu', num_threads=2,
                  gather_threads=2, de=DEParams(backend='pdex'), cache_real=str(OUT / 'real_cache_v5'))
    w = v5._row_weights(len(targets), len(genes), v5.excl_cols(targets, genes))
    rows = []
    for variant in args.variants.split(','):
        rho = float(variant.split('_rho')[1]) if '_rho' in variant else 1.0
        k = int(variant.split('_k')[1].split('_')[0]) if '_k' in variant else 275
        t0 = time.time()
        mn = int(variant.split('_net')[1].split('_')[0]) if '_net' in variant else 1
        DUALFAIL.clear(); X, codes, jk, info = build(variant, targets, genes, ctrl, b, src, V, rho=rho, k_calls=k, min_net=mn, base_p=p)
        info['dual_fallbacks'] = sorted(set(DUALFAIL))
        capped = float((np.minimum(jk, jk_real) * w).sum())
        info.update(claim_capped_by_real=capped, claim_real_max=claim_real_max,
                    r_effective=min(1.0, info['B_final'] / capped) if capped > 0 else None)
        nt = real.X[labels == NT]
        pred = ad.AnnData(sparse.vstack([X, nt], format='csr'),
                          obs=pd.DataFrame({'target_gene': list(targets[codes]) + [NT] * nt.shape[0]},
                                           index=[f'p{i}' for i in range(X.shape[0] + nt.shape[0])]),
                          var=real.var.copy())
        raw = aggregate_metrics(compute_metrics(pred, real, config=cfg)).to_pandas().set_index('metric')['mean'].to_dict()
        scaled, overall = scale_record(raw, base, rep, SIX)
        rec = dict(info, raw=raw, scaled={m: float(v) for m, v in scaled.items()}, overall=float(overall),
                   seconds=time.time() - t0)
        runtime.save(OUT / f'v5_{variant}.json', rec)
        rows.append(dict(variant=variant, overall=float(overall), raw_mse_norm=raw[SIX[1]],
                         **{m: float(v) for m, v in scaled.items()},
                         alpha=info.get('alpha'), B=info['B_final'], claim_pred=info['claim_pred'],
                         claim_capped=capped, r=info['r_effective']))
        print(json.dumps(rows[-1]), flush=True)
    df = pd.DataFrame(rows); print(df.to_string(index=False))
    df.to_csv(OUT / f'v5_scores_{int(time.time())}.csv', index=False)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--variants', default='v4a,v5a_noalpha,v5a'); ap.add_argument('--view', default='native', choices=['native', 'depth20k'])
    main(ap.parse_args())

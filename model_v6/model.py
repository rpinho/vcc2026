"""model_v6: v4 heads (model_v4/model.py, unchanged numerics) made panel- and context-agnostic.

What changed vs model_v5/model.py (a byte copy of model_v4/model.py):
  * the bank directory, promoter-pair table and excluded source families are runtime settings
    (`configure`), not module constants; nothing reads data/pert_counts.csv or assumes A/B/C;
  * every bank / the CD4 DE table may be absent (a destination family held out, or a panel with no
    rows in that source): heads degrade to the remaining sources instead of failing;
  * `coverage_fill` fills targets that no source measures (fallback chain, see that function).
For the official panel with model_v4 banks and no exclusions, heads() is identical to v4a/v4b
(checked by verify_equivalence.py).
AtlasShift and Sidechain kernels are MIT; see vendor/*_LICENSE.
"""
import runtime, os
import numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse
from vendor import atlas_model as atlas
from vendor.adaptive_shrink import adaptive_shrink

WEIGHTS = {'K562_GWPS_raw': 2, 'K562_GWPS_bulk': 2, 'H1': 2, 'HCT116': 1, 'HEK293T': 1, 'HepG2': 1, 'Jurkat': 1,
           'K562': 1, 'RPE1': 1, 'CD4': .5}
FAMILY = {'K562_GWPS_raw': 'K562', 'K562_GWPS_bulk': 'K562', 'K562': 'K562', 'H1': 'H1', 'HCT116': 'HCT116',
          'HEK293T': 'HEK293T', 'HepG2': 'HepG2', 'Jurkat': 'Jurkat', 'RPE1': 'RPE1'}
SHARED_BANKS = runtime.VCC / 'public_data/model_v4_signature_banks'
CFG = dict(banks=SHARED_BANKS, exclude=frozenset(), promoters=runtime.ROOT / 'artifacts/promoter_pairs_v4.csv')


def configure(banks=None, exclude=(), promoters=None):
    CFG['banks'] = Path(banks) if banks else SHARED_BANKS
    CFG['exclude'] = frozenset(exclude)
    if promoters:
        CFG['promoters'] = Path(promoters)
    _POOL_CACHE.clear()


def family(name):
    return 'CD4' if name.startswith('CD4') else FAMILY.get(name, name.split('_')[0])


def allowed(name):
    return family(name) not in CFG['exclude']


def bank(name):
    return CFG['banks'] / f'{name}.npz'


def load_banks(exclude=None, expanded=True, validation_extra=False):
    B = CFG['banks']
    raw = bank('K562_GWPS_raw').exists() and allowed('K562_GWPS_raw')
    names = ['K562_GWPS_raw' if raw else 'K562_GWPS_bulk', 'H1']
    if expanded:
        names += ['HCT116', 'HEK293T']
    if validation_extra:
        names += ['HepG2', 'Jurkat', 'RPE1', 'K562']
    return [n for n in names if bank(n).exists() and allowed(n) and not (exclude and (n == exclude or n.startswith(exclude + '_')))]


def cd4_de_path():
    p = bank('CD4_DE_statistics')
    return p if p.exists() and allowed('CD4') else None


def aligned_bank(name, targets, genes, space, center=True):
    path = bank(name); src = atlas.load_xatlas(path, prior_counts=100000)
    full, mask = atlas.aligned_effect(src, src.targets, genes, space=space, common_subtract=int(center))
    with np.load(path) as d:
        dv = d['delta_cpm_variance'].astype(float)
        ctrl = d['matched_control_mean_cpm'].astype(float); pert = d['target_mean_cpm'].astype(float)
        pseudocount = 1 if space == 'log2fc' else 20
        if 'target_mean_cpm_variance' in d:
            var = d['target_mean_cpm_variance'].astype(float) / (pert + pseudocount)**2 + d['matched_control_mean_cpm_variance'].astype(float) / (ctrl + pseudocount)**2
        else:
            var = dv / np.maximum(np.minimum(ctrl, pert) + pseudocount, 1e-12)**2
        if space == 'log2fc':
            var /= np.log(2)**2
        var = np.maximum(var, 1e-6)
    row = pd.Index(src.targets).get_indexer(targets); col = pd.Index(src.genes).get_indexer(genes)
    r = np.flatnonzero(row >= 0); c = np.flatnonzero(col >= 0)
    eff = np.zeros((len(targets), len(genes)), np.float32); m = np.zeros_like(eff, bool); v = np.full(eff.shape, np.inf)
    eff[r] = full[row[r]]; m[r] = mask[row[r]]; v[np.ix_(r, c)] = var[np.ix_(row[r], col[c])]; v[~m] = np.inf
    return eff, m, v


def cd4_family(targets, genes, space, center=True):
    if not allowed('CD4'):
        return None
    names = [p.stem for p in CFG['banks'].glob('CD4_*.npz') if p.stem != 'CD4_DE_statistics']
    values = []; masks = []; variances = []
    for name in sorted(names):
        a, m, v = aligned_bank(name, targets, genes, space, center); values.append(a); masks.append(m); variances.append(v)
    if not values:
        return None
    a = np.stack(values); m = np.stack(masks); v = np.stack(variances); n = m.sum(0)
    effect = np.divide(np.where(m, a, 0).sum(0), n, out=np.zeros(a.shape[1:]), where=n > 0)
    var = np.divide(np.where(m, v, 0).sum(0), n, out=np.full(a.shape[1:], np.inf), where=n > 0)
    path = cd4_de_path()
    if space == 'log2fc' and path is not None:
        de, ok = atlas.aligned_cd4(path, targets, genes, common_subtract=int(center), center_scope='source')
        with np.load(path) as d:
            ri = pd.Index(d['targets']).get_indexer(targets); ci = pd.Index(d['genes']).get_indexer(genes)
            rr = np.flatnonzero(ri >= 0); cc = np.flatnonzero(ci >= 0)
            vv = d['lfcSE'].astype(float)**2; av = d['available'] & d['quality_pass'] & (d['n_cells'] >= 20)
            num = np.where(av[:, :, None], vv, 0).sum(0); den = av.sum(0)[:, None]
            vs = np.divide(num, den, out=np.full(num.shape, np.inf), where=den > 0)
            aligned = np.full(effect.shape, np.inf); aligned[np.ix_(rr, cc)] = vs[np.ix_(ri[rr], ci[cc])]
        effect = np.where(ok, de, effect); var = np.where(ok, np.maximum(aligned, 1e-6), var); n = np.maximum(n, ok)
    return effect, n > 0, var


def _ash_task(item):
    x, var = item; stats = {}; out = adaptive_shrink(x, var, stats=stats); return out, stats


def _pool(targets, genes, names, space, method='atlas', shrink='none', use_cd4=True, stats=None):
    entries = [(n, aligned_bank(n, targets, genes, space)) for n in names]
    cd = cd4_family(targets, genes, space) if use_cd4 else None
    if cd is not None:
        entries.append(('CD4', cd))
    num = np.zeros((len(targets), len(genes))); den = np.zeros_like(num); coverage = np.zeros_like(num)
    for name, (effect, mask, var) in entries:
        if method == 'iv':
            w = np.where(mask & np.isfinite(var), 1 / np.maximum(var, 1e-6), 0)
        else:
            w = WEIGHTS.get(name, 1) * mask
        num += w * effect; den += w; coverage += mask
    result = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    variance = np.divide(1, den, out=np.full_like(den, np.inf), where=den > 0)
    if shrink != 'none':
        assert method == 'iv'
        workers = min(3, int(os.environ.get('VCC_ASH_WORKERS', '1')))
        items = []
        for i, t in enumerate(targets):
            usable = (den[i] > 0) & (genes != t); items.append((result[i, usable], variance[i, usable]))
        if shrink == 'ash' and workers > 1 and len(targets) > 32:
            import multiprocessing
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as ex:
                outs = list(ex.map(_ash_task, items))
        else:
            outs = []
            for x, v in items:
                if shrink == 'ash':
                    outs.append(_ash_task((x, v)))
                elif shrink == 'js':
                    k = len(x); energy = float(np.sum(x**2)); noise = float(np.sum(v))
                    outs.append((x * max(0, 1 - (max(k - 2, 0) / max(k, 1)) * noise / max(energy, 1e-30)), {}))
                else:
                    raise ValueError(shrink)
        for i, (values, local) in enumerate(outs):
            usable = (den[i] > 0) & (genes != targets[i])
            if usable.any():
                result[i, usable] = values
            if stats is not None:
                for k, value in local.items():
                    stats[k] = stats.get(k, 0) + value
    return result, coverage, variance


_POOL_CACHE = {}
def pool(targets, genes, names, space, method='atlas', shrink='none', use_cd4=True, stats=None):
    key = (tuple(targets), tuple(genes), tuple(names), space, method, shrink, use_cd4, str(CFG['banks']), CFG['exclude'])
    if key not in _POOL_CACHE:
        if len(_POOL_CACHE) >= 2:
            _POOL_CACHE.clear()
        _POOL_CACHE[key] = _pool(targets, genes, names, space, method, shrink, use_cd4, stats)
    return _POOL_CACHE[key]


def desired_cpm(control, effects, amplitude=.6):
    effect = np.clip(amplitude * effects, -3, 3); p = np.maximum((control[None, :] * 1e6 + 1) * np.exp2(effect) - 1, 0)
    p /= p.sum(1, keepdims=True); return p


def anchor_effects(targets, genes, mean, bulk, single=True, restricted=False):
    """v4 'anchor*' fusion: AtlasShift source-centred fusion + CD4 family (weight .5), log2fc space.
    Returns the fused effect per head and its coverage (denominator) before desired_mean."""
    names = load_banks(None, True)
    sources = [atlas.load_xatlas(bank(n), prior_counts=100000) for n in names]
    spaces = ['log2fc'] if single else ['log2fc', 'bulk_delta']
    effects = {}
    for s in spaces:
        if sources:
            effects[s] = atlas.fuse_source_centered(sources, [WEIGHTS[n] for n in names], targets, genes, space=s, common_subtract=1)
        else:
            effects[s] = (np.zeros((len(targets), len(genes)), np.float32), np.zeros((len(targets), len(genes)), np.float32))
    path = cd4_de_path()
    if path is not None:
        cd, available = atlas.aligned_cd4(path, targets, genes, common_subtract=1, center_scope='source')
    out = []
    for space, control in [('log2fc', mean), ('bulk_delta', bulk)]:
        sp = 'log2fc' if single else space
        e, m = effects[sp]
        if not restricted and path is not None:
            e, m = atlas.add_cd4_family(e, m, cd, available, .5, space=sp, control_probability=control)
        out.append((e, m))
    return out, names


def heads(targets, genes, mean, bulk, variant, exclude=None, stats=None, fill=None, info=None):
    """v4 heads. `fill` (dict, optional) enables coverage_fill for targets no source measures."""
    targets = np.asarray(targets); genes = np.asarray(genes)
    expanded = not variant.startswith('old_') and variant != 'anchor_restricted'
    names = load_banks(exclude, expanded)
    single = '_single_' in variant
    method = 'iv' if variant.startswith('iv_') else 'atlas'
    shrink = 'ash' if '_ash' in variant else 'js' if '_js' in variant else 'none'
    if variant.startswith('anchor'):
        fused, names = anchor_effects(targets, genes, mean, bulk, single, variant == 'anchor_restricted')
        result = []; cov = []
        for (e, m), (space, control, amp) in zip(fused, [('log2fc', mean, .6), ('bulk_delta', bulk, .3)]):
            sp = 'log2fc' if single else space
            if fill is not None:
                e, m = coverage_fill(e, m, targets, genes, control, sp, info=info, **fill)
            result.append(atlas.desired_mean(control, e, space=sp, amplitude=(.6 if single else amp), clip=3)); cov.append(m)
    else:
        fc, fc_cov, _ = pool(targets, genes, names, 'log2fc', method, shrink, use_cd4=expanded, stats=stats)
        if fill is not None:
            fc, fc_cov = coverage_fill(fc, fc_cov, targets, genes, mean, 'log2fc', info=info, **fill)
        p = desired_cpm(mean, fc) if ('smooth' in variant or method == 'iv') else atlas.desired_mean(mean, fc, space='log2fc', amplitude=.6, clip=3)
        if single:
            b = atlas.desired_mean(bulk, fc, space='log2fc', amplitude=.6, clip=3); bc = fc_cov
        else:
            bd, bc, _ = pool(targets, genes, names, 'bulk_delta', method, shrink, use_cd4=expanded, stats=stats)
            if fill is not None:
                bd, bc = coverage_fill(bd, bc, targets, genes, bulk, 'bulk_delta', info=None, **fill)
            b = atlas.desired_mean(bulk, bd, space='bulk_delta', amplitude=.3, clip=3)
        result = [p, b]; cov = [fc_cov, bc]
    if (variant.startswith('anchor') and variant != 'anchor_nopromoter') or variant.endswith('_guide'):
        result = [promoters(q, c, targets, genes, coverage=cov[i], guide_aware=not variant.startswith('anchor')) for i, (q, c) in enumerate(zip(result, [mean, bulk]))]
    return result, cov


def promoters(prob, control, targets, genes, coverage, guide_aware=True):
    path = CFG['promoters']
    if not Path(path).exists():
        return prob
    if not guide_aware:
        return atlas.apply_promoter_prior(prob, control, targets, genes, path, .15)[0]
    pairs = pd.read_csv(path); out = prob.copy(); ti = {t: i for i, t in enumerate(targets)}; gi = {g: i for i, g in enumerate(genes)}
    for row in pairs.itertuples():
        if row.target not in ti or row.neighbor not in gi:
            continue
        i, j = ti[row.target], gi[row.neighbor]
        if coverage[i, j] >= 2:
            continue
        distance = float(row.guide_distance) if row.guide_supported else float(row.distance)
        if not row.guide_supported and not row.divergent:
            continue
        ramp = np.clip(np.log(max(distance, 500) / 500) / np.log(10), 0, 1); remain = .15 + .85 * ramp
        if not row.guide_supported:
            remain = np.sqrt(remain)
        out[i, j] = min(out[i, j], control[j] * remain)
    out /= out.sum(1, keepdims=True); return out


# ---------------------------------------------------------------------------------------------
# Coverage fill: targets with no direct source measurement
# ---------------------------------------------------------------------------------------------
def control_embedding(ctrl_X, genes, k=32, max_cells=4000, seed=0):
    """Label-free gene embedding from the destination context's OWN control cells: top-k right
    singular vectors of standardised log1p(CP10K) expression (genes x k), unit-normalised rows.
    Genes expressed in < 1% of cells get a zero row (no neighbours)."""
    X = sparse.csr_matrix(ctrl_X, dtype=np.float64)
    rng = np.random.default_rng(seed)
    if X.shape[0] > max_cells:
        X = X[np.sort(rng.choice(X.shape[0], max_cells, replace=False))]
    d = np.asarray(X.sum(1)).ravel(); Y = sparse.diags(1e4 / d) @ X; Y.data = np.log1p(Y.data)
    detect = np.asarray((Y > 0).mean(0)).ravel()
    Y = Y.toarray().astype(np.float32); mu = Y.mean(0); sd = Y.std(0) + 1e-6; Y = (Y - mu) / sd
    from sklearn.utils.extmath import randomized_svd
    _, s, vt = randomized_svd(Y, n_components=k, n_iter=4, random_state=seed); del Y
    E = (vt.T * s[None, :]).astype(np.float64)
    E[detect < .01] = 0
    n = np.linalg.norm(E, axis=1, keepdims=True)
    return np.divide(E, n, out=np.zeros_like(E), where=n > 0)


def coverage_fill(effect, cov, targets, genes, control, space, embedding=None, k_neighbors=10, shrink=.5,
                  min_sim=.2, mean_amp=0.0, knockdown=.13, info=None):
    """Fallback chain for targets whose fused coverage is zero on every gene (no source row):
      1. cross-screen evidence: already used by the fusion (any source family measuring the target);
      2. gene-embedding neighbours: similarity-weighted mean of the fused effects of the k most similar
         COVERED panel targets in the context's control co-expression embedding (cos >= min_sim),
         times `shrink`, own-gene columns of neighbours excluded;
      3. mean response: `mean_amp` x the across-covered-target mean effect (default 0: the v3b
         lesson, a shared average shift on targets wrecked fidelity; kept as an option);
    plus, for every filled target, CRISPRi self-knockdown of the target gene to `knockdown` x control
    (VCC panel median fold_expr 0.13 in K562 GWPS). Covered targets are unchanged."""
    effect = np.array(effect, dtype=np.float64, copy=True); cov = np.array(cov, copy=True)
    covered = (cov > 0).any(1)
    miss = np.flatnonzero(~covered)
    gi = {g: i for i, g in enumerate(genes)}
    rec = dict(space=space, targets=len(targets), covered=int(covered.sum()), filled_neighbors=0, filled_mean=0, knockdown_only=0)
    if len(miss) == 0:
        if info is not None:
            info.setdefault('coverage_fill', []).append(rec)
        return effect.astype(np.float32), cov
    cidx = np.flatnonzero(covered)
    own = np.array([gi.get(t, -1) for t in targets])
    mean_eff = np.zeros(len(genes))
    if cidx.size:
        E = effect[cidx].copy()
        for r, i in enumerate(cidx):
            if own[i] >= 0:
                E[r, own[i]] = np.nan
        mean_eff = np.nan_to_num(np.nanmean(E, 0))
    kd = np.log2(knockdown) if space == 'log2fc' else None
    for i in miss:
        filled = False
        if embedding is not None and own[i] >= 0 and cidx.size and np.any(embedding[own[i]]):
            cand = cidx[own[cidx] >= 0]
            sims = embedding[own[cand]] @ embedding[own[i]]
            ok = sims >= min_sim
            if ok.any():
                order = np.argsort(-sims[ok])[:k_neighbors]; nb = cand[ok][order]; w = sims[ok][order]
                M = effect[nb].copy()
                for r, j in enumerate(nb):
                    if own[j] >= 0:
                        M[r, own[j]] = 0.0
                effect[i] = shrink * (w[:, None] * M).sum(0) / w.sum()
                rec['filled_neighbors'] += 1; filled = True
        if not filled and mean_amp > 0:
            effect[i] = mean_amp * mean_eff; rec['filled_mean'] += 1; filled = True
        if not filled:
            effect[i] = 0.0; rec['knockdown_only'] += 1
        if own[i] >= 0:
            if space == 'log2fc':
                # amplitude .6 is applied downstream by desired_mean; pre-divide so the target lands at `knockdown`
                effect[i, own[i]] = kd / .6
            else:
                c = 50000 * control[own[i]]
                effect[i, own[i]] = (np.log1p(knockdown * c) - np.log1p(c)) / .6
        cov[i] = 1e-3  # marks synthetic coverage (promoter prior may still act: < 2)
    if info is not None:
        info.setdefault('coverage_fill', []).append(rec)
    return effect.astype(np.float32), cov


def template(a, genes, cells=400, pool_k=4, seed=20260910):
    assert np.array_equal(np.asarray(a.var_names), np.asarray(genes))
    x = sparse.csr_matrix(a.X, dtype=float); depth = np.asarray(x.sum(1)).ravel(); assert np.all(depth > 0)
    p = np.asarray(x.multiply((1 / depth)[:, None]).mean(0)).ravel(); b = np.asarray(x.sum(0)).ravel(); b /= b.sum()
    assert len(depth) >= cells * pool_k, f'need >= {cells * pool_k} control cells, have {len(depth)}'
    selected = np.random.default_rng(seed).choice(len(depth), cells * pool_k, replace=False); selected = selected[np.argsort(depth[selected], kind='stable')]
    q = x[selected].toarray() / depth[selected, None]; q = q.reshape(cells, pool_k, len(genes)).mean(1)
    d = np.rint(depth[selected].reshape(cells, pool_k).mean(1)).astype(np.int64)
    return q, d, p, b


def mean_counts(q, p, b, depths, seed):
    x = q.astype(float).copy(); x = .999 * x + .001 * b; desired = depths.sum() * b
    for it in range(100):
        col = depths @ x; x *= np.divide(desired, col, out=np.zeros_like(col), where=col > 0); x /= np.maximum(x.sum(1, keepdims=True), 1e-30)
        if it % 5 == 4 and np.abs(depths @ x / depths.sum() - b).sum() < 2e-4:
            break
    exp = x * depths[:, None]; integer = np.floor(exp).astype(np.int32); rng = np.random.default_rng(seed)
    for i in range(len(x)):
        rem = int(depths[i] - integer[i].sum())
        if rem:
            cum = np.cumsum(exp[i] - integer[i]); cum *= rem / cum[-1]; loc = np.searchsorted(cum, np.arange(rem) + rng.random(), side='right'); np.add.at(integer[i], loc, 1)
    return atlas.repair_bulk_totals(integer, exp)


def emit(q, p, b, depths, seed, variant):
    if variant.endswith('_mean'):
        return mean_counts(q, p, b, depths, seed)
    return atlas.dual_moment_counts(q, p, b, depths=depths, seed=seed)

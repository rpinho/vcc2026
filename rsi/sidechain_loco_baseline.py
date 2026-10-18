"""Fixed available-moment Sidechain kernel baseline; NOT submitted SER17.

MIT release kernels are unmodified. The source adapter uses observed float32
matched-control mean/variance summaries, not reconstructed PseudobulkSums.
The frozen spec records the omitted Poisson floor, neighbours and scatter table.
This Phase1 reproduction never selects a recipe or opens HOLDOUT.
"""
from pathlib import Path
import gc
import hashlib
import sys
import time

import anndata as ad
import numpy as np
from scipy import sparse
import runtime as rt

RELEASE = Path('vendor/sidechain_release')
sys.path.insert(0, str((RELEASE/'src').resolve()))
from sidechain.models.count_emitters import ContextProfile, PoissonEmitter
from sidechain.submit.build import pooled_delta, TARGET_SELF_LOG2FC

SPEC = Path('artifacts/sidechain_loco_spec.json')
AUDIT = Path('artifacts/sidechain_loco_audit.json')
BASELINE_ID = 'sidechain_available_moments_loco_r25'
TASKS = ['H1_rich_search', 'CD4_depth192_search']
FIELDS = ['target_mean_cpm', 'matched_control_mean_cpm',
          'target_mean_cpm_variance', 'matched_control_mean_cpm_variance']


class MomentSource:
    """Public effect-source interface, with observed matched-arm moments only."""
    def __init__(self, path, targets):
        with np.load(path, allow_pickle=False) as data:
            assert bool(data['arithmetic_cpm_observed']), path
            measured = data['measured_genes'].astype(bool)
            self.genes = data['genes'][measured]
            rows = np.flatnonzero(np.isin(data['targets'], targets) & (data['n_cells'] >= 20))
            self.targets = data['targets'][rows]
            self.pos = {t: i for i, t in enumerate(self.targets)}
            self.values = {}
            for key in FIELDS:
                value = data[key][rows][:, measured].astype(np.float64)
                assert value.shape == (len(rows), len(self.genes)), (path, key)
                assert np.isfinite(value).all() and np.all(value >= 0), (path, key)
                self.values[key] = value
        assert len(set(self.genes)) == len(self.genes)
        assert len(set(self.targets)) == len(self.targets)

    def effect(self, target):
        if target not in self.pos:
            return None
        i = self.pos[target]
        p, c, vp, vc = [self.values[key][i] for key in FIELDS]
        fc = np.log2((p + 1.0) / (c + 1.0))
        var = (vp / (p + 1.0)**2 + vc / (c + 1.0)**2) / np.log(2.0)**2
        return fc, var


def read_spec():
    spec = rt.json.loads(SPEC.read_text())
    assert spec['baseline_id'] == BASELINE_ID
    for group in ['implementation_sha256', 'protocol_sha256', 'source_sha256']:
        for path, digest in spec[group].items():
            assert rt.sha(path) == digest, path
    return spec


def control_profile(ctrl, destination):
    # Call the original public method with the input-control object only.
    original = ad.read_h5ad
    ad.read_h5ad = lambda path: ctrl
    try:
        return ContextProfile.from_controls('input-controls-only', name=destination,
                                             min_libsize=1000.0, keep_cells=True)
    finally:
        ad.read_h5ad = original


def shifts(task, spec, genes, stats=None):
    entry = spec['tasks'][task]
    sources = [MomentSource(Path(p), entry['targets']) for p in entry['sources']]
    result = {}
    for target in entry['targets']:
        value = pooled_delta(target, sources, genes, shrinkage=True,
            var_floor='none', gamma=1.0, coverage_tiers=None, similarity_beta=0.0,
            log_bias_correct=False, shrink_k=1.0, shrink_stage='pooled',
            shrink_rule='adaptive', stats=stats)
        assert value is not None, (task, target, 'No covered source: no guessed fallback')
        assert np.isfinite(value).all()
        value *= 1.35
        self_pos = np.flatnonzero(genes == target)
        if len(self_pos):
            value[self_pos[0]] = TARGET_SELF_LOG2FC
        result[target] = value
    del sources
    gc.collect()
    return result


def seed_for(task, target):
    # Explicit local convention, not an assertion about the submitted RNG.
    raw = f'2026101725:{task}:{target}'.encode()
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], 'little')


def predictions(ctrl, task, spec, diagnostics=None):
    entry = spec['tasks'][task]
    genes = ctrl.var_names.to_numpy(str)
    profile = control_profile(ctrl, entry['destination'])
    pool_stats = {}
    vectors = shifts(task, spec, genes, pool_stats)
    records = []
    for target in entry['targets']:
        emitter = PoissonEmitter(profile, seed=seed_for(task, target), lam=0.5,
                                  bulk_anchor='pooled', libsize_quantiles=(0.0, 1.0))
        vector = vectors.pop(target)
        block = emitter.emit_dual(400, vector, vector, shape=True, scatter=None,
            iterations=100, tolerance=2e-4, max_projection=0.03, on_fail='raise')
        block = sparse.csr_matrix(block)
        assert emitter.last_dual == 'dual' and emitter.last_shaped is True
        depth = np.asarray(block.sum(axis=1)).ravel()
        assert block.shape == (400, len(genes))
        assert np.isfinite(block.data).all() and np.all(block.data >= 0)
        assert np.equal(block.data, np.floor(block.data)).all()
        assert np.all((depth > 0) & (depth <= 1_000_000))
        records.append({'target': target, 'dual': emitter.last_dual,
                        'shaped': bool(emitter.last_shaped),
                        'depth_min': int(depth.min()), 'depth_max': int(depth.max())})
        yield target, block
    if diagnostics is not None:
        diagnostics.update({'pooling': pool_stats, 'emissions': records,
                            'retained_input_controls': profile.n_cells,
                            'source_moment_floor': 'none; observed variance only'})


def score(task):
    import benchmark
    import bounded_de
    import bounded_score
    import official_scoring as scoring
    import stat_cache
    spec = read_spec()
    assert task in TASKS and spec['tasks'][task]['split'] == 'search'
    audit = rt.json.loads(AUDIT.read_text())
    assert audit['status'] == 'PASS' and audit['spec_sha256'] == rt.sha(SPEC)
    assert audit['verifier_sha256'] == rt.sha('verify_sidechain_loco.py')
    output = Path('artifacts/baseline_scores')/BASELINE_ID/task
    output.mkdir(parents=True, exist_ok=True)
    path = output/'s2026101725.json'
    if path.exists():
        prior = rt.json.loads(path.read_text())
        assert prior['spec_sha256'] == rt.sha(SPEC) and prior['audit_sha256'] == rt.sha(AUDIT)
        assert prior['anchors_sha256'] == rt.sha(f'artifacts/scores/{task}/anchors.json')
        return prior
    start = time.monotonic()
    contract, ctrl, real = benchmark.setup(task)
    assert contract['split'] == 'search' and contract['context'] == spec['tasks'][task]['destination']
    anchors_path = Path('artifacts/scores')/task/'anchors.json'
    anchors = rt.json.loads(anchors_path.read_text())
    assert anchors['real_sha256'] == scoring.matrix_id(real)
    assert anchors['n_splits'] == 5
    assert anchors['driver_sha256'] == rt.sha('bounded_score.py')
    assert anchors['scorer_sha256'] == rt.sha('official_scoring.py')
    assert scoring.provenance() == spec['official_environment']
    bounded_de.install()
    stat_cache.install()
    diagnostics = {}
    stats = {}
    cfg = scoring.config(cache_real=str(Path('artifacts/scores')/task/'real_cache'))
    genes = ctrl.var_names.to_numpy(str)

    def iterator():
        for target, block in predictions(ctrl, task, spec, diagnostics):
            stats[target] = bounded_score.stat(bounded_score.make(block, target, genes), cfg)
            yield target, block

    tidy, digest = bounded_score.score_iter(real, iterator(), cfg)
    raw = benchmark.aggregate(tidy)
    scaled, overall = scoring.scale(raw, anchors['baseline'], anchors['replicate'])
    # Prediction-only diagnostic: neither a fitted amplitude nor a score target.
    from cell_eval2.metrics.delta import _across_pert_budget
    pb, pm = bounded_score.assemble(stats)
    labels, means = pb['bulk_lognorm']
    exclusions = {i: int(np.flatnonzero(genes == t)[0])
                  for i, t in enumerate(labels) if np.any(genes == t)}
    weights = np.array([1.0/(len(genes)-(i in exclusions)) for i in range(len(labels))])
    budget = float(_across_pert_budget(means, exclusions, weights))
    claim = float(np.dot(weights, pm['bulk_lognorm'].jk))
    diagnostics['prediction_only_B_selfclaim'] = {
        'B': budget, 'full_prediction_jackknife_claim': claim,
        'ratio': budget/claim if claim > 0 else None,
        'claim_scope': 'Upper bound on sum(w*min(jk_pred,k*jk_real)); no reference used'}
    tidy.write_csv(output/'tidy.csv')
    result = {'baseline_id': BASELINE_ID, 'task': task, 'seed': 2026101725,
        'raw': raw, 'scaled': scaled, 'overall': overall, 'wall_s': time.monotonic()-start,
        'prediction_sha256': digest, 'spec_sha256': rt.sha(SPEC), 'audit_sha256': rt.sha(AUDIT),
        'reference_sha256': anchors['real_sha256'], 'anchors_sha256': rt.sha(anchors_path),
        'diagnostics': diagnostics, 'external_official_label': None,
        'status': 'Measured available-moment Sidechain kernels; NOT SER17 or promotion',
        'utc': __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat()}
    rt.save(path, result)
    print('BASELINE SCORE', task, overall, scaled, flush=True)
    return result

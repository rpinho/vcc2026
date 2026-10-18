"""Frozen Phase1 v5b sign-gated, control-axis-deflated B/claim port, using input controls and archived kernels.

This is an available-input LOCO adaptation, not a hidden-panel reproduction or
a replacement for the historical calibration gate. No reference enters profile
construction, the forty-target claim probe, or amplitude selection.
"""
from contextlib import contextmanager
from pathlib import Path
import gc
import importlib.util
import time

import runtime as rt
import numpy as np
from scipy import sparse
import recipes

SNAPSHOT = Path('artifacts/reproduction_inputs/v5b_round33')
SPEC = Path('artifacts/v5b_loco_spec.json')
AUDIT = Path('artifacts/v5b_loco_audit.json')
BASELINE_ID = 'v5b_available_loco_r33'
TASKS = ['H1_rich_search', 'CD4_depth192_search']
SEED = 20260910
VARIANT = 'anchor_single_mean'


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


helpers = load_module('frozen_v5b_helpers', SNAPSHOT/'v5.py')


def read_spec():
    spec = rt.json.loads(SPEC.read_text())
    assert spec['baseline_id'] == BASELINE_ID
    for group in ['implementation_sha256', 'protocol_sha256', 'source_sha256']:
        for path, digest in spec[group].items():
            assert rt.sha(path) == digest, path
    import official_scoring
    assert official_scoring.provenance() == spec['official_environment']
    return spec


@contextmanager
def local_v4(destination):
    """Exact source/absent-CD4 adaptation used by recipes.predict(v4a)."""
    v4 = recipes.v4
    old = v4.BANKS, v4.load_banks, v4.cd4_family, v4.atlas.aligned_cd4
    try:
        if destination == 'CD4':
            v4.BANKS = rt.ROOT/'artifacts/cd4_banks'
            v4.cd4_family = lambda *a, **k: None
            v4.atlas.aligned_cd4 = lambda path, targets, genes, *a, **k: (
                np.zeros((len(targets), len(genes)), np.float32),
                np.zeros((len(targets), len(genes)), bool))
        v4.load_banks = lambda *a, **k: ['K562_GWPS_raw', 'HCT116']
        yield v4
    finally:
        v4.BANKS, v4.load_banks, v4.cd4_family, v4.atlas.aligned_cd4 = old


def heads(ctrl, task, spec):
    entry = spec['tasks'][task]
    assert task in TASKS and entry['split'] == 'search'
    assert set(ctrl.obs.evaluation_role.astype(str)) == {'input_control'}
    genes = ctrl.var_names.to_numpy(str)
    assert list(genes) == entry['genes']
    with local_v4(entry['destination']) as v4:
        q, depth, mean, bulk = v4.template(ctrl, genes, pool_k=4, seed=SEED)
        (_, b), _ = v4.heads(np.array(entry['targets']), genes, mean, bulk,
                            VARIANT, exclude=entry['destination'], stats={})
    assert np.isfinite(b).all() and np.all(b >= 0)
    return q, depth, b, bulk


def emit(q, depth, b, destination, target):
    seed = recipes.v4.atlas.seed_for(destination+':'+target, SEED)
    original_repair = recipes.v4.atlas.repair_bulk_totals

    def checked_repair(integer, expected):
        fitted = expected.sum(0)/depth.sum()
        error = float(np.abs(fitted-b).sum())
        assert np.isfinite(error) and error < 2e-4, ('Unconverged archived fit', error)
        return original_repair(integer, expected)

    recipes.v4.atlas.repair_bulk_totals = checked_repair
    try:
        block = recipes.v4.emit(q, None, b, depth, seed, VARIANT)
    finally:
        recipes.v4.atlas.repair_bulk_totals = original_repair
    assert block.dtype.kind in 'iu' and np.all(block >= 0)
    assert np.array_equal(block.sum(1), depth)
    assert np.all((depth > 0) & (depth <= 1_000_000))
    return sparse.csr_matrix(block)


def profile(ctrl, task, spec):
    """The original arm-b profile, retaining the actual uncapped gate."""
    entry = spec['tasks'][task]
    targets = np.array(entry['targets'])
    genes = ctrl.var_names.to_numpy(str)
    q, depth, b, bulk = heads(ctrl, task, spec)
    L = helpers.logbulk(b)
    Lc = helpers.logbulk(bulk)[None, :]
    V, singular = helpers.control_axes(ctrl.X, k=50, seed=entry['probe_rng_seed'])
    L = Lc + helpers.deflate(L-Lc, V)
    with local_v4(entry['destination']) as model:
        names = model.load_banks(None, True)
        sources = [model.aligned_bank(n, targets, genes, 'log2fc') for n in names]
        cd = model.cd4_family(targets, genes, 'log2fc')
        if cd is not None:
            sources.append(cd)
        mask, net = helpers.de_call_mask(sources, L-Lc, k=18533)
    L = Lc + np.where(mask, L-Lc, 0.0)
    info = {'B_v4a_profile': float(helpers.budget(b, targets, genes)),
        'calls_k': 18533, 'calls_per_target': mask.sum(1).tolist(),
        'calls_zero_targets': int((mask.sum(1) == 0).sum()),
        'deflation_k': 50, 'control_sv_top5': singular[:5].tolist(),
        'gate_scope': 'Predicted bulk mask; actual Wilcoxon call count is not constrained',
        'emitter_shaping': 'None; original arm-b emit_worker ignores sharp'}
    return q, depth, L, Lc, mask, info


def calibrated_profile(ctrl, task, spec):
    entry = spec['tasks'][task]
    targets = np.array(entry['targets'])
    genes = ctrl.var_names.to_numpy(str)
    q, depth, L, Lc, mask, profile_info = profile(ctrl, task, spec)
    b0 = helpers.from_logbulk(L)
    probe = np.array(entry['probe_indices'], dtype=int)
    jks = []
    for number, i in enumerate(probe):
        block = emit(q, depth, b0[i], entry['destination'], targets[i])
        jk = helpers.jackknife_correction(block, np.zeros(400, np.intp), 1,
                                          helpers.TS, chunk=32)[0]
        jks.append(float(jk))
        if number % 10 == 0:
            print(task, 'fixed claim probe', number+1, '/', len(probe), flush=True)
    weights = helpers._row_weights(len(targets), len(genes), helpers.excl_cols(targets, genes))
    claim = float(np.mean(jks)*weights.sum())
    assert np.isfinite(claim) and claim > 0
    alpha, B = helpers.solve_alpha(L, targets, genes, claim, center=Lc)
    # Do not silently accept the archived helper's saturated hi=8 return.
    assert 0 < alpha < 8 and np.isclose(B, claim, rtol=1e-10, atol=1e-13), (alpha, B, claim)
    bf = helpers.from_logbulk(helpers.rescale(L, alpha, Lc))
    assert np.isfinite(bf).all() and np.all(bf >= 0)
    diagnostics = {'probe_indices': probe.tolist(), 'probe_targets': targets[probe].tolist(),
        'probe_jackknife': jks, 'claim_est_probe': claim, 'alpha': float(alpha),
        'rho': 1.0, 'B_v4a_profile': profile_info['B_v4a_profile'],
        'B_target_profile': float(B), 'centre': 'input-control pooled bulk only',
        'probe_scope': 'Own prediction full jackknife; before reference-side min cap'}
    diagnostics.update(profile_info)
    print(task, 'frozen amplitude', alpha, 'B', B, 'claim_est', claim, flush=True)
    return q, depth, bf, diagnostics


def predictions(ctrl, task, spec, diagnostics=None):
    q, depth, bf, info = calibrated_profile(ctrl, task, spec)
    entry = spec['tasks'][task]
    for i, target in enumerate(entry['targets']):
        yield target, emit(q, depth, bf[i], entry['destination'], target)
    if diagnostics is not None:
        diagnostics.update(info)


def score(task):
    import benchmark
    import bounded_de
    import bounded_score
    import official_scoring as scoring
    import stat_cache
    spec = read_spec()
    assert task in TASKS
    audit = rt.json.loads(AUDIT.read_text())
    assert audit['status'] == 'PASS' and audit['spec_sha256'] == rt.sha(SPEC)
    assert audit['verifier_sha256'] == rt.sha('verify_v5b_loco.py')
    assert set(audit['tasks']) == set(TASKS)
    output = Path('artifacts/baseline_scores')/BASELINE_ID/task
    output.mkdir(parents=True, exist_ok=True)
    path = output/'s20260910.json'
    anchors_path = Path('artifacts/scores')/task/'anchors.json'
    if path.exists():
        record = rt.json.loads(path.read_text())
        assert record['spec_sha256'] == rt.sha(SPEC) and record['audit_sha256'] == rt.sha(AUDIT)
        assert record['anchors_sha256'] == rt.sha(anchors_path)
        return record
    start = time.monotonic()
    contract, ctrl, real = benchmark.setup(task)
    assert contract['split'] == 'search' and contract['context'] == spec['tasks'][task]['destination']
    anchors = rt.json.loads(anchors_path.read_text())
    assert anchors['real_sha256'] == scoring.matrix_id(real) and anchors['n_splits'] == 5
    assert anchors['driver_sha256'] == rt.sha('bounded_score.py')
    assert anchors['scorer_sha256'] == rt.sha('official_scoring.py')
    bounded_de.install()
    stat_cache.install()
    cfg = scoring.config(cache_real=str(Path('artifacts/scores')/task/'real_cache'))
    diagnostics, stats = {}, {}
    genes = ctrl.var_names.to_numpy(str)

    def iterator():
        for target, block in predictions(ctrl, task, spec, diagnostics):
            stats[target] = bounded_score.stat(bounded_score.make(block, target, genes), cfg)
            yield target, block

    tidy, digest = bounded_score.score_iter(real, iterator(), cfg, release=benchmark.release_models)
    raw = benchmark.aggregate(tidy)
    scaled, overall = scoring.scale(raw, anchors['baseline'], anchors['replicate'])
    pb, pm = bounded_score.assemble(stats)
    labels, means = pb['bulk_lognorm']
    w = helpers._row_weights(len(labels), len(genes), helpers.excl_cols(labels, genes))
    B = float(helpers.budget(means, labels, genes, is_log=True))
    claim = float(np.dot(w, pm['bulk_lognorm'].jk))
    diagnostics['prediction_only_B_selfclaim'] = {'B': B,
        'full_prediction_jackknife_claim': claim, 'ratio': B/claim,
        'claim_scope': 'Upper bound on sum(w*min(jk_pred,k*jk_real)); no reference used'}
    tidy.write_csv(output/'tidy.csv')
    record = {'baseline_id': BASELINE_ID, 'task': task, 'seed': SEED,
        'raw': raw, 'scaled': scaled, 'overall': overall, 'wall_s': time.monotonic()-start,
        'prediction_sha256': digest, 'spec_sha256': rt.sha(SPEC), 'audit_sha256': rt.sha(AUDIT),
        'reference_sha256': anchors['real_sha256'], 'anchors_sha256': rt.sha(anchors_path),
        'diagnostics': diagnostics, 'external_official_label': None,
        'status': 'Measured fixed v5b available-input LOCO adaptation; not a promotion',
        'utc': __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat()}
    rt.save(path, record)
    read_spec()
    print('BASELINE SCORE', task, overall, scaled, flush=True)
    gc.collect()
    return record

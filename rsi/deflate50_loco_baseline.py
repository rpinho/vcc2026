"""Fixed Phase1 description-based control-axis port on audited AtlasShift.

Primary cached Becker/Hao descriptions specify bulk-channel-only projection.
Not an exact reproduction of Becker c6d81cc or its submitted input banks.
"""
from pathlib import Path
import gc
import importlib.util
import numpy as np
from scipy import sparse
import runtime as rt
import atlas_loco_baseline as atlas

SPEC = Path('artifacts/deflate50_loco_spec.json')
AUDIT = Path('artifacts/deflate50_loco_audit.json')
BASELINE_ID = 'atlasshift_deflate50_available_loco_r45'
TASKS = atlas.TASKS


def read_spec():
    spec = rt.json.loads(SPEC.read_text())
    assert spec['baseline_id'] == BASELINE_ID
    atlas.read_spec()
    for group in ['implementation_sha256', 'protocol_sha256']:
        for path, digest in spec[group].items():
            assert rt.sha(path) == digest, path
    import official_scoring
    assert official_scoring.provenance() == spec['official_environment']
    return spec


def input_controls(task):
    import anndata as ad
    assert task in TASKS
    a = ad.read_h5ad(f'artifacts/tasks/{task}.h5ad', backed='r')
    try:
        roles = a.obs.evaluation_role.astype(str).to_numpy()
        rows = np.flatnonzero(roles == 'input_control')
        assert set(a.obs.source_row.iloc[rows]).isdisjoint(set(a.obs.source_row.iloc[np.flatnonzero(roles != 'input_control')]))
        return a[rows].to_memory()
    finally:
        a.file.close()


def control_axes(counts, seed, k=50):
    """Same fixed bounded float32 randomized PCA convention as archived v5b."""
    from sklearn.utils.extmath import randomized_svd
    x = sparse.csr_matrix(counts, dtype=np.float64)
    if x.shape[0] > 4000:
        chosen = np.sort(np.random.default_rng(seed).choice(x.shape[0], 4000, replace=False))
        x = x[chosen]
    depth = np.asarray(x.sum(axis=1)).ravel()
    assert np.all(depth > 0) and min(x.shape) > k
    y = sparse.diags(50000/depth) @ x
    y.data = np.log1p(y.data)
    y = y.toarray().astype(np.float32)
    y -= y.mean(axis=0, keepdims=True)
    _, s, vt = randomized_svd(y, n_components=k, n_iter=5, random_state=seed)
    return vt.T.astype(np.float64), s


def project(effect, axes):
    return effect - (effect @ axes) @ axes.T


def heads(ctrl, task, spec, diagnostics=None):
    assert task in TASKS and spec['tasks'][task]['split'] == 'search'
    assert set(ctrl.obs.evaluation_role.astype(str)) == {'input_control'}
    assert set(ctrl.obs.target_gene.astype(str)) == {'non-targeting'}
    axes, singular = control_axes(ctrl.X, spec['axis_seed'])
    model, _ = atlas.release_modules()
    original = model.desired_mean
    captures = {}
    # This narrow hook leaves every released source, mask, CD4 vote, mean
    # channel, promoter and emission operation unchanged. Only bulk effect
    # AFTER family fusion / BEFORE desired_mean amplitude+clip is projected.
    def desired(control, effect, *, space, amplitude, clip):
        if space == 'bulk_delta':
            before = float(np.square(effect).sum())
            effect = project(effect, axes)
            captures.update({'bulk_effect_squared_norm_before': before,
                'bulk_effect_squared_norm_after': float(np.square(effect).sum())})
        return original(control, effect, space=space, amplitude=amplitude, clip=clip)
    model.desired_mean = desired
    try:
        result = atlas.heads(ctrl, task, spec)
    finally:
        model.desired_mean = original
    if diagnostics is not None:
        diagnostics.update(captures)
        diagnostics.update({'axis_count': 50, 'eta': 1, 'axis_seed': spec['axis_seed'],
            'control_singular_top5': singular[:5].tolist(), 'axis_scope': 'input controls only',
            'projection_location': 'bulk_delta after add_cd4_family before desired_mean',
            'mean_CPM_channel_unchanged': True, 'B_claim_or_call_gating': None})
    return result


def predictions(ctrl, task, spec, diagnostics=None):
    model, _ = atlas.release_modules()
    q, depths, desired, _ = heads(ctrl, task, spec, diagnostics)
    for i, target in enumerate(spec['tasks'][task]['targets']):
        seed = model.seed_for(spec['tasks'][task]['destination']+':'+target, spec['seed'])
        block = model.dual_moment_counts(q, desired[0][i], desired[1][i], depths=depths, seed=seed)
        assert np.array_equal(block.sum(axis=1), depths)
        yield target, sparse.csr_matrix(block)


def score(task):
    import time
    import benchmark
    import bounded_de
    import bounded_score
    import official_scoring as scoring
    import stat_cache
    spec = read_spec()
    audit = rt.json.loads(AUDIT.read_text())
    assert audit['status'] == 'PASS' and audit['spec_sha256'] == rt.sha(SPEC)
    assert audit['verifier_sha256'] == rt.sha('verify_deflate50_loco.py')
    assert set(audit['tasks']) == set(TASKS)
    output = Path('artifacts/baseline_scores')/BASELINE_ID/task
    output.mkdir(parents=True, exist_ok=True)
    path = output/'s20260910.json'
    if path.exists():
        record = rt.json.loads(path.read_text())
        assert record['spec_sha256'] == rt.sha(SPEC) and record['audit_sha256'] == rt.sha(AUDIT)
        return record
    start = time.monotonic()
    contract, ctrl, real = benchmark.setup(task)
    anchors_path = Path(f'artifacts/scores/{task}/anchors.json')
    anchors = rt.json.loads(anchors_path.read_text())
    assert contract['split'] == 'search' and anchors['n_splits'] == 5
    assert anchors['real_sha256'] == scoring.matrix_id(real)
    assert anchors['driver_sha256'] == rt.sha('bounded_score.py')
    assert anchors['scorer_sha256'] == rt.sha('official_scoring.py')
    scoring.scale(anchors['baseline'], anchors['baseline'], anchors['replicate'])
    bounded_de.install(); stat_cache.install()
    details = {}
    tidy, digest = bounded_score.score_iter(real, predictions(ctrl, task, spec, details),
        scoring.config(cache_real=f'artifacts/scores/{task}/real_cache'), release=benchmark.release_models)
    raw = benchmark.aggregate(tidy)
    scaled, overall = scoring.scale(raw, anchors['baseline'], anchors['replicate'])
    tidy.write_csv(output/'tidy.csv')
    record = {'baseline_id': BASELINE_ID, 'task': task, 'seed': spec['seed'],
        'raw': raw, 'scaled': scaled, 'overall': overall, 'wall_s': time.monotonic()-start,
        'prediction_sha256': digest, 'spec_sha256': rt.sha(SPEC), 'audit_sha256': rt.sha(AUDIT),
        'reference_sha256': anchors['real_sha256'], 'anchors_sha256': rt.sha(anchors_path),
        'diagnostics': details, 'external_official_label': None,
        'utc': __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat()}
    rt.save(path, record)
    read_spec()
    print('BASELINE SCORE', task, overall, scaled, flush=True)
    gc.collect()
    return record

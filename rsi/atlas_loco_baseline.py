"""Scope-limited AtlasShift LOCO baseline, using the unmodified MIT release.

This is Phase 1 reproduction, not SER17, a new variant, or a reproduction of
AtlasShift's official score. The frozen spec names every omitted source and
available-panel difference. Historical recipe/provenance files stay untouched.
"""
import runtime as rt
import gc
import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

SPEC = Path('artifacts/atlas_loco_spec.json')
RELEASE = Path('vendor/atlas_release')
TASKS = ['H1_rich_search', 'CD4_depth192_search']
BASELINE_ID = 'atlasshift_available_loco_r22'
REQUIRED = ['source', 'targets', 'genes', 'n_cells', 'target_count_sums',
            'target_mean_cpm', 'matched_control_probability',
            'matched_control_mean_cpm', 'global_control_probability',
            'global_control_mean_cpm', 'measured_genes']


def release_modules():
    """Load local archived public code without editing it or keeping aliases."""
    name = '_rsi_atlas_release_model'
    if name in sys.modules:
        return sys.modules[name], sys.modules['_rsi_atlas_release_predict']
    spec = importlib.util.spec_from_file_location(name, RELEASE/'model.py')
    model = importlib.util.module_from_spec(spec)
    sys.modules[name] = model
    spec.loader.exec_module(model)
    previous = sys.modules.get('model')
    try:
        sys.modules['model'] = model
        spec = importlib.util.spec_from_file_location('_rsi_atlas_release_predict', RELEASE/'predict.py')
        predict = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = predict
        spec.loader.exec_module(predict)
    finally:
        if previous is None:
            sys.modules.pop('model', None)
        else:
            sys.modules['model'] = previous
    return model, predict


def validate_bank(path):
    """Reject missing moments; float32 is the *published AtlasShift* schema."""
    with np.load(path, allow_pickle=False) as data:
        absent = set(REQUIRED)-set(data.files)
        if absent:
            raise ValueError(f'{path}: missing {sorted(absent)}')
        targets, genes = data['targets'], data['genes']
        assert len(set(targets)) == len(targets) and len(set(genes)) == len(genes)
        assert data['n_cells'].shape == (len(targets),)
        assert np.all(data['n_cells'] >= 0)
        assert data['measured_genes'].shape == (len(genes),)
        schema = {}
        for key in REQUIRED:
            value = data[key]
            schema[key] = {'shape': list(value.shape), 'dtype': str(value.dtype)}
            if value.dtype.kind in 'fiu':
                assert np.isfinite(value).all() and np.all(value >= 0), (path, key)
        for key in ['target_count_sums', 'target_mean_cpm', 'matched_control_probability',
                    'matched_control_mean_cpm']:
            assert data[key].shape == (len(targets), len(genes)), (path, key)
        measured = data['measured_genes'].astype(bool)
        assert measured.any()
        return {'sha256': rt.sha(path), 'bytes': path.stat().st_size,
                'schema': schema, 'targets': len(targets), 'genes': len(genes),
                'measured_genes': int(measured.sum()),
                'target_order_sha256': __import__('hashlib').sha256('\n'.join(targets).encode()).hexdigest(),
                'exact_complete_published_source_preparation': False}


def read_spec(check_files=True):
    spec = rt.json.loads(SPEC.read_text())
    assert spec['baseline_id'] == BASELINE_ID
    if check_files:
        for path, digest in spec['implementation_sha256'].items():
            assert rt.sha(path) == digest, path
        for path, digest in spec['protocol_sha256'].items():
            assert rt.sha(path) == digest, path
        for path, record in spec['source_files'].items():
            assert rt.sha(path) == record['sha256'], path
    return spec


def heads(ctrl, task, spec):
    """Only the input control arm and allowed public source summaries enter."""
    entry = spec['tasks'][task]
    targets = np.asarray(entry['targets'])
    genes = ctrl.var_names.to_numpy(str)
    model, published = release_modules()
    # The public helper accepts a file; retain its exact source and operation
    # order by providing this control-only object through its module binding.
    original = published.ad.read_h5ad
    published.ad.read_h5ad = lambda path: ctrl
    try:
        q, depths, mean, bulk = published.control_template(
            'control-only-memory-input', genes, 400, 4, spec['seed'])
    finally:
        published.ad.read_h5ad = original
    sources = [model.load_xatlas(Path(p), prior_counts=100000)
               for p in entry['sources']]
    effects = {space: model.fuse_source_centered(
        sources, entry['weights'], targets, genes, space=space, common_subtract=1)
        for space in ['log2fc', 'bulk_delta']}
    del sources
    if entry['cd4_de'] is not None:
        cd4, available = model.aligned_cd4(Path(entry['cd4_de']), targets, genes,
                                          common_subtract=1, center_scope='source')
    desired, coverage = [], []
    for space, control, amplitude in [('log2fc', mean, .6), ('bulk_delta', bulk, .3)]:
        effect, covered = effects[space]
        if entry['cd4_de'] is not None:
            effect, covered = model.add_cd4_family(effect, covered, cd4, available, .5,
                                                  space=space, control_probability=control)
        probability = model.desired_mean(control, effect, space=space,
                                          amplitude=amplitude, clip=3)
        probability, _ = model.apply_promoter_prior(
            probability, control, targets, genes, Path('artifacts/promoter_pairs.csv'), .15)
        desired.append(probability)
        coverage.append(covered)
    del effects
    gc.collect()
    return q, depths, desired, coverage


def predictions(ctrl, task, spec, diagnostics=None):
    model, _ = release_modules()
    q, depths, desired, coverage = heads(ctrl, task, spec)
    entry = spec['tasks'][task]
    if diagnostics is not None:
        diagnostics.update({
            'targets_with_any_transfer_evidence': int(np.sum(np.any(coverage[0] > 0, axis=1))),
            'readouts_with_transfer_evidence': int(np.sum(coverage[0] > 0)),
            'total_target_readout_pairs': int(coverage[0].size),
            'depth_min': int(depths.min()), 'depth_max': int(depths.max())})
    for index, target in enumerate(entry['targets']):
        seed = model.seed_for(entry['destination']+':'+target, spec['seed'])
        counts = model.dual_moment_counts(q, desired[0][index], desired[1][index],
                                          depths=depths, seed=seed)
        assert np.array_equal(counts.sum(axis=1), depths)
        yield target, sparse.csr_matrix(counts)


def score(task):
    import time
    import benchmark
    import bounded_de
    import bounded_score
    import official_scoring as scoring
    import stat_cache
    spec = read_spec()
    assert task in TASKS and spec['tasks'][task]['split'] == 'search'
    audit = rt.json.loads(Path('artifacts/atlas_loco_audit.json').read_text())
    assert audit['status'] == 'PASS' and audit['spec_sha256'] == rt.sha(SPEC)
    assert audit['verifier_sha256'] == rt.sha('verify_atlas_loco.py')
    output = Path('artifacts/baseline_scores')/BASELINE_ID/task
    output.mkdir(parents=True, exist_ok=True)
    path = output/'s20260910.json'
    if path.exists():
        prior = rt.json.loads(path.read_text())
        assert prior['spec_sha256'] == rt.sha(SPEC)
        return prior
    start = time.monotonic()
    contract, ctrl, real = benchmark.setup(task)
    assert contract['split'] == 'search' and contract['context'] == spec['tasks'][task]['destination']
    anchor_path = Path('artifacts/scores')/task/'anchors.json'
    anchors = rt.json.loads(anchor_path.read_text())
    assert anchors['real_sha256'] == scoring.matrix_id(real)
    assert anchors['driver_sha256'] == rt.sha('bounded_score.py')
    assert anchors['scorer_sha256'] == rt.sha('official_scoring.py')
    assert anchors['n_splits'] == 5
    scoring.scale(anchors['baseline'], anchors['baseline'], anchors['replicate'])
    assert scoring.provenance() == spec['official_environment']
    bounded_de.install()
    stat_cache.install()
    details = {}
    tidy, digest = bounded_score.score_iter(
        real, predictions(ctrl, task, spec, details),
        scoring.config(cache_real=str(Path('artifacts/scores')/task/'real_cache')))
    raw = benchmark.aggregate(tidy)
    scaled, overall = scoring.scale(raw, anchors['baseline'], anchors['replicate'])
    tidy.write_csv(output/'tidy.csv')
    result = {'baseline_id': BASELINE_ID, 'task': task, 'seed': spec['seed'],
              'raw': raw, 'scaled': scaled, 'overall': overall,
              'wall_s': time.monotonic()-start, 'prediction_sha256': digest,
              'spec_sha256': rt.sha(SPEC), 'audit_sha256': rt.sha('artifacts/atlas_loco_audit.json'),
              'reference_sha256': anchors['real_sha256'], 'anchors_sha256': rt.sha(anchor_path),
              'diagnostics': details, 'external_official_label': None,
              'status': 'Measured scope-limited public-context baseline; no external score pairing or promotion',
              'utc': __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat()}
    rt.save(path, result)
    print('BASELINE SCORE', task, overall, scaled, flush=True)
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', choices=TASKS, required=True)
    score(parser.parse_args().task)

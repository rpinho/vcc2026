"""Fixed available-input port of Gerami's MIT norm-restored transfer.

Released model/data/emitter kernels stay unmodified. Four existing source lines,
symbol-native axes and frozen source summaries replace the full eight-source
preparation. This Phase1 baseline has no paired external official label.
"""
import runtime as rt
import gc
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

RELEASE = Path('vendor/zsp_release')
SPEC = Path('artifacts/zsp_loco_spec.json')
AUDIT = Path('artifacts/zsp_loco_audit.json')
BASELINE_ID = 'zsp_equal_norm_available_loco_r41'
TASKS = ['H1_rich_search', 'CD4_depth192_search']
SOURCES = {
    'k562': 'artifacts/sources/K562_GWPS_bulk.npz',
    'hepg2': 'artifacts/sources/HepG2.npz',
    'jurkat': 'artifacts/sources/Jurkat.npz',
    'hct116': '../public_data/model_v4_signature_banks/HCT116.npz',
}


def kernels():
    src = str((RELEASE/'src').resolve())
    if src not in sys.path:
        sys.path.insert(0, src)
    from zsp import data, models, emit
    assert Path(models.__file__).resolve() == (RELEASE/'src/zsp/models.py').resolve()
    return data, models, emit


def read_spec(check_files=True):
    spec = rt.json.loads(SPEC.read_text())
    assert spec['baseline_id'] == BASELINE_ID
    if check_files:
        for group in ['implementation_sha256', 'protocol_sha256', 'input_sha256']:
            for path, digest in spec[group].items():
                assert rt.sha(path) == digest, path
    return spec


def input_controls(task):
    """Read only SEARCH input-control expression, never reference/response rows."""
    assert task in TASKS
    a = ad.read_h5ad(f'artifacts/tasks/{task}.h5ad', backed='r')
    try:
        rows = np.flatnonzero(a.obs.evaluation_role.astype(str).to_numpy() == 'input_control')
        ctrl = a[rows].to_memory()
    finally:
        a.file.close()
    assert set(ctrl.obs.target_gene.astype(str)) == {'non-targeting'}
    return ctrl


def source_contexts():
    """Use actual pooled raw summaries, not rounded signatures/mean-cell CPM."""
    data, _, _ = kernels()
    contexts = []
    for name, path in SOURCES.items():
        with np.load(path, allow_pickle=False) as bank:
            if name == 'hct116':
                measured = bank['measured_genes'].astype(bool)
                genes = pd.Index(bank['genes'][measured].astype(str))
                targets = bank['targets'].astype(str)
                n = bank['n_cells'].astype(np.float64)
                means = bank['target_count_sums'][:, measured].astype(np.float64) / np.maximum(n[:, None], 1)
                # Only normalized global control probabilities are available.
                # decompose() normalizes this again; its absolute depth cancels.
                control = bank['global_control_probability'][measured].astype(np.float64)
            else:
                genes = pd.Index(bank['genes'].astype(str))
                targets = bank['targets'].astype(str)
                labels = bank['raw_labels'].astype(str)
                assert np.array_equal(labels, np.concatenate((['non-targeting'], targets)))
                means = bank['raw_means'][1:].astype(np.float64)
                control = bank['raw_means'][0].astype(np.float64)
                n = bank['n_cells'].astype(np.float64)
            assert genes.is_unique and len(set(targets)) == len(targets)
            assert np.isfinite(means).all() and (means >= 0).all()
            assert np.isfinite(control).all() and (control >= 0).all() and control.sum() > 0
            assert (means[n >= 10].sum(axis=1) > 0).all()
            contexts.append(data.Context(name, genes, pd.Series(genes, index=genes),
                control, pd.DataFrame(means, index=targets, columns=genes),
                pd.Series(n, index=targets)))
    return data.align_union(contexts)


def heads(ctrl, task, spec, diagnostics=None):
    data, model, _ = kernels()
    contexts = source_contexts()
    sources = [model.decompose(c, min_cells=10) for c in contexts]
    genes = contexts[0].genes
    panel = pd.Index(ctrl.var_names.astype(str))
    positions = genes.get_indexer(panel)
    hit = positions >= 0
    counts = sparse.csr_matrix(ctrl.X)
    # Same sparse float32 mean/order as scripts/05_make_submission.py.
    mean = np.asarray(counts.mean(0)).ravel()
    basal = np.zeros(len(genes))
    basal[positions[hit]] = data.cp10k(mean)[0][hit]
    # Infinite temperature is exactly equal weights. No LOSO/refitting/tuning.
    prediction = model.NormRestoredTransfer(temperature=np.inf)(
        sources, basal, spec['tasks'][task]['targets'], genes, genes)
    lfc = np.zeros((len(prediction.lfc), len(panel)))
    lfc[:, hit] = prediction.lfc.to_numpy()[:, positions[hit]]
    desired = model.lfc_to_counts(lfc, mean)
    assert np.isfinite(desired).all() and (desired >= 0).all()
    if diagnostics is not None:
        target_list = spec['tasks'][task]['targets']
        diagnostics.update({'sources': list(SOURCES), 'readouts_with_any_source': int(hit.sum()),
            'source_union_readouts': len(genes),
            'target_coverage_by_source': {s.name: sum(t in s.lfc.index for t in target_list) for s in sources},
            'lfc_abs_max': float(np.abs(lfc).max()),
            'mean_count_max': float(desired.max())})
    del sources, contexts, prediction
    gc.collect()
    return mean, desired


def predictions(ctrl, task, spec, diagnostics=None):
    _, _, emit = kernels()
    mean, desired = heads(ctrl, task, spec, diagnostics)
    cells = sparse.csr_matrix(ctrl.X)
    # One seed0 RNG for each separately run context, in frozen manifest order.
    rng = np.random.default_rng(spec['seed'])
    for index, target in enumerate(spec['tasks'][task]['targets']):
        x = emit.emit_cells(cells, mean, desired[index], 400, rng)
        assert x.shape == (400, ctrl.n_vars)
        assert np.isfinite(x.data).all() and (x.data >= 0).all()
        assert np.equal(x.data, np.floor(x.data)).all()
        yield target, x


def score(task):
    import time
    import benchmark
    import bounded_de
    import bounded_score
    import official_scoring as scoring
    import stat_cache
    spec = read_spec()
    assert task in TASKS
    audit = rt.json.loads(AUDIT.read_text())
    assert audit['status'] == 'PASS' and audit['spec_sha256'] == rt.sha(SPEC)
    assert audit['verifier_sha256'] == rt.sha('verify_zsp_loco.py')
    output = Path('artifacts/baseline_scores')/BASELINE_ID/task
    output.mkdir(parents=True, exist_ok=True)
    path = output/'s0.json'
    if path.exists():
        prior = rt.json.loads(path.read_text())
        assert prior['spec_sha256'] == rt.sha(SPEC)
        return prior
    start = time.monotonic()
    contract, ctrl, real = benchmark.setup(task)
    assert contract['split'] == 'search'
    anchor_path = Path(f'artifacts/scores/{task}/anchors.json')
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
    tidy, digest = bounded_score.score_iter(real, predictions(ctrl, task, spec, details),
        scoring.config(cache_real=f'artifacts/scores/{task}/real_cache'))
    raw = benchmark.aggregate(tidy)
    scaled, overall = scoring.scale(raw, anchors['baseline'], anchors['replicate'])
    tidy.write_csv(output/'tidy.csv')
    result = {'baseline_id': BASELINE_ID, 'task': task, 'seed': spec['seed'],
        'raw': raw, 'scaled': scaled, 'overall': overall,
        'wall_s': time.monotonic()-start, 'prediction_sha256': digest,
        'spec_sha256': rt.sha(SPEC), 'audit_sha256': rt.sha(AUDIT),
        'reference_sha256': anchors['real_sha256'], 'anchors_sha256': rt.sha(anchor_path),
        'diagnostics': details, 'external_official_label': None,
        'status': 'Measured fixed Phase1 four-source public-method baseline; no promotion'}
    rt.save(path, result)
    print('BASELINE SCORE', task, overall, scaled, flush=True)
    return result

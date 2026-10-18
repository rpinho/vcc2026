"""Independent control-only audit of the exact frozen v5b recipe extension.

Scalar source-sign votes, projection algebra, all40 independent jackknifes,
independent amplitude bisection and unchanged integer emission are checked.
No response matrix, reference-control expression or biological score is read.
"""
from datetime import datetime, timezone
from pathlib import Path
import gc
import hashlib
import math
import anndata as ad
import numpy as np
from scipy import sparse
import runtime as rt
import recipes
import v5b_loco_baseline as baseline
import verify_v5a_loco as independent


def scalar_gate(sources, delta, k):
    """Per-coordinate voting without production vector arithmetic/sorting."""
    n, g = delta.shape
    net = np.zeros((n, g), np.float32)
    mask = np.zeros((n, g), bool)
    for i in range(n):
        candidates = []
        for j in range(g):
            sign = int(delta[i, j] > 0)-int(delta[i, j] < 0)
            vote, strength = 0, 0.0
            for effect, available, variance in sources:
                if not available[i, j] or not math.isfinite(float(variance[i, j])):
                    continue
                z = float(effect[i, j])/math.sqrt(max(float(variance[i, j]), 1e-12))
                direction = int(z > 0)-int(z < 0)
                if direction == sign:
                    strength += abs(z)
                if abs(z) > 1.96 and sign:
                    if direction == sign:
                        vote += 1
                    elif direction == -sign:
                        vote -= 1
            net[i, j] = vote
            if vote >= 1:
                candidates.append((-(vote*1e6+min(strength, 1e5)), j))
        for _, j in sorted(candidates)[:k]:
            mask[i, j] = True
    return mask, net


def synthetic():
    rng = np.random.default_rng(333)
    delta = np.array([[1, -1, 0, 1, 1, 1, 1, 1]], float)
    e = np.array([[2, -2, 8, 1.96, -8, 5, 2, 2]], float)
    m = np.ones_like(e, bool); m[0, 5] = False
    v = np.ones_like(e); v[0, 4] = np.inf
    sources = [(e, m, v)]
    for k in [2, 18533]:
        expected = scalar_gate(sources, delta, k)
        actual = baseline.helpers.de_call_mask(sources, delta, k=k)
        assert all(np.array_equal(a, b) for a, b in zip(actual, expected))
    assert np.flatnonzero(expected[0]).tolist() == [0, 1, 6, 7]
    # Deflation may produce negative log-bulk coordinates. The archived
    # inverse clips expm1 at zero before row normalization.
    log_fixture = np.array([[-0.4, 0.0, 0.7, 1.5]], np.float64)
    inverse = np.maximum(np.expm1(log_fixture), 0.0)
    inverse /= inverse.sum(1, keepdims=True)
    np.testing.assert_array_equal(baseline.helpers.from_logbulk(log_fixture), inverse)
    assert inverse[0, 0] == 0 and inverse[0, 1] == 0
    # Exact PCA projector on a small, separated-spectrum count fixture.
    counts = rng.poisson(np.geomspace(1, 35, 12), size=(90, 12))
    V, _ = baseline.helpers.control_axes(sparse.csr_matrix(counts), k=4, seed=1)
    y = np.log1p(50000*counts/counts.sum(1, keepdims=True)).astype(np.float32)
    y -= y.mean(0, keepdims=True)
    _, _, vt = np.linalg.svd(y, full_matrices=False)
    np.testing.assert_allclose(V@V.T, vt[:4].T@vt[:4], rtol=0, atol=3e-6)
    return {'scalar_gate_boundary_missing_variance_tie': 'PASS', 'small_exact_PCA_projector': 'PASS',
            'negative_logbulk_inverse_clip': 'PASS'}


def audit_task(task, spec):
    entry = spec['tasks'][task]
    path = Path('artifacts/tasks')/(task+'.h5ad')
    assert rt.sha(path) == rt.json.loads(path.with_suffix('.json').read_text())['sha256']
    a = ad.read_h5ad(path, backed='r')
    selected = a.obs.evaluation_role.astype(str).to_numpy() == 'input_control'
    ctrl = a[selected].to_memory()
    assert set(ctrl.obs.source_row).isdisjoint(set(a.obs.source_row[~selected]))
    a.file.close(); del a
    before = independent.fingerprint(ctrl)
    targets, genes = np.array(entry['targets']), ctrl.var_names.to_numpy(str)
    assert set(ctrl.obs.target_gene.astype(str)) == {'non-targeting'}
    assert set(ctrl.obs.evaluation_role.astype(str)) == {'input_control'}
    assert list(genes) == entry['genes']
    probes = np.sort(np.random.default_rng(entry['probe_rng_seed']).choice(79, 40, replace=False))
    assert probes.tolist() == entry['probe_indices']
    with independent.no_h5ad_reads():
        # Separately loaded archived model with explicit LOCO bank roots.
        model, q, depth, b, bulk = independent.original_heads(ctrl, task, spec)
        Lc = np.log1p(50000*bulk)[None, :]
        delta = np.log1p(50000*b)-Lc
        V, singular = baseline.helpers.control_axes(ctrl.X, k=50, seed=entry['probe_rng_seed'])
        np.testing.assert_allclose(V.T@V, np.eye(50), rtol=0, atol=3e-6)
        projected = delta-(delta@V)@V.T
        projection_error_bound = (np.linalg.norm(delta@V)*np.linalg.norm(np.eye(50)-V.T@V)
                                  + 1e-9*np.linalg.norm(delta))
        assert np.linalg.norm(projected@V) <= projection_error_bound
        np.testing.assert_allclose(baseline.helpers.deflate(delta, V), projected, rtol=0, atol=0)
        allowed = ['K562_GWPS_raw', 'HCT116']
        sources = [model.aligned_bank(n, targets, genes, 'log2fc') for n in allowed]
        cd = model.cd4_family(targets, genes, 'log2fc')
        if cd is not None:
            sources.append(cd)
        assert len(sources) == (2 if entry['destination'] == 'CD4' else 3)
        # Preserve the original add/subtract round trip before consensus signs.
        raw_L = Lc+projected
        expected_mask, expected_net = scalar_gate(sources, raw_L-Lc, 18533)
        actual_mask, actual_net = baseline.helpers.de_call_mask(sources, raw_L-Lc, k=18533)
        assert np.array_equal(actual_mask, expected_mask) and np.array_equal(actual_net, expected_net)
        L = Lc+np.where(expected_mask, raw_L-Lc, 0)
        del sources, raw_L, projected, delta
        qq, dd, actual_L, actual_Lc, mask, profile_info = baseline.profile(ctrl, task, spec)
        for actual, expected in [(qq,q), (dd,depth), (actual_L,L), (actual_Lc,Lc), (mask,expected_mask)]:
            assert np.array_equal(actual, expected), (task, 'original arm-b profile')
        assert profile_info['calls_k'] == 18533 and profile_info['deflation_k'] == 50
        np.testing.assert_array_equal(profile_info['control_sv_top5'], singular[:5])
        negative_logbulk_coordinates = int(np.count_nonzero(L < 0))
        b0 = np.maximum(np.expm1(L), 0.0)
        b0 /= b0.sum(1, keepdims=True)
        np.testing.assert_array_equal(b0, baseline.helpers.from_logbulk(L))
        observed = []
        original_emit = baseline.emit
        def observe_emit(q_, d_, profile, destination, target):
            i = probes[len(observed)]
            assert target == targets[i] and destination == entry['destination']
            assert np.array_equal(profile, b0[i])
            block = original_emit(q_, d_, profile, destination, target)
            observed.append(independent.independent_jk(block))
            if len(observed) <= 2:
                expected = model.emit(q, None, profile, depth,
                    model.atlas.seed_for(destination+':'+target, 20260910), 'anchor_single_mean')
                assert np.array_equal(block.toarray(), expected)
            assert block.shape == (400, len(genes)) and block.dtype.kind in 'iu'
            assert np.array_equal(np.asarray(block.sum(1)).ravel(), depth)
            return block
        baseline.emit = observe_emit
        try:
            qq, dd, final, info = baseline.calibrated_profile(ctrl, task, spec)
        finally:
            baseline.emit = original_emit
        assert len(observed) == 40
        np.testing.assert_allclose(observed, info['probe_jackknife'], rtol=2e-11, atol=1e-10)
        weights = independent.independent_weights(targets, genes)
        claim = math.fsum(observed)/40*math.fsum(weights)
        np.testing.assert_allclose(info['claim_est_probe'], claim, rtol=2e-11, atol=1e-12)
        alpha, expected_final = independent.independent_solve(L, Lc, targets, genes, info['claim_est_probe'])
        np.testing.assert_allclose(info['alpha'], alpha, rtol=2e-10, atol=1e-12)
        np.testing.assert_allclose(final, expected_final, rtol=2e-10, atol=1e-14)
        np.testing.assert_allclose(independent.independent_budget(np.log1p(50000*final), targets, genes),
                                   info['B_target_profile'], rtol=2e-10, atol=1e-12)
        assert before == independent.fingerprint(ctrl), 'Control mutation'
    recipes.v4._POOL_CACHE.clear()
    gc.collect()
    return {'input_fingerprint': before, 'input_cells': ctrl.n_obs,
            'reference_expression_read': False, 'source_families': allowed+(['CD4'] if cd is not None else []),
            'scalar_vote_coordinates': len(targets)*len(genes), 'original_profile_bitwise': True,
            'gate_mask_sha256': hashlib.sha256(mask.tobytes()).hexdigest(),
            'negative_logbulk_coordinates': negative_logbulk_coordinates,
            'independent_probe_count': len(observed), 'archived_emission_checks': 2,
            'orthonormal_projection_check': 'PASS', 'diagnostics': info}


def main():
    assert not baseline.AUDIT.exists(), 'Never overwrite an audit receipt'
    spec = baseline.read_spec()
    fixtures = synthetic()
    records = {}
    for task in baseline.TASKS:
        records[task] = audit_task(task, spec)
        print('PASS control-only v5b audit', task, flush=True)
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, e in index['entries'].items():
        assert rt.sha(path) == e['sha256'], path
    assert rt.sha(index['archive']) == index['archive_sha256']
    baseline.read_spec()
    assert rt.json.loads(Path('artifacts/calibration_gate.json').read_text())['status'] == 'CLOSED'
    assert not Path('artifacts/holdout_events.jsonl').exists()
    rt.save(baseline.AUDIT, {'utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS',
        'spec_sha256': rt.sha(baseline.SPEC), 'verifier_sha256': rt.sha(__file__),
        'fixtures': fixtures, 'tasks': records, 'frozen_metadata_members_unchanged': len(index['entries']),
        'scope': 'Fixed v5b available-input LOCO port; controls only; no biological score',
        'external_official_label': None, 'holdout_scored': False, 'phase2_allowed': False})
    print('PASS frozen v5b adapter; eligible for fixed Phase1 SEARCH reproduction only', flush=True)


if __name__ == '__main__':
    main()

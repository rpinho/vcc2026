"""Independent numerical/identity/isolation audit before biological scoring.

Replays the hash-bound original model with explicit source restrictions; checks
full-panel heads, two original local v4a emissions, all forty probe jackknifes,
separate weighted-centering/bisection arithmetic and original integer repair.
The audit reads SEARCH controls only, with model-time H5AD reads forbidden.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import gc
import hashlib
import math

import runtime as rt
import anndata as ad
import numpy as np
from scipy import sparse
import provenance
import recipes
import v5a_loco_baseline as baseline


def fingerprint(ctrl):
    x = sparse.csr_matrix(ctrl.X)
    h = hashlib.sha256()
    for value in [x.data, x.indices, x.indptr]:
        h.update(value.tobytes())
    h.update(ctrl.obs.to_json().encode())
    h.update('\0'.join(ctrl.var_names).encode())
    return h.hexdigest()


def independent_weights(targets, genes):
    axis = set(genes)
    return np.array([1/(len(genes)-(t in axis)) for t in targets])


def independent_budget(L, targets, genes):
    # Exclude each row's target BEFORE per-gene weighted centering. Include
    # mixed resolved/unresolved target weights; no scorer/helper imports here.
    weights = independent_weights(targets, genes)
    W = np.broadcast_to(weights[:, None], L.shape).copy()
    pos = {g: i for i, g in enumerate(genes)}
    for i, target in enumerate(targets):
        if target in pos:
            W[i, pos[target]] = 0
    center = np.sum(W*L, axis=0)/W.sum(0)
    return float(np.sum(W*np.square(L-center)))


def independent_jk(block):
    # Direct delete-one variance, two passes over dense13-row chunks. This
    # neither calls the official routine nor shares its shifted s2-s1^2 form.
    x = sparse.csr_matrix(block, dtype=float)
    counts = np.asarray(x.sum(0)).ravel()
    depths = np.asarray(x.sum(1)).ravel()
    total, n = float(counts.sum()), x.shape[0]
    mean = np.zeros(x.shape[1])
    for lo in range(0, n, 13):
        val = np.log1p(50000*(counts-x[lo:lo+13].toarray())/(total-depths[lo:lo+13, None]))
        mean += val.sum(0)/n
    answer = 0.0
    for lo in range(0, n, 13):
        val = np.log1p(50000*(counts-x[lo:lo+13].toarray())/(total-depths[lo:lo+13, None]))
        answer += float(np.sum(np.square(val-mean)))
    return (n-1)/n*answer


def independent_solve(L, Lc, targets, genes, claim):
    def profile(alpha):
        e = np.maximum(np.expm1(Lc+alpha*(L-Lc)), 0)
        return e/e.sum(1, keepdims=True)
    lo, hi = 0.0, 8.0
    assert independent_budget(np.log1p(50000*profile(hi)), targets, genes) >= claim
    for _ in range(60):
        mid = (lo+hi)/2
        if independent_budget(np.log1p(50000*profile(mid)), targets, genes) < claim:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2, profile((lo+hi)/2)


@contextmanager
def no_h5ad_reads():
    old = ad.read_h5ad
    def forbidden(*args, **kwargs):
        raise AssertionError('Model construction attempted a file read beyond supplied input controls')
    ad.read_h5ad = forbidden
    try:
        yield
    finally:
        ad.read_h5ad = old


def original_heads(ctrl, task, spec):
    """Separate archived module; explicit restriction instead of adapter context."""
    model = baseline.load_module('v5a_independent_archived_model', baseline.SNAPSHOT/'model.py')
    destination = spec['tasks'][task]['destination']
    model.BANKS = rt.ROOT/('artifacts/cd4_banks' if destination == 'CD4' else
                         '../public_data/model_v4_signature_banks')
    model.load_banks = lambda *a, **k: ['K562_GWPS_raw', 'HCT116']
    old_align = model.atlas.aligned_cd4
    try:
        if destination == 'CD4':
            model.cd4_family = lambda *a, **k: None
            model.atlas.aligned_cd4 = lambda path, targets, genes, *a, **k: (
                np.zeros((len(targets), len(genes)), np.float32),
                np.zeros((len(targets), len(genes)), bool))
        genes = ctrl.var_names.to_numpy(str)
        q, d, mean, bulk = model.template(ctrl, genes, seed=20260910)
        heads, _ = model.heads(np.array(spec['tasks'][task]['targets']), genes, mean, bulk,
                              'anchor_single_mean', exclude=destination, stats={})
        return model, q, d, heads[1], bulk
    finally:
        model.atlas.aligned_cd4 = old_align


def audit_task(task, spec):
    print('Auditing frozen v5a', task, flush=True)
    entry = spec['tasks'][task]
    path = Path('artifacts/tasks')/(task+'.h5ad')
    task_receipt = rt.json.loads(path.with_suffix('.json').read_text())
    assert rt.sha(path) == task_receipt['sha256']
    a = ad.read_h5ad(path, backed='r')
    mask = a.obs.evaluation_role.astype(str).to_numpy() == 'input_control'
    ctrl = a[mask].to_memory()
    reference_rows = set(a.obs.source_row[~mask])
    assert reference_rows.isdisjoint(set(ctrl.obs.source_row))
    a.file.close()
    del a
    assert set(ctrl.obs.target_gene.astype(str)) == {'non-targeting'}
    assert set(ctrl.obs.evaluation_role.astype(str)) == {'input_control'}
    before = fingerprint(ctrl)
    genes, targets = ctrl.var_names.to_numpy(str), np.array(entry['targets'])
    assert list(genes) == entry['genes']
    expected_probe = np.sort(np.random.default_rng(entry['probe_rng_seed']).choice(79, 40, replace=False))
    assert expected_probe.tolist() == entry['probe_indices']
    with no_h5ad_reads():
        model, q, d, b, bulk = original_heads(ctrl, task, spec)
        qq, dd, bb, bc = baseline.heads(ctrl, task, spec)
        for actual, expected in [(qq, q), (dd, d), (bb, b), (bc, bulk)]:
            assert np.array_equal(actual, expected), (task, 'archived full-panel bridge')
        del qq, dd, bb, bc
        # Two predetermined emissions reproduce recipes.predict's actual
        # original b (without a log/expm1 round trip).
        base = recipes.predict(ctrl, list(targets), entry['destination'], 'v4a', seed=42)
        base_checks = []
        try:
            for i in range(2):
                target, actual = next(base)
                expected = model.emit(q, None, b[i], d,
                    model.atlas.seed_for(entry['destination']+':'+target, 20260910), 'anchor_single_mean')
                assert target == targets[i] and np.array_equal(actual.toarray(), expected)
                base_checks.append(target)
        finally:
            base.close()
        L = np.log1p(50000*b)
        Lc = np.log1p(50000*bulk)[None, :]
        b0 = np.expm1(L); b0 /= b0.sum(1, keepdims=True)
        observed_jk, repair_checks, probe_checks = [], [], []
        original_emit = baseline.emit
        original_repair = recipes.v4.atlas.repair_bulk_totals

        def audit_repair(integer, expected):
            columns = expected.sum(0)
            totals = np.array([math.floor(float(v)) for v in columns], dtype=np.int64)
            remaining = int(integer.sum()-sum(map(int, totals)))
            order = sorted(range(len(totals)), key=lambda j: (float(columns[j]-totals[j]), j))
            if remaining:
                totals[order[-remaining:]] += 1
            assert 0 <= remaining <= len(totals)
            row_totals = integer.sum(1).copy()
            result = original_repair(integer, expected)
            assert np.array_equal(result.sum(0), totals)
            assert np.array_equal(result.sum(1), row_totals)
            repair_checks.append(True)
            return result

        def observe_emit(q_, depth, profile, destination, target):
            index = len(observed_jk)
            i = expected_probe[index]
            assert target == targets[i] and destination == entry['destination']
            assert np.array_equal(profile, b0[i])
            actual = original_emit(q_, depth, profile, destination, target)
            independent = independent_jk(actual)
            official = float(baseline.helpers.jackknife_correction(actual,
                np.zeros(400, np.intp), 1, 50000, chunk=32)[0])
            np.testing.assert_allclose(independent, official, rtol=2e-11, atol=1e-10)
            if index < 2:
                expected = model.emit(q, None, profile, d,
                    model.atlas.seed_for(destination+':'+target, 20260910), 'anchor_single_mean')
                assert np.array_equal(actual.toarray(), expected)
                probe_checks.append(target)
            observed_jk.append(independent)
            return actual

        baseline.emit = observe_emit
        recipes.v4.atlas.repair_bulk_totals = audit_repair
        try:
            qq, dd, bf, info = baseline.calibrated_profile(ctrl, task, spec)
        finally:
            baseline.emit = original_emit
            recipes.v4.atlas.repair_bulk_totals = original_repair
        assert len(observed_jk) == 40 and len(repair_checks) == 42
        assert np.array_equal(qq, q) and np.array_equal(dd, d)
        w = independent_weights(targets, genes)
        assert np.array_equal(w, baseline.helpers._row_weights(79, len(genes),
            baseline.helpers.excl_cols(targets, genes)))
        claim = math.fsum(observed_jk)/40*math.fsum(w)
        np.testing.assert_allclose(claim, info['claim_est_probe'], rtol=2e-11, atol=1e-12)
        np.testing.assert_allclose(info['probe_jackknife'], observed_jk, rtol=2e-11, atol=1e-10)
        np.testing.assert_allclose(independent_budget(L, targets, genes),
            info['B_v4a_profile'], rtol=2e-11, atol=1e-12)
        # Independent bisection uses the archived probe values (already checked
        # independently) to avoid converting tiny summation errors to a new rule.
        alpha, independent_bf = independent_solve(L, Lc, targets, genes, info['claim_est_probe'])
        np.testing.assert_allclose(alpha, info['alpha'], rtol=2e-10, atol=1e-12)
        np.testing.assert_allclose(independent_bf, bf, rtol=2e-10, atol=1e-13)
        achieved = independent_budget(np.log1p(50000*bf), targets, genes)
        np.testing.assert_allclose(achieved, claim, rtol=2e-10, atol=1e-12)
        # Final-profile emission at two fixed rows; audit all fitted integer
        # column totals as well as byte equality with the original emitter.
        final_checks = []
        recipes.v4.atlas.repair_bulk_totals = audit_repair
        try:
            for i in range(2):
                actual = original_emit(q, d, bf[i], entry['destination'], targets[i])
                expected = model.emit(q, None, bf[i], d, model.atlas.seed_for(
                    entry['destination']+':'+targets[i], 20260910), 'anchor_single_mean')
                assert np.array_equal(actual.toarray(), expected)
                final_checks.append(targets[i])
        finally:
            recipes.v4.atlas.repair_bulk_totals = original_repair
        assert len(repair_checks) == 46
    assert fingerprint(ctrl) == before
    historical = provenance.identity(task)
    assert historical == entry['historical_implementation_identity']
    for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
        old = rt.json.loads(Path(f'artifacts/scores/{task}/{recipe}_s42.json').read_text())
        assert old['implementation_sha256'] == historical
    del ctrl, q, qq, b, bf, independent_bf
    gc.collect()
    return {'status': 'PASS', 'full_panel_heads_template_depths_archived_bitwise_equal': True,
        'local_v4a_emission_bridge_targets': base_checks, 'probe_kernel_targets': probe_checks,
        'final_kernel_targets': final_checks, 'independent_probe_jackknifes': 40,
        'exact_fitted_integer_repair_checks': len(repair_checks),
        'B_independent': achieved, 'claim_independent': claim, 'alpha_independent': alpha,
        'diagnostics': info, 'input_controls_unchanged': before,
        'model_time_h5ad_reads_forbidden': True, 'input_reference_source_rows_disjoint': True,
        'historical_recipe_identity_unchanged': historical}


def main():
    assert not baseline.AUDIT.exists(), 'Preserve existing audit receipts'
    spec = baseline.read_spec()
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    assert rt.sha(index['archive']) == index['archive_sha256']
    # Mixed target resolution is essential because measured axes need not
    # contain every target. Test own-target coordinate exclusion explicitly.
    targets = np.array(['a', 'missing', 'b'])
    genes = np.array(['a', 'b', 'c', 'd'])
    L = np.array([[1., 2., 3., 4.], [4., 5., 6., 7.], [2., 4., 3., 8.]])
    fixture = []
    for changed in [L, L+np.array([[999, 0, 0, 0], [0, 0, 0, 0], [0, 777, 0, 0]])]:
        exact = Fraction(0)
        for j, gene in enumerate(genes):
            rows = [(Fraction(float(changed[i, j])), Fraction(1, len(genes)-(t in set(genes))))
                    for i, t in enumerate(targets) if t != gene]
            center = sum(x*w for x, w in rows)/sum(w for x, w in rows)
            exact += sum(w*(x-center)**2 for x, w in rows)
        ours = independent_budget(changed, targets, genes)
        published = baseline.helpers.budget(changed, targets, genes, is_log=True)
        np.testing.assert_allclose(ours, float(exact), rtol=1e-14)
        # The official shifted reducer subtracts excluded outliers from its
        # accumulators. Bound that cancellation on this adversarial fixture;
        # do not apply a looser tolerance to the biological profile checks.
        bound = 16*np.finfo(float).eps*changed.size*float(np.max(np.abs(changed)))**2
        assert abs(published-float(exact)) <= bound
        fixture.append({'exact_rational': str(exact), 'independent': ours,
            'official': float(published), 'official_abs_error': abs(published-float(exact)),
            'conservative_float64_cancellation_bound': bound})
    assert fixture[0]['exact_rational'] == fixture[1]['exact_rational']
    records = {task: audit_task(task, spec) for task in baseline.TASKS}
    baseline.read_spec()
    assert rt.json.loads(Path('artifacts/calibration_gate.json').read_text())['status'] == 'CLOSED'
    assert not Path('artifacts/holdout_events.jsonl').exists()
    rt.save(baseline.AUDIT, {'utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS',
        'spec_sha256': rt.sha(baseline.SPEC), 'verifier_sha256': rt.sha(__file__),
        'tasks': records, 'frozen_metadata_members_unchanged': len(index['entries']),
        'mixed_target_resolution_fixture': 'PASS own-target exclusion and unequal row weights',
        'mixed_target_resolution_exact_rational_checks': fixture,
        'scope': 'Frozen available-input v5a LOCO; no biological scores or external pairing',
        'external_official_label': None, 'holdout_scored': False})
    print('PASS frozen v5a adapter; biological scoring now eligible within Phase1 only', flush=True)


if __name__ == '__main__':
    main()

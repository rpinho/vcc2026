"""Independent mean/variance adapter and public-kernel bridge audit.

No biological model score is used. Checks all selected source rows, scalar
inverse-variance pooling before adaptive shrinkage, control-profile arithmetic
and two predetermined complete integer emissions per context.
"""
from datetime import datetime, timezone
from pathlib import Path
from itertools import islice
import gc
import math
import tempfile
import numpy as np
import runtime as rt
import benchmark
import provenance
import sidechain_loco_baseline as baseline
from sidechain.models.adaptive_shrink import adaptive_shrink
from sidechain.submit.build import pooled_delta


def check_source(source):
    max_fc, max_var = 0.0, 0.0
    for target in source.targets:
        i = source.pos[target]
        p, c, vp, vc = [source.values[k][i] for k in baseline.FIELDS]
        fc, var = source.effect(target)
        # Independent expression order and scalar samples of the delta method.
        expected_fc = (np.log(p+1.0)-np.log(c+1.0))/math.log(2.0)
        np.testing.assert_allclose(fc, expected_fc, rtol=1e-12, atol=5e-15)
        max_fc = max(max_fc, float(np.max(np.abs(fc-expected_fc))))
        for g in sorted(set([0, len(p)//2, len(p)-1])):
            expected = (float(vp[g])/(float(p[g])+1.0)**2
                      + float(vc[g])/(float(c[g])+1.0)**2)/math.log(2.0)**2
            np.testing.assert_allclose(var[g], expected, rtol=2e-15, atol=1e-18)
            max_var = max(max_var, abs(float(var[g])-expected))
    assert source.effect('__definitely_missing_target__') is None
    return {'source_rows': len(source.targets), 'measured_genes': len(source.genes),
            'maximum_fc_reexpression_error': max_fc, 'maximum_variance_scalar_error': max_var}


def check_pool(sources, target, genes):
    # Independent coordinate alignment and scalar sum over sources. Adaptive
    # shrinkage itself is the unmodified separately-tested public kernel.
    indexes = [{name: i for i, name in enumerate(src.genes)} for src in sources]
    effects = [src.effect(target) for src in sources]
    expected, variances, active = np.zeros(len(genes)), np.zeros(len(genes)), []
    for g, gene in enumerate(genes):
        numerator = denominator = 0.0
        for index, effect in zip(indexes, effects):
            if effect is None or gene not in index:
                continue
            fc, var = effect
            j = index[gene]
            w = 1.0/max(float(var[j]), 1e-6)
            numerator += float(fc[j])*w
            denominator += w
        if denominator:
            expected[g] = numerator/denominator
            variances[g] = 1.0/denominator
            active.append(g)
    active = np.array(active, dtype=int)
    expected[active] = adaptive_shrink(expected[active], variances[active])
    actual = pooled_delta(target, sources, genes, shrinkage=True, var_floor='none',
                          shrink_stage='pooled', shrink_rule='adaptive', shrink_k=1.0)
    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-14)
    return {'target': target, 'active_readouts': len(active),
            'maximum_pooling_error': float(np.max(np.abs(actual-expected)))}


def failure_fixture():
    with tempfile.TemporaryDirectory(dir='artifacts', prefix='sidechain_adapter_fixture_') as tmp:
        p = Path(tmp)/'missing.npz'
        np.savez(p, arithmetic_cpm_observed=True, measured_genes=np.array([True]),
                 genes=np.array(['g']), targets=np.array(['t']), n_cells=np.array([20]),
                 target_mean_cpm=np.ones((1, 1)))
        try:
            baseline.MomentSource(p, ['t'])
        except KeyError:
            return 'PASS required missing moments rejected'
        raise AssertionError('Missing variance was silently synthesized')


def main():
    spec = baseline.read_spec()
    failure = failure_fixture()
    frozen = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, rec in frozen['entries'].items():
        assert rt.sha(path) == rec['sha256'], path
    assert rt.sha(frozen['archive']) == frozen['archive_sha256']
    records = {}
    for task in baseline.TASKS:
        print('Checking Sidechain available-moment bridge', task, flush=True)
        entry = spec['tasks'][task]
        before = provenance.identity(task)
        contract, ctrl, real = benchmark.setup(task)
        assert contract['split'] == 'search'
        del real
        genes = ctrl.var_names.to_numpy(str)
        sources = [baseline.MomentSource(Path(p), entry['targets']) for p in entry['sources']]
        source_checks = {p: check_source(s) for p, s in zip(entry['sources'], sources)}
        pool_checks = [check_pool(sources, t, genes) for t in entry['targets'][:2]]
        del sources
        profile = baseline.control_profile(ctrl, entry['destination'])
        x = ctrl.X.astype(np.float64)
        lib = np.asarray(x.sum(axis=1)).ravel()
        keep = lib > 1000.0
        from scipy import sparse
        mean = np.asarray((sparse.diags(1.0/lib[keep]) @ x[keep]).mean(axis=0)).ravel()
        summed = np.asarray(x[keep].sum(axis=0)).ravel()
        np.testing.assert_allclose(profile.fraction, mean/mean.sum(), rtol=2e-13, atol=1e-16)
        assert np.array_equal(profile.bulk_fraction, summed/summed.sum())
        assert np.array_equal(profile.libsizes, lib[keep])
        assert profile.n_cells == int(keep.sum())
        difference = profile.cells - ctrl.X[keep]
        difference.eliminate_zeros()
        assert difference.nnz == 0
        # Compare the entire yielded path to direct released-emitter invocation.
        # Pooling is independent per target. Keep the complete source-row audit
        # above, then evaluate exactly the two predetermined emission cases.
        import copy
        pair_spec = copy.deepcopy(spec)
        pair_spec['tasks'][task]['targets'] = entry['targets'][:2]
        vectors = baseline.shifts(task, pair_spec, genes)
        actual = baseline.predictions(ctrl, task, pair_spec)
        comparisons = []
        for target, block in islice(actual, 2):
            em = baseline.PoissonEmitter(profile, seed=baseline.seed_for(task, target),
                lam=0.5, bulk_anchor='pooled', libsize_quantiles=(0.0, 1.0))
            # Capture what the unchanged fit actually asks the integer repair
            # to pin. Its requested profile may be projected and the fitted
            # moments have a tolerance; neither equals exact integer totals.
            import sidechain.models.count_emitters as emit_module
            old_repair = emit_module.repair_bulk_totals
            fitted = {}
            def capture_repair(integer, expected):
                fitted['columns'] = expected.sum(axis=0, dtype=np.float64)
                fitted['row_totals'] = integer.sum(axis=1, dtype=np.int64)
                return old_repair(integer, expected)
            emit_module.repair_bulk_totals = capture_repair
            try:
                original = em.emit_dual(400, vectors[target], vectors[target], shape=True,
                    scatter=None, iterations=100, tolerance=2e-4, max_projection=0.03,
                    on_fail='raise')
            finally:
                emit_module.repair_bulk_totals = old_repair
            difference = block-original
            difference.eliminate_zeros()
            assert difference.nnz == 0
            total = int(block.sum(dtype=np.int64))
            expected = fitted['columns']
            wanted = np.floor(expected).astype(np.int64)
            remaining = total-int(wanted.sum())
            assert 0 <= remaining <= len(genes)
            if remaining:
                wanted[np.argsort(expected-wanted, kind='stable')[-remaining:]] += 1
            columns = np.asarray(block.sum(axis=0, dtype=np.int64)).ravel()
            assert np.array_equal(columns, wanted)
            assert np.array_equal(np.asarray(block.sum(axis=1, dtype=np.int64)).ravel(),
                                  fitted['row_totals'])
            column_error = float(np.max(np.abs(columns-expected)))
            assert column_error <= 1.00001, column_error
            requested = em._fraction(vectors[target], bulk=True)
            fitted_bulk_error = float(np.abs(expected/total-requested).sum())
            assert fitted_bulk_error <= 0.03100001, fitted_bulk_error
            comparisons.append({'target': target, 'direct_release_counts_bitwise_equal': True,
                                 'integer_totals_equal_fitted_largest_remainder': True,
                                 'maximum_pinned_fitted_column_error_counts': column_error,
                                 'fitted_vs_requested_bulk_L1': fitted_bulk_error,
                                 'dual': em.last_dual, 'shaped': bool(em.last_shaped)})
        actual.close()
        for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
            old = rt.json.loads(Path(f'artifacts/scores/{task}/{recipe}_s42.json').read_text())
            assert old['implementation_sha256'] == before
        assert provenance.identity(task) == before
        records[task] = {'status': 'PASS', 'source_moments': source_checks,
                         'scalar_pooling_before_public_adaptive': pool_checks,
                         'input_control_profile': 'PASS independent arithmetic/cells/depths',
                         'predetermined_emission_checks': comparisons,
                         'historical_recipe_identity_unchanged': before}
        del ctrl, x, profile, vectors, actual, original, block, em
        gc.collect()
    baseline.read_spec()
    rt.save(baseline.AUDIT, {'utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS',
        'spec_sha256': rt.sha(baseline.SPEC), 'verifier_sha256': rt.sha(__file__),
        'tasks': records, 'missing_moment_fixture': failure,
        'frozen_metadata_members_unchanged': len(frozen['entries']),
        'external_official_label': None, 'holdout_scored': False,
        'scope': 'Observed float32-moment adapter and released pooling/shape kernels only; NOT exact SER17'})
    print('PASS available-moment Sidechain bridge; NOT SER17', flush=True)


if __name__ == '__main__':
    main()

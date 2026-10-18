"""Independent available-input reproduction check before Phase1 baseline scoring.

Compare the standalone unmodified release against the separately archived v4
atlas_adapted path, full source centering and every head on both SEARCH panels.
The release sums dense control-normalized blocks; v4 uses sparse reductions.
CD4's archived path also executes a float32 multiply/divide for an absent CD4
vote. Audit that operation exactly, then bypass it ONLY in the comparison path.
With identical control summaries, heads must match bitwise. Actual released vs
sparse-control heads retain the original float64 tolerance (2e-12 relative,
2e-15 absolute). Templates/coverage remain bitwise checks. Record count
differences on two predetermined emissions against the unchanged archived path. Neither
path reads reference responses for prediction. No benchmark identity changes.
"""
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path
import gc
import hashlib
import runtime as rt
import numpy as np
import atlas_loco_baseline as baseline
import benchmark
import provenance
import recipes


def main():
    spec = baseline.read_spec()
    records = {}
    frozen = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, entry in frozen['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    assert rt.sha(frozen['archive']) == frozen['archive_sha256']
    assert rt.json.loads(Path('artifacts/calibration_gate.json').read_text())['status'] == 'CLOSED'
    for task in baseline.TASKS:
        print('Checking standalone AtlasShift against archived adaptation', task, flush=True)
        m, ctrl, real = benchmark.setup(task)
        del real  # Predictions receive the input controls only.
        before = provenance.identity(task)
        q, depths, desired, coverage = baseline.heads(ctrl, task, spec)
        v4 = recipes.v4
        old_banks, old_load, old_align, old_cd4 = v4.BANKS, v4.load_banks, v4.atlas.aligned_cd4, v4.cd4_family
        old_add = v4.atlas.add_cd4_family
        absent_vote_checks = []
        try:
            if m['context'] == 'CD4':
                v4.BANKS = rt.ROOT/'artifacts/cd4_banks'
                v4.cd4_family = lambda *a, **k: None
                v4.atlas.aligned_cd4 = lambda path, targets, genes, *a, **k: (
                    np.zeros((len(targets), len(genes)), np.float32),
                    np.zeros((len(targets), len(genes)), bool))
            v4.load_banks = lambda *a, **k: ['K562_GWPS_raw', 'HCT116']
            qq, dd, mean, bulk = v4.template(ctrl, ctrl.var_names.to_numpy(str), seed=20260910)
            reference, reference_coverage = v4.heads(np.array(m['targets']), ctrl.var_names.to_numpy(str),
                                                     mean, bulk, 'anchor', exclude=m['context'], stats={})
            archived_reference = reference
            if m['context'] == 'CD4':
                def absent_vote_bridge(effect, covered, cd4, available, weight, **kwargs):
                    assert not available.any() and not cd4.any() and weight == .5
                    rounded, denominator = old_add(effect, covered, cd4, available, weight, **kwargs)
                    independent = np.divide(effect * covered, covered,
                                            out=np.zeros_like(effect), where=covered > 0)
                    assert np.array_equal(rounded, independent)
                    assert np.array_equal(denominator, covered)
                    absent_vote_checks.append({
                        'space': kwargs['space'], 'effect_dtype': str(effect.dtype),
                        'coverage_dtype': str(covered.dtype),
                        'entries_changed_by_absent_vote_roundtrip': int(np.count_nonzero(rounded != effect)),
                        'max_abs_effect_change': float(np.max(np.abs(rounded-effect))),
                        'independent_float32_roundtrip_bitwise_equal': True})
                    return effect, covered
                v4.atlas.add_cd4_family = absent_vote_bridge
                reference, reference_coverage = v4.heads(np.array(m['targets']), ctrl.var_names.to_numpy(str),
                                                         mean, bulk, 'anchor', exclude=m['context'], stats={})
                assert len(absent_vote_checks) == 2
            assert np.array_equal(q, qq) and np.array_equal(depths, dd)
            # Isolate control accumulation from every later operation: the
            # released runner with the archived summaries must match exactly.
            _, published = baseline.release_modules()
            old_template = published.control_template
            try:
                published.control_template = lambda *a, **k: (qq, dd, mean, bulk)
                _, _, common_control_heads, common_coverage = baseline.heads(ctrl, task, spec)
            finally:
                published.control_template = old_template
            head_checks = []
            for index in range(2):
                assert np.array_equal(common_control_heads[index], reference[index]), (task, index, 'common-control heads')
                assert np.array_equal(common_coverage[index], reference_coverage[index])
                np.testing.assert_allclose(desired[index], reference[index], rtol=2e-12, atol=2e-15)
                assert np.array_equal(coverage[index], reference_coverage[index]), (task, index, 'coverage')
                head_checks.append({'head': ['mean_cpm', 'pseudobulk'][index],
                    'bitwise_equal': bool(np.array_equal(desired[index], reference[index])),
                    'max_abs_error': float(np.max(np.abs(desired[index]-reference[index]))),
                    'same_control_summary_heads_bitwise_equal': True,
                    'original_archived_head_max_abs_difference': float(np.max(np.abs(desired[index]-archived_reference[index]))),
                    'relative_tolerance': 2e-12, 'absolute_tolerance': 2e-15})
        finally:
            v4.BANKS, v4.load_banks, v4.atlas.aligned_cd4, v4.cd4_family = old_banks, old_load, old_align, old_cd4
            v4.atlas.add_cd4_family = old_add
        del qq, q, reference, archived_reference, common_control_heads, common_coverage, desired, coverage
        gc.collect()
        actual = baseline.predictions(ctrl, task, spec)
        expected = recipes.predict(ctrl, m['targets'], m['context'], 'atlas_adapted', seed=42)
        emitted = []
        for (ta, xa), (tb, xb) in islice(zip(actual, expected, strict=True), 2):
            xa, xb = xa.tocsr(), xb.tocsr()
            assert ta == tb and xa.shape == xb.shape == (400, ctrl.n_vars)
            difference = xa.astype(np.int64)-xb.astype(np.int64)
            difference.eliminate_zeros()
            assert np.array_equal(xa.sum(axis=1), xb.sum(axis=1)), (task, ta, 'depths')
            assert np.isfinite(xa.data).all() and np.all(xa.data >= 0)
            assert np.equal(xa.data, np.floor(xa.data)).all()
            digest = hashlib.sha256()
            for value in [xa.data, xa.indices, xa.indptr]:
                digest.update(value.tobytes())
            emitted.append({'target': ta, 'prediction_sha256': digest.hexdigest(),
                            'bitwise_equal': difference.nnz == 0,
                            'count_entries_differing': difference.nnz,
                            'max_abs_count_difference': int(np.max(np.abs(difference.data))) if difference.nnz else 0})
        actual.close()
        expected.close()
        benchmark.release_models()
        assert provenance.identity(task) == before
        for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
            prior = rt.json.loads(Path(f'artifacts/scores/{task}/{recipe}_s42.json').read_text())
            assert prior['implementation_sha256'] == before
        records[task] = {'status': 'PASS', 'full_panel_head_checks': head_checks,
                         'absent_cd4_vote_auditor_bridge': absent_vote_checks,
                         'control_template_and_depths_bitwise_equal': True,
                         'two_predetermined_emission_comparisons': emitted,
                         'historical_recipe_identity_unchanged': before}
        del actual, expected, ctrl
        gc.collect()
    rt.save('artifacts/atlas_loco_audit.json', {
        'utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS',
        'spec_sha256': rt.sha(baseline.SPEC), 'verifier_sha256': rt.sha(__file__),
        'tasks': records, 'frozen_metadata_members_unchanged': len(frozen['entries']),
        'scope': 'Unchanged released kernels match the separately archived available-input adaptation; shared source preparation is not independently reconstructed.',
        'external_official_label': None, 'holdout_scored': False})
    print('PASS available-input AtlasShift reproduction', flush=True)


if __name__ == '__main__':
    main()

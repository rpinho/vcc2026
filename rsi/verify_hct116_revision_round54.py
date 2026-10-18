"""Independent archived-call-path and source-centering checks before scoring."""
from datetime import datetime, timezone
from pathlib import Path
import gc
import hashlib
import runtime as rt
import numpy as np
from scipy import sparse
import hct116_revision_round54 as revision
import recipes


def digest(x):
    return hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()


def controls(task):
    # Read only the frozen input controls; no destination reference response.
    import anndata as ad
    path = Path('artifacts/tasks')/(task+'.h5ad')
    assert rt.sha(path) == revision.load(path.with_suffix('.json'))['sha256']
    a = ad.read_h5ad(path, backed='r')
    try:
        return a[a.obs.evaluation_role.astype(str).to_numpy() == 'input_control'].to_memory()
    finally:
        a.file.close()


def original_call(ctrl, entry, recipe):
    v4 = recipes.v4
    template, heads, emit = v4.template, v4.heads, v4.emit
    saved = {}
    def capture_template(*a, **k):
        value = template(*a, **k)
        saved['template'] = value
        return value
    def capture_heads(*a, **k):
        value = heads(*a, **k)
        saved['heads'] = value
        return value
    v4.template, v4.heads = capture_template, capture_heads
    generator = recipes.predict(ctrl, entry['targets'], entry['context'], recipe, 42)
    try:
        target, matrix = next(generator)
        assert target == entry['targets'][0]
        saved['emission'] = matrix.toarray()
    finally:
        generator.close()
        v4.template, v4.heads, v4.emit = template, heads, emit
        v4._POOL_CACHE.clear()
    return saved


def center_audit(bank, genes, panel):
    atlas = recipes.v4.atlas
    source = atlas.load_xatlas(bank, prior_counts=100000)
    assert len(source.targets) == 509 and len(set(source.targets)) == 509
    assert 'HMGXB4' not in source.targets
    result = {}
    for space in ['log2fc', 'bulk_delta']:
        raw, mask = atlas.aligned_effect(source, source.targets, genes, space=space, common_subtract=0)
        centered, cmask = atlas.aligned_effect(source, source.targets, genes, space=space, common_subtract=1)
        np.testing.assert_array_equal(mask, cmask)
        # Independent per-gene float64 reduction, with original self exclusion.
        reference = raw.astype(float)
        eligible = mask.copy()
        eligible &= source.targets[:, None] != genes[None, :]
        denominator = eligible.sum(0)
        common = np.divide(np.where(eligible, reference, 0).sum(0), denominator,
                           out=np.zeros(len(genes)), where=denominator > 0)
        reference = np.where(mask, reference-common, 0)
        np.testing.assert_allclose(centered, reference, rtol=3e-5, atol=3e-6)
        # The source threshold is unchanged, including every old low-count row.
        rows = {t:i for i,t in enumerate(source.targets)}
        expected = np.array([source.n_cells[rows[t]] >= 20 if t in rows else False for t in panel])
        aligned, available = atlas.aligned_effect(source, np.array(panel), genes, space=space, common_subtract=0)
        np.testing.assert_array_equal(available.any(1), expected)
        assert sum(expected) == 78 and not expected[panel.index('HMGXB4')]
        full, cov = atlas.fuse_source_centered([source], [1], np.array(panel), genes,
                                               space=space, common_subtract=1)
        for i,t in enumerate(panel):
            if t in rows:
                np.testing.assert_array_equal(full[i], centered[rows[t]])
        panel_centered, _ = atlas.aligned_effect(source, np.array(panel), genes, space=space, common_subtract=1)
        separation = float(np.max(np.abs(full-panel_centered)))
        assert separation > 1e-5, 'Fixture must distinguish source from panel centering'
        result[space] = dict(source_population=509, eligible_source_rows=int(np.sum(source.n_cells >= 20)),
            float64_reference_max_abs_error=float(np.max(np.abs(centered-reference))),
            panel_centering_max_abs_difference=separation, centered_sha256=digest(centered),
            HMGXB4_abstains=True, source_wide_fusion_bitwise=True)
    return result


def main():
    assert not (revision.ROOT/'audit.json').exists()
    plan = revision.check_plan()
    assert __import__('official_scoring').provenance() == plan['official_environment']
    # Prior moment auditor certified raw union; replay all current bindings.
    from audit_hct116_extension_round53 import bindings
    bindings()
    records = {}
    for task in revision.TASKS:
        entry = plan['tasks'][task]
        ctrl = controls(task)
        assert ctrl.n_obs >= 1600 and ctrl.n_vars > 17000
        old_folder = rt.ROOT/'artifacts/cd4_banks' if entry['context'] == 'CD4' else recipes.v4.BANKS
        checks = {}
        for recipe in revision.RECIPES:
            print('Archived-path equality',task,recipe,flush=True)
            old = original_call(ctrl, entry, recipe)
            q,d,desired,cov = revision.heads(ctrl, entry['targets'], entry['context'], recipe, old_folder)
            for a,b in zip([q,d],old['template'][:2]):
                np.testing.assert_array_equal(a,b)
            for a,b in zip(desired,old['heads'][0]):
                np.testing.assert_array_equal(a,b)
            for a,b in zip(cov,old['heads'][1]):
                np.testing.assert_array_equal(a,b)
            target = entry['targets'][0]
            counts = recipes.v4.emit(q, desired[0][0], desired[1][0], d,
                recipes.seed_for(entry['context']+':'+target,20260910),revision.VARIANTS[recipe])
            np.testing.assert_array_equal(counts, old['emission'])
            checks[recipe] = dict(all79_heads_and_coverage_bitwise=True, template_depths_bitwise=True,
                first_emission_bitwise=True, first_target=target, emission_sha256=digest(counts))
            del old,q,d,desired,cov,counts
            gc.collect()
        records[task] = checks
        del ctrl
        gc.collect()
    ctrl = controls('CD4_depth192_search')
    centering = center_audit(revision.ROOT/'banks/CD4_depth192_search/HCT116.npz',
                             ctrl.var_names.to_numpy(str),plan['tasks']['CD4_depth192_search']['targets'])
    assert all(__import__('provenance').identity(t) == h for t,h in plan['original_implementation_sha256'].items())
    revision.check_plan()
    rt.save(revision.ROOT/'audit.json',dict(status='PASS archived paths, source-wide centering and source isolation',
        utc=datetime.now(timezone.utc).isoformat(),plan_sha256=rt.sha(revision.PLAN),
        verifier_sha256=rt.sha(__file__),archived_path_checks=records,centering=centering,
        bank_sha256=rt.sha('artifacts/HCT116_extension_round53/bank.npz'),
        original_metadata_members_unchanged=38, historical_scores_unchanged=10,
        no_destination_response_used_for_prediction=True,holdout_scored=False,
        phase2_allowed=False,scope='Numerical adapter equality, not a biological gain or exact external official-score reproduction'))
    print('PASS revision prerequisites',flush=True)


if __name__ == '__main__':
    main()

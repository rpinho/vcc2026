"""Replay the archived v4b kernel on SEARCH controls, bypassing recipes.predict.

No destination response, metric, or HOLDOUT outcome enters model construction.
The LOCO source exclusions and context seed mapping are held fixed, not searched.
"""
from datetime import datetime, timezone
import gc
import hashlib
import importlib.util
from pathlib import Path

import runtime as rt
import numpy as np
from scipy import sparse


def main():
    task = 'H1_rich_search'
    contract = rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())[task]
    assert contract['split'] == 'search' and contract['context'] == 'H1'
    record_path = Path('artifacts/scores')/task/'v4b_s42.json'
    record = rt.json.loads(record_path.read_text())
    import provenance
    assert provenance.identity(task) == record['implementation_sha256']
    historical = rt.ROOT.parent/'model_v4'
    manifest = rt.json.loads((historical/'out/prediction_v4b.manifest.json').read_text())
    identity = rt.json.loads((historical/'artifacts/source_identity.json').read_text())
    files = {}
    for current, archived, expected in [
        ('vendor/v4_model.py', historical/'model.py', manifest['model_sha256']),
        ('vendor/atlas_model.py', historical/'vendor/atlas_model.py', identity['code']['atlas_model.py']),
        ('vendor/adaptive_shrink.py', historical/'vendor/adaptive_shrink.py', identity['code']['adaptive_shrink.py'])]:
        assert rt.sha(current) == rt.sha(archived) == expected
        files[current] = expected
    assert manifest['variant'] == 'iv_ash_dual_guide' and manifest['seed'] == 20260910
    spec = importlib.util.spec_from_file_location('archived_v4b_replay', historical/'model.py')
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    shared = rt.ROOT.parent/'public_data/model_v4_signature_banks'
    original.BANKS = shared
    # This is the documented H1 SEARCH adaptation, rather than the hidden-panel
    # source pool: H1 destination and globally sealed HEK293T/RPE1 excluded.
    names = ['K562_GWPS_raw', 'HCT116']
    original.load_banks = lambda *args, **kwargs: names.copy()
    for name in names + ['CD4_Rest', 'CD4_Stim48hr', 'CD4_Stim8hr']:
        path = shared/(name+'.npz')
        expected = next(r['sha256'] for r in identity['banks'] if r['source'] == name)
        assert rt.sha(path) == expected
        files['shared_sources/'+path.name] = expected
    files['shared_sources/CD4_DE_statistics.npz'] = rt.sha(shared/'CD4_DE_statistics.npz')
    files['artifacts/promoter_pairs.csv'] = rt.sha('artifacts/promoter_pairs.csv')
    import task_io
    control, unused_real = task_io.split(Path('artifacts/tasks')/(task+'.h5ad'))
    del unused_real
    gc.collect()
    genes = control.var_names.to_numpy(str)
    targets = np.array(contract['targets'])
    q, depths, mean, bulk = original.template(control, genes, seed=20260910)
    del control
    shrink_calls = []
    adaptive = original.adaptive_shrink

    def observe(x, variance, **kwargs):
        out = adaptive(x, variance, **kwargs)
        shrink_calls.append(dict(input_norm=float(np.linalg.norm(x)),
            output_norm=float(np.linalg.norm(out)), coordinates=len(x)))
        return out

    original.adaptive_shrink = observe
    stats = {}
    (p, b), coverage = original.heads(targets, genes, mean, bulk,
        manifest['variant'], exclude='H1', stats=stats)
    print('Archived heads complete', stats, flush=True)
    original._POOL_CACHE.clear()
    digest = hashlib.sha256()
    counts = 0
    for i, target in enumerate(targets):
        rt.guard(512*1024**2)
        x = sparse.csr_matrix(original.emit(q, p[i], b[i], depths,
            original.atlas.seed_for('H1:'+target, 20260910), manifest['variant']))
        x.eliminate_zeros()
        assert x.shape == (400, len(genes)) and np.array_equal(np.asarray(x.sum(1)).ravel(), depths)
        for value in [x.data, x.indices, x.indptr]:
            digest.update(value.tobytes())
        digest.update(target.encode())
        counts += 1
        if counts % 10 == 0:
            print('Archived v4b replay', counts, '/', len(targets), flush=True)
    assert counts == record['n_targets'] == 79
    assert digest.hexdigest() == record['prediction_sha256'], (digest.hexdigest(), record['prediction_sha256'])
    ratios = [r['output_norm']/r['input_norm'] for r in shrink_calls if r['input_norm'] > 0]
    rt.save('artifacts/v4b_port_audit.json', dict(
        status='PASS archived v4b kernel and full SEARCH count stream match',
        utc=datetime.now(timezone.utc).isoformat(), verifier_sha256=rt.sha(__file__),
        task=task, targets=counts, prediction_sha256=digest.hexdigest(),
        score_record_sha256=rt.sha(record_path), implementation_sha256=record['implementation_sha256'],
        historical_manifest_sha256=rt.sha(historical/'out/prediction_v4b.manifest.json'),
        files=files, source_names=names+['CD4 family'], seed=20260910,
        adaptive_stats=stats, shrink_vector_norm_ratio_min_median_max=[float(f(ratios)) for f in [np.min, np.median, np.max]],
        mean_head_L1_from_control=float(np.abs(p-mean).sum(1).mean()),
        bulk_head_L1_from_control=float(np.abs(b-bulk).sum(1).mean()),
        holdout_scored=False,
        limitation='Tests the documented LOCO adaptation, not hidden-panel equivalence. Destination control distributions, readout axes and available source pool differ from official validation. No result-based model correction or rank repair.'))
    print('PASS all 79 v4b counts match the completed score digest', flush=True)


if __name__ == '__main__':
    main()

"""Independent source/mask/norm/inversion/emission audit, no response reads.

Uses separate raw-bank arithmetic, full-panel comparisons, adversarial NaN and
unseen-target fixtures, direct RNG replay and selected unmodified release tests.
"""
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path
import gc
import hashlib
import subprocess
import sys
import os

import runtime as rt
import numpy as np
import pandas as pd
from scipy import sparse
import zsp_loco_baseline as baseline


def independent_sources(union):
    """Reconstruct LFC from native pooled counts before union projection."""
    result = []
    for name, path in baseline.SOURCES.items():
        with np.load(path, allow_pickle=False) as d:
            if name == 'hct116':
                mask = d['measured_genes'].astype(bool)
                genes = d['genes'][mask].astype(str)
                raw = d['target_count_sums'][:, mask].astype(float)
                control = d['global_control_probability'][mask].astype(float)
            else:
                genes = d['genes'].astype(str)
                raw = d['raw_means'][1:].astype(float)
                control = d['raw_means'][0].astype(float)
            keep = d['n_cells'] >= 10
            targets = d['targets'][keep].astype(str)
            raw = raw[keep]
            # Sum normalization cancels target_count_sums' cell-count factor.
            expected = np.log2((10000*raw/raw.sum(axis=1)[:, None]+0.1) /
                              (10000*control/control.sum()+0.1))
            generic = np.median(expected, axis=0)
            aligned = np.full((len(targets), len(union)), np.nan)
            positions = union.get_indexer(genes)
            assert (positions >= 0).all()
            aligned[:, positions] = expected
            gen = np.full(len(union), np.nan)
            gen[positions] = generic
            result.append((name, targets, aligned, gen))
    return result


def independent_prediction(sources, targets, width):
    stack = np.full((len(sources), len(targets), width), np.nan)
    for i, (_, labels, lfc, _) in enumerate(sources):
        lookup = {t: j for j, t in enumerate(labels)}
        for j, target in enumerate(targets):
            if target in lookup:
                stack[i, j] = lfc[lookup[target]]
    finite = np.isfinite(stack)
    count = finite.sum(axis=0)
    prediction = np.divide(np.nansum(stack, axis=0), count,
        out=np.zeros((len(targets), width)), where=count > 0)
    gen = np.stack([s[3] for s in sources])
    gcnt = np.isfinite(gen).sum(axis=0)
    fallback = np.divide(np.nansum(gen, axis=0), gcnt,
        out=np.zeros(width), where=gcnt > 0)
    for j in range(len(targets)):
        seen = finite[:, j].any(axis=1)
        if not seen.any():
            prediction[j] = fallback
            continue
        norms = [np.sqrt(np.nansum(stack[i, j]**2)) for i in np.flatnonzero(seen)]
        radius = sum(norms)/len(norms)
        norm = np.sqrt(np.square(prediction[j]).sum())
        if norm > 0:
            prediction[j] *= radius/norm
    return prediction


def synthetic_fixtures():
    _, model, _ = baseline.kernels()
    genes = pd.Index(['A', 'B', 'C', 'D'])
    src = []
    # Includes measured zero, gene absent from one source, exact cancellation,
    # per-source-only targets, partial readout gaps, and wholly unseen target.
    fixtures = [(['P', 'ONLY1'], [[1., -2., 0., np.nan], [2., 0., 3., np.nan]], [1.5,-1.,1.5,np.nan]),
                (['P', 'ONLY2'], [[-1., 2., 4., 2.], [0., -3., 1., 0.]], [-.5,-.5,2.5,1.])]
    independent = []
    for i, (labels, vals, gen) in enumerate(fixtures):
        frame = pd.DataFrame(vals, index=labels, columns=genes)
        generic = np.asarray(gen)
        src.append(model.SourceEffects(str(i), frame, generic, frame-generic,
                                       np.arange(1., 5.), pd.Series(dtype=float)))
        independent.append((str(i), np.asarray(labels), np.asarray(vals), generic))
    targets = ['P','ONLY1','ONLY2','UNSEEN']
    got = model.NormRestoredTransfer(np.inf)(src, np.arange(1.,5.), targets, genes).lfc.to_numpy()
    expected = independent_prediction(independent, targets, 4)
    np.testing.assert_allclose(got, expected, rtol=1e-13, atol=1e-13)
    assert got[0,0] == got[0,1] == 0
    assert got[1,3] == 0  # readout absent for this target is not generic-filled
    assert got[3,3] == 1  # entirely unseen target gets generic fallback
    assert got[0,2] > 2  # norm is restored after consensus dilution
    return {'targets': targets, 'prediction': got.tolist(), 'status': 'PASS'}


def main():
    spec = baseline.read_spec()
    assert rt.json.loads(Path('artifacts/calibration_gate.json').read_text())['status'] == 'CLOSED'
    frozen = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, entry in frozen['entries'].items():
        assert rt.sha(path) == entry['sha256']
    assert rt.sha(frozen['archive']) == frozen['archive_sha256']
    receipt = rt.json.loads(Path('artifacts/zsp_release_receipt.json').read_text())
    tree = rt.json.loads(Path('../research/sources/MMatinGerami__zero-shot-perturbation__tree.json').read_text())
    assert receipt['commit'] == tree['sha'] == '2a198d5fa165179633ed0c0becc0d2b050b27187'
    blobs = {i['path']: i for i in tree['tree'] if i['type'] == 'blob'}
    for path, item in receipt['files'].items():
        data = (baseline.RELEASE/path).read_bytes()
        git = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        assert git == blobs[path]['sha'] == item['git_blob']
    # Run only released tests of kernels used here; no remote I/O or models fit.
    env = dict(os.environ, PYTHONPATH=str((baseline.RELEASE/'src').resolve()))
    selection = 'log_fold_change or cp10k or align_intersects or decompose_removes or mean_transfer_falls or lfc_to_counts_roundtrip or emit_cells_matches or norm_restored_and_median'
    test = subprocess.run([sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider',
        'tests/test_data_models.py', '-k', selection], cwd=baseline.RELEASE,
        env=env, capture_output=True, text=True)
    Path('artifacts/zsp_release_kernel_tests_round41.txt').write_text(test.stdout+test.stderr)
    print(test.stdout, flush=True)
    assert test.returncode == 0, test.stderr
    fixture = synthetic_fixtures()
    contexts = baseline.source_contexts()
    union = contexts[0].genes
    sources = independent_sources(union)
    _, model, _ = baseline.kernels()
    for ctx, (_, targets, expected, generic) in zip(contexts, sources, strict=True):
        released = model.decompose(ctx, min_cells=10)
        assert list(released.lfc.index) == list(targets)
        np.testing.assert_allclose(released.lfc.to_numpy(), expected, rtol=2e-11, atol=2e-12, equal_nan=True)
        np.testing.assert_allclose(released.generic, generic, rtol=2e-11, atol=2e-12, equal_nan=True)
    del contexts, released
    records = {}
    for task in baseline.TASKS:
        print('Independent full-panel prediction audit', task, flush=True)
        ctrl = baseline.input_controls(task)
        assert set(ctrl.obs.evaluation_role.astype(str)) == {'input_control'}
        diagnostics = {}
        mean, got = baseline.heads(ctrl, task, spec, diagnostics)
        predicted = independent_prediction(sources, spec['tasks'][task]['targets'], len(union))
        panel = pd.Index(ctrl.var_names.astype(str))
        pos = union.get_indexer(panel)
        lfc = np.zeros((len(predicted),len(panel)))
        lfc[:, pos >= 0] = predicted[:, pos[pos >= 0]]
        basal = mean.astype(float)/mean.astype(float).sum()*10000
        expected = np.maximum((basal+.1)*np.exp2(lfc)-.1, 0)*(mean.sum()/10000)
        np.testing.assert_allclose(got, expected, rtol=3e-11, atol=3e-9)
        # Independently replay two predetermined emissions from seed0.
        rng = np.random.default_rng(0)
        emissions = []
        iterator = baseline.predictions(ctrl, task, spec)
        for i, (target, x) in enumerate(islice(iterator, 2)):
            ratio = np.ones_like(got[i])
            valid = mean > 0
            ratio[valid] = got[i,valid]/mean[valid]
            chosen = rng.integers(0, ctrl.n_obs, 400)
            raw = sparse.csr_matrix(ctrl.X)[chosen].toarray().astype(float)*ratio
            floor = np.floor(raw)
            direct = (floor+(rng.random(raw.shape)<raw-floor)).astype(np.float32)
            assert np.array_equal(x.toarray(), direct)
            emissions.append({'target': target, 'sha256_dense_float32': hashlib.sha256(direct.tobytes()).hexdigest(),
                'integer_nonnegative': True, 'RNG_replay_bitwise_equal': True})
        iterator.close()
        records[task] = {'status': 'PASS', 'input_controls': ctrl.n_obs,
            'source_diagnostics': diagnostics, 'full_mean_max_abs_error': float(np.abs(got-expected).max()),
            'rtol': 3e-11, 'atol': 3e-9, 'emissions': emissions}
        del ctrl, got, expected, predicted, lfc, raw, direct
        gc.collect()
    baseline.read_spec()
    rt.save(baseline.AUDIT, {'utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS',
        'spec_sha256': rt.sha(baseline.SPEC), 'verifier_sha256': rt.sha(__file__),
        'release_kernel_tests_sha256': rt.sha('artifacts/zsp_release_kernel_tests_round41.txt'),
        'synthetic_fixture': fixture, 'tasks': records,
        'frozen_metadata_members_unchanged': len(frozen['entries']), 'external_official_label': None,
        'holdout_scored': False, 'scope': 'Actual released kernels with four fixed available source lines; independent numerical and source/mask audit, not biological validation or exact eight-source release.'})
    print('PASS fixed four-source Gerami adapter', flush=True)


if __name__ == '__main__':
    main()

"""Check restored SEARCH identities and retained scorer certificates before reuse."""
from datetime import datetime, timezone
from pathlib import Path
import gc

import runtime as rt
import benchmark
import official_scoring as scoring


def main():
    current = scoring.provenance()
    assert current['version'] == '0.18.0'
    certificates = [
        ('artifacts/bounded_fixture.json', 'adapter_sha256', 'bounded_de.py', 'cell_eval2'),
        ('artifacts/bounded_score_fixture.json', 'driver_sha256', 'bounded_score.py', 'official'),
    ]
    for file, key, implementation, environment in certificates:
        receipt = rt.json.loads(Path(file).read_text())
        assert receipt['status'] == 'PASS bitwise equality'
        assert receipt[key] == rt.sha(implementation)
        assert receipt[environment] == current, f'Scoring environment changed: {file}'
    import stat_cache
    stat_cache.install()
    assert rt.json.loads(Path('artifacts/task_io_fixture.json').read_text())['implementation_sha256'] == rt.sha('task_io.py')
    baseline_fixture = rt.json.loads(Path('artifacts/baseline_profile_fixture.json').read_text())
    assert baseline_fixture['status'] == 'PASS bitwise profile equality and unchanged reference'
    assert baseline_fixture['builder_sha256'] == rt.sha('baseline_protocol.py')
    atlas = rt.json.loads(Path('artifacts/atlas_reproduction.json').read_text())
    assert atlas['source_sha256'] == rt.sha('vendor/atlas_model.py')
    tasks = {}
    for task in ['H1_rich_search', 'CD4_matched_search']:
        contract, control, real = benchmark.setup(task)
        assert contract['split'] == 'search' and contract['eligible_for_primary_calibration']
        reference = scoring.matrix_id(real)
        folder = Path('artifacts/scores')/task
        reused = []
        baseline = folder/'baseline_mean_response.json'
        if baseline.exists():
            description = rt.json.loads((folder/'baseline_mean_response_protocol.json').read_text())
            assert description['reference_sha256'] == reference
            assert description['builder_sha256'] in [baseline_fixture['builder_sha256']] + baseline_fixture['compatible_legacy_builder_sha256']
            reused.append(str(baseline))
        assert not (folder/'anchors.json').exists(), 'Review any retained primary anchors explicitly'
        assert not list(folder.glob('*_s[0-9]*.json')), 'Review retained completed scores explicitly'
        assert not list((Path('artifacts/raw_scores')/task).glob('*.json')), 'Review retained raw scores explicitly'
        tasks[task] = dict(reference_sha256=reference, shape=list(real.shape), input_controls=control.n_obs,
            retained_baseline=reused, implementation_sha256=__import__('provenance').identity(task))
        del control, real
        gc.collect()
    rt.save('artifacts/calibration_restore_preflight.json', dict(status='PASS',
        utc=datetime.now(timezone.utc).isoformat(), tasks=tasks, official=current,
        verifier_sha256=rt.sha(__file__), holdout_scored=False,
        limitation='Identity/certificate audit, no new metric equality fixture or model score; legacy sensitivity excluded'))
    print('PASS SEARCH calibration restore preflight', flush=True)


if __name__ == '__main__':
    main()

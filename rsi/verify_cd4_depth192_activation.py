"""Read-only protocol/setup/cache audit; computes no biological scores."""
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt
import benchmark
import official_scoring as scoring


def main():
    revision = Path('artifacts/protocol_revisions/CD4_depth192')
    activation = rt.json.loads((revision/'activation.json').read_text())
    assert activation['status'].startswith('ACTIVE count-audited protocol;')
    for path, entry in activation['after'].items():
        assert rt.sha(path) == rt.sha(entry['archive']) == entry['sha256']
    old_manifest = rt.json.loads((revision/'manifest.json.before').read_text())
    new_manifest = rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())
    assert set(new_manifest) == set(old_manifest) | {'CD4_depth192_search'}
    for task, contract in old_manifest.items():
        if task != 'CD4_depth120_search':
            assert new_manifest[task] == contract, task
    assert new_manifest['CD4_depth120_search']['eligible_for_primary_calibration'] is False
    roles = rt.json.loads(Path('artifacts/benchmark_roles.json').read_text())
    old_roles = rt.json.loads((revision/'benchmark_roles.json.before').read_text())
    assert roles['primary_search'] == ['H1_rich_search', 'CD4_depth192_search']
    assert roles['primary_holdout'] == old_roles['primary_holdout']
    for key in ['depth_sensitivity', 'depth_holdout_sensitivity', 'diagnostic_holdout']:
        assert roles[key] == old_roles[key]
    policy = rt.json.loads(Path('artifacts/anchor_policy.json').read_text())
    assert policy['CD4_depth192_search'] == 5
    assert new_manifest['CD4_depth192_search']['replicate_splits'] == 5
    current = scoring.provenance()
    for file, key, implementation, environment in [
        ('artifacts/bounded_fixture.json', 'adapter_sha256', 'bounded_de.py', 'cell_eval2'),
        ('artifacts/bounded_score_fixture.json', 'driver_sha256', 'bounded_score.py', 'official')]:
        receipt = rt.json.loads(Path(file).read_text())
        assert receipt['status'] == 'PASS bitwise equality'
        assert receipt[key] == rt.sha(implementation) and receipt[environment] == current
    import stat_cache
    stat_cache.install()
    contract, control, real = benchmark.setup('CD4_depth192_search')
    assert contract['n_real'] == 192 and real.n_obs == 15168 + 768
    assert control.n_obs == 1920 and real.n_vars == 17772
    reference = scoring.matrix_id(real)
    folder = Path('artifacts/scores/CD4_depth192_search')
    assert not folder.exists() or not list(folder.iterdir()), 'Must audit before new score/cache creation'
    import provenance
    h1_identity = provenance.identity('H1_rich_search')
    scores = list(Path('artifacts/scores/H1_rich_search').glob('*_s42.json'))
    assert len(scores) == 5
    for path in scores:
        assert rt.json.loads(path.read_text())['implementation_sha256'] == h1_identity, path
    rt.save('artifacts/CD4_depth192_activation_audit.json', dict(status='PASS',
        utc=datetime.now(timezone.utc).isoformat(), verifier_sha256=rt.sha(__file__),
        activation_sha256=rt.sha(revision/'activation.json'),
        task='CD4_depth192_search', reference_sha256=reference,
        implementation_sha256=provenance.identity('CD4_depth192_search'),
        shape=list(real.shape), input_controls=control.n_obs, anchor_splits=5,
        unchanged_H1_completed_score_identities=len(scores), holdout_contracts_unchanged=True,
        holdout_scored=False, limitation='Setup, identities and retained official equality certificates only; new scales and two-context ranks UNVERIFIED'))
    print('PASS activated CD4 setup, five splits, official certificates, five unchanged H1 score identities and HOLDOUT contracts', flush=True)


if __name__ == '__main__':
    main()

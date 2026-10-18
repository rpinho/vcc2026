"""Freeze the supplied v5b recipe extension before any new biological score."""
from datetime import datetime, timezone
from pathlib import Path
import runtime as rt
import v5a_loco_baseline as parent
import v5b_loco_baseline as baseline


def main():
    assert not baseline.SPEC.exists(), 'Never replace a frozen contract'
    original = parent.read_spec()
    diagnostic = rt.json.loads(Path('artifacts/v5a_mse_diagnostic_round32.json').read_text())
    plan = rt.json.loads(Path('artifacts/v5a_mse_diagnostic_plan_round32.json').read_text())
    assert diagnostic['status'].startswith('PASS same-count-stream')
    assert diagnostic['auditor_sha256'] == plan['auditor_sha256'] == rt.sha('audit_v5a_mse.py')
    assert diagnostic['plan_sha256'] == rt.sha('artifacts/v5a_mse_diagnostic_plan_round32.json')
    assert diagnostic['completed_score_audit_sha256'] == rt.sha('artifacts/v5a_completed_audit_round32.json')
    for path, digest in plan['preserve_sha256'].items():
        assert rt.sha(path) == digest, path
    events = [rt.json.loads(x) for x in Path('ledger/steps.jsonl').read_text().splitlines()]
    assert any(e.get('name') == 'audit_v5a_mse_round32' and e.get('event') == 'end'
               and e['rc'] == 0 for e in events)
    build_path = baseline.SNAPSHOT/'artifacts/v5b_build.json'
    build = rt.json.loads(build_path.read_text())
    assert build['arm'] == 'b'
    for filename, key in [('v5.py', 'v5_sha256'), ('model.py', 'model_sha256'),
                          ('build_v5.py', 'builder_sha256')]:
        assert rt.sha(baseline.SNAPSHOT/filename) == build[key]
        assert rt.sha(baseline.SNAPSHOT/filename) == rt.sha(parent.SNAPSHOT/filename)
    assert rt.sha('vendor/v4_model.py') == build['model_sha256']
    assert all(c['calls_k'] == 18533 and c['deflation_k'] == 50 and c['rho'] == 1
               for c in build['contexts'])
    package = rt.json.loads((baseline.SNAPSHOT/'artifacts/v5b_package.json').read_text())
    package_path = Path('../model_v5/v5b.vcc')
    assert package_path.stat().st_size == package['package_bytes']
    assert rt.sha(package_path) == package['package_sha256']
    tracked = ['v5b_loco_baseline.py', 'freeze_v5b_loco.py', 'verify_v5b_loco.py',
               'run_v5b_loco.py', 'v5a_loco_baseline.py', 'verify_v5a_loco.py']
    tracked += [str(p) for p in baseline.SNAPSHOT.rglob('*')
                if p.is_file() and '__pycache__' not in p.parts]
    implementation = dict(original['implementation_sha256'])
    implementation.update({p: rt.sha(p) for p in tracked})
    protocol = dict(original['protocol_sha256'])
    for p in ['artifacts/v5a_loco_spec.json', 'artifacts/v5a_loco_audit.json',
              'artifacts/v5a_loco_results.json', 'artifacts/v5a_completed_audit_round32.json',
              'artifacts/v5a_mse_diagnostic_plan_round32.json',
              'artifacts/v5a_mse_diagnostic_round32.json']:
        protocol[p] = rt.sha(p)
    rt.save(baseline.SPEC, {
        'utc': datetime.now(timezone.utc).isoformat(), 'round': 33,
        'baseline_id': baseline.BASELINE_ID,
        'phase': '1 fixed supplied-recipe reproduction extension; no variant selection',
        'tasks': original['tasks'], 'implementation_sha256': implementation,
        'protocol_sha256': protocol, 'source_sha256': original['source_sha256'],
        'official_environment': original['official_environment'],
        'recipe': {**original['recipe'], 'arm': 'b', 'parent': parent.BASELINE_ID,
                   'callset_control': {'k': 18533, 'z_sig': 1.96, 'min_net': 1,
                      'ordering': 'Original net significant sign votes then same-sign zsum',
                      'scope': 'Uncapped bulk-response mask, not actual Wilcoxon call-set enforcement'},
                   'deflation': {'k': 50, 'eta': 1, 'max_cells': 4000,
                                 'seed': 'Fixed context index0/1',
                                 'input': 'Only provided input controls, log1p CP50K'},
                   'order': 'v4a bulk -> deflate -> source sign gate -> original own-claim amplitude -> unchanged emission'},
        'adaptations': original['adaptations'] + [
            'Supplied v5b original build uses k18533, not helper default275. Keep all source-supported sign votes.',
            'Original emit_worker ignores sharp; keep the original pooled-moment emitter, without shape control.',
            'Mask zeroes uncalled log-bulk deltas before probability renormalization; it does not guarantee zero DE on other genes.',
            'Sources and destination-family exclusions identical to parent port, including source-specific variance masks.',
            'Seeded approximate top50 PCA is fixed before any score; no k, eta, call-count or amplitude sweep.',
            'Supplied official0.1528 remains unpaired; neither parent nor extension can replace failed historical rank gate.'],
        'original_package_verified': {'path': str(package_path), 'bytes': package['package_bytes'],
                                      'sha256': package['package_sha256']},
        'external_official_label': None, 'predicted_official_score': None,
        'phase2_allowed': False, 'holdout_scoring_allowed': False})
    print('FROZEN', baseline.SPEC, rt.sha(baseline.SPEC), flush=True)


if __name__ == '__main__':
    main()

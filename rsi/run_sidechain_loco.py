"""Serial fixed Phase1 baseline, official scoring and exact saved arithmetic."""
from datetime import datetime, timezone
from pathlib import Path
import runtime as rt
import sidechain_loco_baseline as baseline


def main():
    spec = baseline.read_spec()
    gate_before = rt.sha('artifacts/calibration_gate.json')
    from audit_phase1_calibration import exact_scaled
    import provenance
    results = {}
    for task in baseline.TASKS:
        result = baseline.score(task)
        anchors = rt.json.loads(Path(f'artifacts/scores/{task}/anchors.json').read_text())
        assert result['reference_sha256'] == anchors['real_sha256']
        assert result['anchors_sha256'] == rt.sha(f'artifacts/scores/{task}/anchors.json')
        exact_scaled(result, anchors)
        signature = provenance.identity(task)
        for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
            old = rt.json.loads(Path(f'artifacts/scores/{task}/{recipe}_s42.json').read_text())
            assert old['implementation_sha256'] == signature
        results[task] = result
    assert rt.sha('artifacts/calibration_gate.json') == gate_before
    baseline.read_spec()
    record = {'utc': datetime.now(timezone.utc).isoformat(),
        'status': 'PASS completed available-moment Sidechain kernels and exact score arithmetic',
        'baseline_id': baseline.BASELINE_ID, 'spec_sha256': rt.sha(baseline.SPEC),
        'runner_sha256': rt.sha(__file__), 'tasks': results,
        'equal_context_macro': sum(r['overall'] for r in results.values())/len(results),
        'external_official_label': None, 'predicted_official_score': None,
        'historical_score_identities_unchanged': 10, 'gate_unchanged_sha256': gate_before,
        'limitations': spec['adaptations'], 'phase2_allowed': False, 'holdout_scored': False}
    rt.save('artifacts/sidechain_loco_results.json', record)
    print('PASS scoped Sidechain baseline', record['equal_context_macro'], flush=True)


if __name__ == '__main__':
    main()

"""Serial fixed Phase1 baseline scoring; exact saved arithmetic audit."""
from datetime import datetime, timezone
from pathlib import Path
import runtime as rt
import zsp_loco_baseline as baseline


def main():
    from audit_phase1_calibration import exact_scaled
    import provenance
    spec = baseline.read_spec()
    gate = rt.sha('artifacts/calibration_gate.json')
    results = {}
    for task in baseline.TASKS:
        result = baseline.score(task)
        anchors = rt.json.loads(Path(f'artifacts/scores/{task}/anchors.json').read_text())
        exact_scaled(result, anchors)
        identity = provenance.identity(task)
        for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
            prior = rt.json.loads(Path(f'artifacts/scores/{task}/{recipe}_s42.json').read_text())
            assert prior['implementation_sha256'] == identity
        results[task] = result
    assert rt.sha('artifacts/calibration_gate.json') == gate
    baseline.read_spec()
    record = {'utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS',
        'baseline_id': baseline.BASELINE_ID, 'spec_sha256': rt.sha(baseline.SPEC),
        'runner_sha256': rt.sha(__file__), 'tasks': results,
        'equal_context_macro': sum(r['overall'] for r in results.values())/len(results),
        'external_official_label': None, 'predicted_official_score': None,
        'gate_unchanged_sha256': gate, 'historical_score_identities_unchanged': 10,
        'limitations': spec['adaptations'], 'phase2_allowed': False, 'holdout_scored': False}
    rt.save('artifacts/zsp_loco_results.json', record)
    print('PASS fixed Gerami four-source baseline', record['equal_context_macro'], flush=True)


if __name__ == '__main__':
    main()

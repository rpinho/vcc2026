"""Serial Phase1 baseline scoring and saved-score arithmetic audit.

Run through rsi_run.py --detach. Requires the separate input/kernel audit PASS;
never invokes HOLDOUT, modifies the historical calibration, or opens the gate.
"""
from datetime import datetime, timezone
from pathlib import Path
import runtime as rt
import deflate50_loco_baseline as baseline


def main():
    spec = baseline.read_spec()
    gate = Path('artifacts/calibration_gate.json')
    gate_before = rt.sha(gate)
    results = {}
    # Reuse the separately written exact-rational audit, not production scale().
    from audit_phase1_calibration import exact_scaled
    import provenance
    for task in baseline.TASKS:
        result = baseline.score(task)
        anchors = rt.json.loads(Path(f'artifacts/scores/{task}/anchors.json').read_text())
        assert result['reference_sha256'] == anchors['real_sha256']
        assert result['anchors_sha256'] == rt.sha(f'artifacts/scores/{task}/anchors.json')
        exact_scaled(result, anchors)
        signature = provenance.identity(task)
        for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
            historical = rt.json.loads(Path(f'artifacts/scores/{task}/{recipe}_s42.json').read_text())
            assert historical['implementation_sha256'] == signature
        results[task] = result
    assert rt.sha(gate) == gate_before
    baseline.read_spec()
    record = {'utc': datetime.now(timezone.utc).isoformat(),
        'status': 'PASS completed two-context scope-limited baseline and exact score-arithmetic audit',
        'baseline_id': baseline.BASELINE_ID, 'spec_sha256': rt.sha(baseline.SPEC),
        'runner_sha256': rt.sha(__file__), 'tasks': results,
        'equal_context_macro': sum(r['overall'] for r in results.values())/len(results),
        'external_official_label': None, 'predicted_official_score': None,
        'gate_unchanged_sha256': gate_before, 'historical_score_identities_unchanged': 10,
        'limitations': spec['adaptations'], 'phase2_allowed': False, 'holdout_scored': False}
    rt.save('artifacts/deflate50_loco_results.json', record)
    print('PASS Phase1 deflate50 available-input baseline', record['equal_context_macro'], flush=True)


if __name__ == '__main__':
    main()

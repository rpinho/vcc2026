"""Score both frozen SEARCH tasks only after the independent adapter audit."""
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import math
import runtime as rt
import v5a_loco_baseline as baseline


def main():
    spec = baseline.read_spec()
    from audit_phase1_calibration import exact_scaled
    import provenance
    results = {}
    for task in baseline.TASKS:
        result = baseline.score(task)
        anchors = rt.json.loads(Path(f'artifacts/scores/{task}/anchors.json').read_text())
        assert result['reference_sha256'] == anchors['real_sha256']
        exact_scaled(result, anchors)
        current = provenance.identity(task)
        assert current == spec['tasks'][task]['historical_implementation_identity']
        for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
            old = rt.json.loads(Path(f'artifacts/scores/{task}/{recipe}_s42.json').read_text())
            assert old['implementation_sha256'] == current
        results[task] = result
    baseline.read_spec()
    macro = float(sum(Fraction(r['overall']) for r in results.values())/len(results))
    assert math.isfinite(macro)
    rt.save('artifacts/v5a_loco_results.json', {'utc': datetime.now(timezone.utc).isoformat(),
        'status': 'PASS completed fixed v5a LOCO port and exact score arithmetic',
        'baseline_id': baseline.BASELINE_ID, 'spec_sha256': rt.sha(baseline.SPEC),
        'runner_sha256': rt.sha(__file__), 'tasks': results, 'equal_context_macro': macro,
        'external_official_label': None, 'predicted_official_score': None,
        'historical_score_identities_unchanged': 10,
        'gate_unchanged_sha256': rt.sha('artifacts/calibration_gate.json'),
        'limitations': spec['adaptations'], 'phase2_allowed': False, 'holdout_scored': False})
    print('PASS scoped v5a baseline', macro, flush=True)


if __name__ == '__main__':
    main()

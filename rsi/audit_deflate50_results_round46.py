"""Audit completed fixed public baseline; no scoring, fitting or task changes."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import gc
import json
import math
from pathlib import Path

import runtime as rt

PLAN = Path('artifacts/deflate50_results_plan_round46.json')
OUT = Path('artifacts/deflate50_results_audit_round46.json')


def load(path):
    return json.loads(Path(path).read_text())


def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-13), (a, b)


def freeze():
    assert not PLAN.exists()
    spec = load('artifacts/deflate50_loco_spec.json')
    paths = ['audit_deflate50_results_round46.py', 'artifacts/deflate50_loco_spec.json',
             'artifacts/deflate50_loco_audit.json', 'artifacts/deflate50_loco_results.json',
             'artifacts/frozen_metadata_index.json', 'artifacts/frozen_metadata.tar.gz',
             'artifacts/download_allocations.json', 'audit_phase1_calibration.py',
             'artifacts/deflate50_primary_descriptions_round45.json']
    for group in ['implementation_sha256', 'protocol_sha256', 'source_files']:
        paths.extend(spec[group])
    primary = load('artifacts/deflate50_primary_descriptions_round45.json')
    paths.extend([primary['source'], 'artifacts/bounded_score_fixture.json'])
    for task in spec['tasks']:
        folder = Path('artifacts/baseline_scores')/spec['baseline_id']/task
        paths.extend([str(folder/'s20260910.json'), str(folder/'tidy.csv')])
    rt.save(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        scope='Completed fixed Phase1 baseline verification only; no new predictions or selection',
        preserve_sha256={p: rt.sha(p) for p in sorted(set(paths))}))
    print('FROZEN', PLAN, flush=True)


def main():
    import benchmark
    import official_scoring as scoring
    import polars as pl
    import provenance
    import deflate50_loco_baseline as baseline
    from audit_phase1_calibration import exact_scaled

    assert not OUT.exists(), 'Preserve audit evidence'
    plan = load(PLAN)
    for path, digest in plan['preserve_sha256'].items():
        assert rt.sha(path) == digest, path
    spec = baseline.read_spec()
    audit, result = load(baseline.AUDIT), load('artifacts/deflate50_loco_results.json')
    assert result['status'] == 'PASS completed two-context scope-limited baseline and exact score-arithmetic audit'
    assert audit['status'] == 'PASS'
    assert result['baseline_id'] == baseline.BASELINE_ID
    assert result['spec_sha256'] == audit['spec_sha256'] == rt.sha(baseline.SPEC)
    assert result['runner_sha256'] == rt.sha('run_deflate50_loco.py')
    assert audit['verifier_sha256'] == rt.sha('verify_deflate50_loco.py')
    primary = load('artifacts/deflate50_primary_descriptions_round45.json')
    assert primary['source_sha256'] == rt.sha(primary['source'])
    assert len(primary['author_descriptions']) == 2
    for line in primary['author_descriptions']:
        assert line in Path(primary['source']).read_text().splitlines()
    assert scoring.provenance() == spec['official_environment']
    fixture = load('artifacts/bounded_score_fixture.json')
    assert fixture['status'] == 'PASS bitwise equality'
    assert fixture['driver_sha256'] == rt.sha('bounded_score.py')
    assert fixture['official'] == spec['official_environment']
    steps = [json.loads(line) for line in Path('ledger/steps.jsonl').read_text().splitlines()]
    step_ends = {}
    for name in ['verify_deflate50_round45', 'score_deflate50_round45']:
        ends = [s for s in steps if s.get('name') == name and s['event'] == 'end']
        assert ends and ends[-1]['rc'] == 0
        step_ends[name] = ends[-1]
    frozen = load('artifacts/frozen_metadata_index.json')
    for path, entry in frozen['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    assert rt.sha(frozen['archive']) == frozen['archive_sha256']
    assert set(result['tasks']) == set(audit['tasks']) == set(baseline.TASKS)
    gate = load('artifacts/calibration_gate.json')
    assert gate['status'] == 'CLOSED' and not gate['phase2_allowed']
    assert result['gate_unchanged_sha256'] == rt.sha('artifacts/calibration_gate.json')
    assert result['external_official_label'] is result['predicted_official_score'] is None
    assert not result['holdout_scored'] and not result['phase2_allowed']
    records, exact_overalls = {}, []
    for task in baseline.TASKS:
        print('Auditing saved result and current SEARCH reference', task, flush=True)
        folder = Path('artifacts/baseline_scores')/baseline.BASELINE_ID/task
        saved = load(folder/'s20260910.json')
        assert saved == result['tasks'][task]
        assert saved['task'] == task and saved['seed'] == spec['seed'] == spec['axis_seed'] == 20260910
        assert saved['spec_sha256'] == rt.sha(baseline.SPEC)
        assert saved['audit_sha256'] == rt.sha(baseline.AUDIT)
        anchors = load(f'artifacts/scores/{task}/anchors.json')
        assert saved['anchors_sha256'] == rt.sha(f'artifacts/scores/{task}/anchors.json')
        assert saved['reference_sha256'] == anchors['real_sha256']
        contract, ctrl, real = benchmark.setup(task)
        assert contract['split'] == 'search'
        assert list(contract['targets']) == spec['tasks'][task]['targets']
        assert scoring.matrix_id(real) == saved['reference_sha256']
        assert anchors['n_splits'] == 5
        assert set(saved['scaled']) == set(scoring.SIX)
        for metric in scoring.SIX:
            a, b = anchors['baseline'][metric], anchors['replicate'][metric]
            assert (b < a) if metric in ['expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae'] else (b > a)
        exact_overalls.append(exact_scaled(saved, anchors))
        raw = benchmark.aggregate(pl.read_csv(folder/'tidy.csv'))
        assert set(raw) == set(saved['raw'])
        for metric, value in raw.items():
            close(value, saved['raw'][metric])
        identity = provenance.identity(task)
        for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
            prior = load(f'artifacts/scores/{task}/{recipe}_s42.json')
            assert prior['implementation_sha256'] == identity
            assert prior['reference_sha256'] == saved['reference_sha256']
            exact_scaled(prior, anchors)
        records[task] = dict(overall=saved['overall'], scaled=saved['scaled'],
            raw_aggregation_checked=True, reference_verified=True,
            historical_score_identities_and_arithmetic_checked=5,
            deflation_diagnostics=saved['diagnostics'])
        del real, ctrl
        gc.collect()
    macro = sum(exact_overalls, Fraction())/len(exact_overalls)
    close(macro, result['equal_context_macro'])
    for path, digest in plan['preserve_sha256'].items():
        assert rt.sha(path) == digest, path
    rt.save(OUT, dict(utc=datetime.now(timezone.utc).isoformat(), status='PASS',
        plan_sha256=rt.sha(PLAN), verifier_sha256=rt.sha(__file__),
        spec_sha256=rt.sha(baseline.SPEC), result_sha256=rt.sha('artifacts/deflate50_loco_results.json'),
        step_ends=step_ends, tasks=records, equal_context_macro=float(macro),
        frozen_metadata_members_unchanged=len(frozen['entries']),
        historical_score_identities_unchanged=10, calibration_gate='CLOSED',
        phase2_allowed=False, holdout_scored=False, new_predictions_scored=False,
        external_official_label=None, predicted_official_score=None,
        scope='Saved-score arithmetic, tidy aggregation, current reference and immutable provenance audit. Prediction stream not regenerated; prior source/head/RNG audit retained. Description-based available-input deflate50, not exact Becker c6d81cc release or its external 0.1669 score.'))
    print('PASS', OUT, float(macro), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    freeze() if args.freeze else main()

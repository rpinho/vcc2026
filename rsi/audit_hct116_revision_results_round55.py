"""Independent saved-result audit for the fixed HCT116 source revision.

Freeze completed output identities with --freeze, then run without arguments.
Never regenerates predictions, opens HOLDOUT, or changes the historical gate.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
import gc
import math

import runtime as rt

ROOT = Path('artifacts/HCT116_revision_round54')
PLAN = ROOT/'results_audit_plan_round55.json'
OUT = ROOT/'results_audit_round55.json'
TASKS = ['H1_rich_search', 'CD4_depth192_search']
RECIPES = ['v1', 'v2b', 'v3b', 'v4a', 'v4b']
OFFICIAL = dict(zip(RECIPES, [.0877, .0617, -.2482, .1438, .1243]))


def load(path):
    return rt.json.loads(Path(path).read_text())


def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-13), (a, b)


def correlation(values):
    # Independent scalar ranks/moments; no call to the runner's correlation helper.
    xx = [float(values[r]) for r in RECIPES if r in values]
    yy = [OFFICIAL[r] for r in RECIPES if r in values]
    def ranks(x):
        return [1 + sum(v < a for v in x) + (sum(v == a for v in x)-1)/2 for a in x]
    def pearson(x, y):
        mx, my = math.fsum(x)/len(x), math.fsum(y)/len(y)
        a, b = [v-mx for v in x], [v-my for v in y]
        return math.fsum(u*v for u, v in zip(a, b))/math.sqrt(
            math.fsum(u*u for u in a)*math.fsum(v*v for v in b))
    rx, ry = ranks(xx), ranks(yy)
    rho = pearson(rx, ry)
    null = list(permutations(ry))
    p = sum(abs(pearson(rx, z)) >= abs(rho)-1e-12 for z in null)/len(null)
    return dict(n=len(xx), spearman=rho, pearson=pearson(xx, yy),
                exact_two_sided_spearman_p=p)


def successful_steps():
    rows = [rt.json.loads(x) for x in Path('ledger/steps.jsonl').read_text().splitlines()]
    result = {}
    for name in ['freeze_HCT116_revision_round54', 'verify_HCT116_revision_round54',
                 'recalibrate_HCT116_revision_round54']:
        events = [r for r in rows if r.get('name') == name]
        assert events and events[-1]['event'] == 'end' and events[-1]['rc'] == 0, name
        result[name] = events[-1]
    return result


def freeze():
    assert not PLAN.exists(), 'Preserve frozen audit plan'
    successful_steps()
    result = load(ROOT/'results.json')
    assert set(result['tasks']) == set(TASKS)
    paths = {str(ROOT/p) for p in ['plan.json', 'audit.json', 'results.json']}
    paths.update(['audit_hct116_revision_results_round55.py',
                  'artifacts/phase1_calibration_audit.json',
                  'artifacts/bounded_score_fixture.json'])
    paths.update(load(ROOT/'plan.json')['sha256'])
    for task in TASKS:
        for recipe in RECIPES:
            paths.add(result['tasks'][task][recipe]['path'])
        for recipe in ['v4a', 'v4b']:
            paths.add(str(ROOT/'scores'/task/(recipe+'_tidy.csv')))
    rt.save(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        scope='Completed fixed-recipe saved scores only; no new selection or predictions',
        sha256={p: rt.sha(p) for p in sorted(paths)}))
    print('FROZEN', PLAN, flush=True)


def main():
    import numpy as np
    import pandas as pd
    import polars as pl
    import benchmark
    import official_scoring as scoring
    import provenance
    import hct116_revision_round54 as revision
    from audit_phase1_calibration import exact_scaled

    assert not OUT.exists(), 'Preserve independent audit evidence'
    audit_plan = load(PLAN)
    for path, digest in audit_plan['sha256'].items():
        assert rt.sha(path) == digest, path
    steps = successful_steps()
    plan = revision.check_plan()
    result, prerequisite = load(ROOT/'results.json'), load(ROOT/'audit.json')
    assert result['plan_sha256'] == prerequisite['plan_sha256'] == rt.sha(revision.PLAN)
    assert result['runner_sha256'] == rt.sha('run_hct116_revision_round54.py')
    assert prerequisite['status'].startswith('PASS')
    assert prerequisite['verifier_sha256'] == rt.sha('verify_hct116_revision_round54.py')
    assert scoring.provenance() == plan['official_environment']
    fixture = load('artifacts/bounded_score_fixture.json')
    assert fixture['status'] == 'PASS bitwise equality'
    assert fixture['official'] == plan['official_environment']
    assert fixture['driver_sha256'] == rt.sha('bounded_score.py')
    for recipe, value in OFFICIAL.items():
        line, = [s for s in Path('../OFFICIAL_SCORES.md').read_text().splitlines()
                 if s.startswith('| '+recipe+' ')]
        close(float(line.split('|')[2].strip().replace('**', '')), value)
    roles = load('artifacts/benchmark_roles.json')
    assert roles['primary_search'] == TASKS
    policy = load('artifacts/anchor_policy.json')
    tasks, overalls = {}, {}
    for task in TASKS:
        print('Auditing fixed revised SEARCH', task, flush=True)
        contract, ctrl, real = benchmark.setup(task)
        assert contract == plan['tasks'][task] and contract['split'] == 'search'
        assert contract['exact_matched'] and len(contract['targets']) == 79
        original = Path('artifacts/scores')/task
        anchors = load(original/'anchors.json')
        reference = scoring.matrix_id(real)
        identity = provenance.identity(task)
        assert identity == plan['original_implementation_sha256'][task]
        assert anchors['real_sha256'] == reference
        assert anchors['n_splits'] == policy.get(task, policy['default_splits']) == 5
        assert anchors['version'] == '0.18.0'
        assert anchors['driver_sha256'] == rt.sha('bounded_score.py')
        assert anchors['scorer_sha256'] == rt.sha('official_scoring.py')
        assert anchors['baseline'] == load(original/'baseline_mean_response.json')
        assert anchors['baseline_construction']['reference_sha256'] == reference
        splits = pd.read_csv(original/'anchor_splits.csv', float_precision='round_trip')
        summary = pd.read_csv(original/'anchor_summary.csv',
                              float_precision='round_trip').set_index('metric')
        seeds = np.random.SeedSequence(0).generate_state(5).astype(np.int64)
        assert set(scoring.SIX) <= set(splits.metric)
        for metric in scoring.SIX:
            rows = splits.loc[splits.metric == metric].sort_values('split_index')
            assert rows.split_index.tolist() == list(range(5))
            np.testing.assert_array_equal(rows.seed, seeds)
            values = rows.value.to_numpy()
            for key, value in [('replicate', values.mean()),
                               ('replicate_sd', values.std(ddof=0)),
                               ('replicate_min', values.min()),
                               ('replicate_max', values.max())]:
                close(summary.loc[metric, key], value)
            close(anchors['replicate'][metric], values.mean())
            delta = F(anchors['replicate'][metric])-F(anchors['baseline'][metric])
            lower = metric in ['expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae']
            assert delta < 0 if lower else delta > 0
        overalls[task], records = {}, {}
        assert set(result['tasks'][task]) == set(RECIPES)
        for recipe in RECIPES:
            previous = load(original/(recipe+'_s42.json'))
            assert previous['implementation_sha256'] == identity
            assert previous['reference_sha256'] == reference
            exact_scaled(previous, anchors)
            item = result['tasks'][task][recipe]
            path = (original if recipe in RECIPES[:3] else ROOT/'scores'/task)/(recipe+'_s42.json')
            assert item['path'] == str(path) and item['sha256'] == rt.sha(path)
            record = load(path)
            assert record['task'] == task and record['recipe'] == recipe and record['seed'] == 42
            assert record['reference_sha256'] == reference
            assert record['selection_sha256'] == contract['selection_sha256']
            assert record['cells_pred'] == 400 and record['n_targets'] == 79
            assert set(scoring.SIX) <= set(record['raw'])
            assert set(record['scaled']) == set(scoring.SIX)
            assert len(bytes.fromhex(record['prediction_sha256'])) == 32
            exact_scaled(record, anchors)
            close(item['overall'], record['overall'])
            legacy = recipe in RECIPES[:3]
            assert item['reused_unchanged_legacy_result'] is legacy
            if not legacy:
                assert record['effective_seed'] == 20260910
                assert record['revision_id'] == plan['revision_id']
                assert record['plan_sha256'] == rt.sha(revision.PLAN)
                assert record['audit_sha256'] == rt.sha(ROOT/'audit.json')
                assert record['anchors_sha256'] == rt.sha(original/'anchors.json')
                tidy_path = path.with_name(recipe+'_tidy.csv')
                raw = benchmark.aggregate(pl.read_csv(tidy_path))
                assert set(raw) == set(record['raw'])
                for metric in raw:
                    close(raw[metric], record['raw'][metric])
            overalls[task][recipe] = record['overall']
            records[recipe] = dict(overall=record['overall'],
                sha256=rt.sha(path), unchanged_legacy=legacy,
                scaled=record['scaled'], saved_aggregation_checked=not legacy)
        tasks[task] = dict(reference_sha256=reference,
            original_implementation_sha256=identity, records=records,
            correlation=correlation(overalls[task]), six_scales_valid=True)
        del real, ctrl
        gc.collect()
    macro = {r: float(sum(F(overalls[t][r]) for t in TASKS)/2) for r in RECIPES}
    for recipe in RECIPES:
        close(macro[recipe], result['macro'][recipe])
    pairs = [dict(a=a, b=b, matches=(macro[a]-macro[b])*(OFFICIAL[a]-OFFICIAL[b]) > 0)
             for a, b in combinations(RECIPES, 2)]
    assert pairs == result['all_pairs']
    assert result['all_pairs_rank_consistent'] is all(p['matches'] for p in pairs)
    corr, without = correlation(macro), correlation({r: v for r, v in macro.items() if r != 'v3b'})
    for field, expected in [('correlation', corr), ('without_v3b_diagnostic_only', without)]:
        for key, value in expected.items():
            close(value, result[field][key])
    assert result['historical_gate_preserved']
    assert not result['phase2_allowed'] and not result['holdout_scored']
    gate = load('artifacts/calibration_gate.json')
    assert gate['status'] == 'CLOSED' and not gate['phase2_allowed']
    manifest = load('artifacts/tasks/manifest.json')
    old = load('artifacts/protocol_revisions/CD4_depth192/manifest.json.before')
    for task, contract in manifest.items():
        if contract['split'] == 'holdout':
            assert contract == old[task]
    assert not Path('artifacts/holdout_events.jsonl').exists()
    assert not list(ROOT.glob('scores/*holdout*/*'))
    revision.check_plan()
    for path, digest in audit_plan['sha256'].items():
        assert rt.sha(path) == digest, path
    rt.save(OUT, dict(status='PASS saved identities, scales, aggregation and independent rank arithmetic',
        utc=datetime.now(timezone.utc).isoformat(), verifier_sha256=rt.sha(__file__),
        plan_sha256=rt.sha(PLAN), source_revision_plan_sha256=rt.sha(revision.PLAN),
        results_sha256=rt.sha(ROOT/'results.json'), steps=steps, tasks=tasks,
        macro=macro, correlation=corr, all_pairs=pairs,
        without_v3b_diagnostic_only=without,
        original_scores_unchanged=10, original_metadata_members_unchanged=38,
        historical_gate='CLOSED', holdout_scored=False, phase2_allowed=False,
        scope='Saved-result verification; prediction stream not replayed. Revised-source ports are not exact submitted-score pairings. No official forecast or external-baseline reproduction claim.'))
    print('PASS saved revised results', macro, corr, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    freeze() if args.freeze else main()

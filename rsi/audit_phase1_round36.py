"""Certify the completed fixed v5b port; diagnose saved-anchor precision only.

No expression matrices, new predictions, or HOLDOUT files are opened. The five
delete-one means per context are all reported, never chosen as a replacement.
This is a finite sensitivity check, not a bootstrap or confidence interval.
"""
import argparse
import csv
import itertools
import json
import math
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

import runtime as rt

PLAN = Path('artifacts/phase1_round36_plan.json')
OUT = Path('artifacts/phase1_round36_audit.json')
TASKS = ['H1_rich_search', 'CD4_depth192_search']
RECIPES = ['v1', 'v2b', 'v3b', 'v4a', 'v4b']


def load(path):
    return json.loads(Path(path).read_text())


def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-13), (a, b)


def freeze():
    assert not PLAN.exists(), 'Do not overwrite a frozen plan'
    paths = ['audit_phase1_round36.py', 'artifacts/v5b_loco_spec.json',
             'artifacts/v5b_loco_audit.json', 'artifacts/v5b_loco_results.json',
             'artifacts/frozen_metadata_index.json', 'artifacts/frozen_metadata.tar.gz',
             'artifacts/download_allocations.json']
    spec = load('artifacts/v5b_loco_spec.json')
    for group in ['implementation_sha256', 'protocol_sha256', 'source_sha256']:
        paths.extend(spec[group])
    for task in TASKS:
        paths.extend(f'artifacts/scores/{task}/{name}' for name in
                     ['anchor_splits.csv', 'anchor_summary.csv', 'anchors.json'])
        paths.append(f'artifacts/baseline_scores/v5b_available_loco_r33/{task}/s20260910.json')
    rt.save(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        status='FROZEN before anchor sensitivity computation',
        tasks=TASKS, recipes=RECIPES,
        method='All five delete-one split means per task; all25 joint combinations. Fixed raw scores, baselines, six metrics, official clips and equal context weights. Invalid scale excludes whole diagnostic combination.',
        interpretation='Descriptive sensitivity only: overlapping splits are dependent. No model-seed, cell-sampling or hidden-context uncertainty estimate. No replacement split choice or gate update.',
        next_reproduction='Competitor priority2: clean-room synthetic official-Wilcoxon call-set-control fixture, input controls and source-vote ranking only; no biological variant scoring before historical gate opens.',
        preserve_sha256={p: rt.sha(p) for p in sorted(set(paths))},
        holdout_scored=False, phase2_allowed=False))
    print('FROZEN', PLAN, flush=True)


def independent_scale(raw, base, rep, metrics):
    """Exact rational arithmetic on the binary64 scorer inputs."""
    scaled = {}
    for metric in metrics:
        denominator = F(rep[metric])-F(base[metric])
        positive = metric not in ['expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae']
        if not (denominator > 0 if positive else denominator < 0):
            return None
        value = (F(raw[metric])-F(base[metric]))/denominator
        if metric == 'expr_mse_unbiased_capped_norm':
            value = min(F(1), max(F(0), value))
        elif metric == 'de_wilcoxon_lfc_nmae':
            value = max(F(-6), value)
        scaled[metric] = value
    return scaled, sum(scaled.values(), F(0))/len(metrics)


def main():
    from audit_phase1_calibration import exact_scaled, OFFICIAL
    import official_scoring as scoring
    import provenance
    import v5b_loco_baseline as baseline
    import numpy as np

    plan = load(PLAN)
    for path, digest in plan['preserve_sha256'].items():
        assert rt.sha(path) == digest, path
    spec = baseline.read_spec()
    result = load('artifacts/v5b_loco_results.json')
    audit = load(baseline.AUDIT)
    assert result['status'].startswith('PASS completed')
    assert result['spec_sha256'] == audit['spec_sha256'] == rt.sha(baseline.SPEC)
    assert result['runner_sha256'] == rt.sha('run_v5b_loco.py')
    assert audit['status'] == 'PASS' and audit['verifier_sha256'] == rt.sha('verify_v5b_loco.py')
    assert set(result['tasks']) == set(audit['tasks']) == set(TASKS)
    assert not result['holdout_scored'] and not result['phase2_allowed']
    assert result['external_official_label'] is result['predicted_official_score'] is None
    assert result['gate_unchanged_sha256'] == rt.sha('artifacts/calibration_gate.json')
    ends = [json.loads(line) for line in Path('ledger/steps.jsonl').read_text().splitlines()]
    for name in ['verify_v5b_loco_round34', 'score_v5b_loco_after_r34']:
        event, = [e for e in ends if e.get('name') == name and e.get('event') == 'end']
        assert event['rc'] == 0
    index = load('artifacts/frozen_metadata_index.json')
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    assert rt.sha(index['archive']) == index['archive_sha256']
    history = load('artifacts/v5a_completed_audit_round32.json')['historical_sha256']
    for path, digest in history.items():
        assert rt.sha(path) == digest, path
    assert len(history) == 10
    assert load('artifacts/benchmark_roles.json')['primary_search'] == TASKS
    certificates, sensitivity = {}, {}
    for task in TASKS:
        folder = Path('artifacts/scores')/task
        anchors = load(folder/'anchors.json')
        record = result['tasks'][task]
        assert record == load(f'artifacts/baseline_scores/{baseline.BASELINE_ID}/{task}/s20260910.json')
        assert record['spec_sha256'] == rt.sha(baseline.SPEC)
        assert record['audit_sha256'] == rt.sha(baseline.AUDIT)
        assert record['anchors_sha256'] == rt.sha(folder/'anchors.json')
        assert record['reference_sha256'] == anchors['real_sha256']
        exact_scaled(record, anchors)
        identity = provenance.identity(task)
        assert identity == spec['tasks'][task]['historical_implementation_identity']
        old = {r: load(folder/f'{r}_s42.json') for r in RECIPES}
        for rec in old.values():
            assert rec['implementation_sha256'] == identity
            exact_scaled(rec, anchors)
        diag = record['diagnostics']['prediction_only_B_selfclaim']
        close(F(diag['B'])/F(diag['full_prediction_jackknife_claim']), diag['ratio'])
        certificates[task] = dict(raw=record['raw'], scaled=record['scaled'],
            overall=record['overall'], prediction_sha256=record['prediction_sha256'],
            alpha=record['diagnostics']['alpha'],
            probe_claim=record['diagnostics']['claim_est_probe'], final_emitted=diag)
        rows = list(csv.DictReader((folder/'anchor_splits.csv').open()))
        by_metric = {m: sorted([r for r in rows if r['metric'] == m],
                              key=lambda r: int(r['split_index'])) for m in scoring.SIX}
        seeds = np.random.SeedSequence(0).generate_state(5).astype(np.int64).tolist()
        for m, rr in by_metric.items():
            assert [int(r['split_index']) for r in rr] == list(range(5))
            assert [int(r['seed']) for r in rr] == seeds
            close(sum(F(float(r['value'])) for r in rr)/5, anchors['replicate'][m])
        scenarios = []
        for dropped in range(5):
            rep = {m: float(sum(F(float(r['value'])) for i, r in enumerate(rr)
                                if i != dropped)/4) for m, rr in by_metric.items()}
            headroom = {m: rep[m]-anchors['baseline'][m] for m in scoring.SIX}
            values = {}
            for recipe, rec in old.items():
                independent = independent_scale(rec['raw'], anchors['baseline'], rep, scoring.SIX)
                if independent is None:
                    break
                official_scaled, official_overall = scoring.scale(rec['raw'], anchors['baseline'], rep)
                for m in scoring.SIX:
                    close(independent[0][m], official_scaled[m])
                close(independent[1], official_overall)
                values[recipe] = official_overall
            scenarios.append(dict(dropped_split=dropped, headrooms=headroom,
                valid_all_six_scales=len(values) == len(RECIPES), local=values))
        sensitivity[task] = scenarios
    close(sum(F(r['overall']) for r in result['tasks'].values())/2, result['equal_context_macro'])
    joint = []
    for a, b in itertools.product(*[sensitivity[t] for t in TASKS]):
        valid = a['valid_all_six_scales'] and b['valid_all_six_scales']
        item = dict(dropped_splits=[a['dropped_split'], b['dropped_split']], valid=valid)
        if valid:
            macro = {r: float((F(a['local'][r])+F(b['local'][r]))/2) for r in RECIPES}
            mismatches = [[r, s] for r, s in itertools.combinations(RECIPES, 2)
                          if (macro[r]-macro[s])*(OFFICIAL[r]-OFFICIAL[s]) <= 0]
            item.update(local_macro=macro, mismatched_pairs=mismatches,
                        v4a_minus_v4b=macro['v4a']-macro['v4b'])
        joint.append(item)
    valid = [j for j in joint if j['valid']]
    summary = dict(combinations=len(joint), valid=len(valid), invalid=len(joint)-len(valid),
        all_pairs_consistent=sum(not j['mismatched_pairs'] for j in valid),
        v4a_v4b_reversed=sum(j['v4a_minus_v4b'] < 0 for j in valid),
        v4a_minus_v4b_range=([min(j['v4a_minus_v4b'] for j in valid),
                                  max(j['v4a_minus_v4b'] for j in valid)] if valid else None))
    for path, digest in plan['preserve_sha256'].items():
        assert rt.sha(path) == digest, path
    assert not Path('artifacts/holdout_events.jsonl').exists()
    rt.save(OUT, dict(status='PASS saved v5b provenance and exact score arithmetic; completed anchor sensitivity',
        utc=datetime.now(timezone.utc).isoformat(), auditor_sha256=rt.sha(__file__),
        plan_sha256=rt.sha(PLAN), result_sha256=rt.sha('artifacts/v5b_loco_results.json'),
        spec_bindings=sum(len(spec[g]) for g in ['implementation_sha256','protocol_sha256','source_sha256']),
        historical_identities=10, historical_sha256=history, frozen_metadata_members=len(index['entries']),
        tasks=certificates, equal_context_macro=result['equal_context_macro'],
        anchor_delete_one=sensitivity, joint_delete_one=joint, sensitivity_summary=summary,
        limitations=plan['interpretation'], models_changed=False, protocols_changed=False,
        gate_changed=False, holdout_scored=False, phase2_allowed=False))
    print('PASS v5b', result['equal_context_macro'], 'anchor sensitivity', summary, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    freeze() if args.freeze else main()

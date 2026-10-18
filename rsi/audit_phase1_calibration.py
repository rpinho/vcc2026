"""Certify the saved two-context calibration without selecting new data or splits.

Reads SEARCH matrices only, recomputes no biological metric, and never changes
the benchmark, gate, models, or official-score forecast. Split dispersion and
leave-v3b-out correlations are descriptive diagnostics, never selection rules.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
import gc
import math

import runtime as rt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

import benchmark
import official_scoring as scoring
import provenance


RECIPES = ['v1', 'v2b', 'v3b', 'v4a', 'v4b']
TASKS = ['H1_rich_search', 'CD4_depth192_search']
OFFICIAL = dict(zip(RECIPES, [.0877, .0617, -.2482, .1438, .1243]))
SIX = scoring.SIX


def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-13), (a, b)


def correlations(values):
    keys = [k for k in RECIPES if k in values]
    xx = [values[k] for k in keys]
    yy = [OFFICIAL[k] for k in keys]
    rho = float(spearmanr(xx, yy).statistic)
    null = list(permutations(yy))
    p = sum(abs(float(spearmanr(xx, y).statistic)) >= abs(rho)-1e-12
            for y in null) / len(null)
    return dict(n=len(keys), spearman=rho, pearson=float(pearsonr(xx, yy).statistic),
                exact_two_sided_spearman_p=p)


def exact_scaled(record, anchors):
    values = {}
    for metric in SIX:
        raw, base, rep = [F(float(v)) for v in
                          [record['raw'][metric], anchors['baseline'][metric],
                           anchors['replicate'][metric]]]
        value = (raw-base)/(rep-base)
        if metric == 'expr_mse_unbiased_capped_norm':
            value = min(F(1), max(F(0), value))
        if metric == 'de_wilcoxon_lfc_nmae':
            value = max(F(-6), value)
        close(value, record['scaled'][metric])
        values[metric] = value
    overall = sum(values.values(), F(0))/6
    close(overall, record['overall'])
    return overall


def main():
    roles_path = Path('artifacts/benchmark_roles.json')
    manifest_path = Path('artifacts/tasks/manifest.json')
    roles = rt.json.loads(roles_path.read_text())
    manifest = rt.json.loads(manifest_path.read_text())
    policy = rt.json.loads(Path('artifacts/anchor_policy.json').read_text())
    assert roles['primary_search'] == TASKS
    official_path = Path('../OFFICIAL_SCORES.md')
    # Independently bind the expected labels to the source table.
    for recipe, value in OFFICIAL.items():
        line, = [line for line in official_path.read_text().splitlines()
                 if line.startswith('| '+recipe+' ')]
        close(float(line.split('|')[2].strip().replace('**', '')), value)
    tracked = [roles_path, manifest_path, official_path,
               Path('artifacts/anchor_policy.json'), Path('artifacts/calibration.json'),
               Path('artifacts/calibration_gate.json')]
    original = {str(p): rt.sha(p) for p in tracked}
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    assert rt.sha(index['archive']) == index['archive_sha256']
    environment = scoring.provenance()
    fixture = rt.json.loads(Path('artifacts/bounded_score_fixture.json').read_text())
    assert fixture['status'] == 'PASS bitwise equality'
    assert fixture['driver_sha256'] == rt.sha('bounded_score.py')
    assert fixture['official'] == environment
    records, tasks = {}, {}
    for task in TASKS:
        print('Auditing SEARCH', task, flush=True)
        assert manifest[task]['split'] == 'search'
        contract, ctrl, real = benchmark.setup(task)
        reference = scoring.matrix_id(real)
        genes = set(real.var_names)
        unresolved = sorted(set(contract['targets'])-genes)
        implementation = provenance.identity(task)
        folder = Path('artifacts/scores')/task
        anchors = rt.json.loads((folder/'anchors.json').read_text())
        baseline = rt.json.loads((folder/'baseline_mean_response.json').read_text())
        assert anchors['real_sha256'] == reference
        assert anchors['scorer_sha256'] == rt.sha('official_scoring.py')
        assert anchors['driver_sha256'] == rt.sha('bounded_score.py')
        assert anchors['version'] == environment['version'] == '0.18.0'
        assert anchors['n_splits'] == policy.get(task, policy['default_splits']) == 5
        if 'replicate_splits' in contract:
            assert contract['replicate_splits'] == 5
        assert anchors['baseline'] == baseline
        protocol = anchors['baseline_construction']
        assert protocol['reference_sha256'] == reference
        if protocol['builder_sha256'] != rt.sha('baseline_protocol.py'):
            cert_path = Path('artifacts/baseline_profile_fixture.json')
            cert = rt.json.loads(cert_path.read_text())
            assert protocol['equivalence_certificate_sha256'] == rt.sha(cert_path)
            assert cert['status'] == 'PASS bitwise profile equality and unchanged reference'
            assert cert['builder_sha256'] == rt.sha('baseline_protocol.py')
            assert protocol['builder_sha256'] in cert['compatible_legacy_builder_sha256']
        assert not (folder/'invalid_scale.json').exists()
        splits = pd.read_csv(folder/'anchor_splits.csv', float_precision='round_trip')
        summary = pd.read_csv(folder/'anchor_summary.csv', float_precision='round_trip').set_index('metric')
        seeds = np.random.SeedSequence(0).generate_state(5).astype(np.int64)
        for metric, rows in splits.groupby('metric'):
            rows = rows.sort_values('split_index')
            assert rows.split_index.tolist() == list(range(5))
            assert np.array_equal(rows.seed, seeds)
            vv = rows.value.to_numpy()
            for key, value in [('replicate', vv.mean()), ('replicate_sd', vv.std(ddof=0)),
                               ('replicate_min', vv.min()), ('replicate_max', vv.max())]:
                close(summary.loc[metric, key], value)
            close(anchors['replicate'][metric], vv.mean())
        headrooms = {}
        for metric in SIX:
            base, rep = [float(anchors[k][metric]) for k in ['baseline', 'replicate']]
            delta = F(rep)-F(base)
            positive = metric not in ['expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae']
            assert delta > 0 if positive else delta < 0
            vv = splits.loc[splits.metric == metric, 'value'].to_numpy()
            headrooms[metric] = dict(baseline=base, replicate=rep,
                exact_binary64_headroom=str(delta), headroom=float(delta),
                split_population_sd=float(vv.std(ddof=0)),
                splits_in_usable_direction=int(np.sum(vv > base if positive else vv < base)),
                valid_frozen_mean_scale=True)
        records[task] = {}
        for recipe in RECIPES:
            path = folder/f'{recipe}_s42.json'
            record = rt.json.loads(path.read_text())
            assert (record['task'], record['recipe'], record['seed']) == (task, recipe, 42)
            assert record['implementation_sha256'] == implementation
            assert record['reference_sha256'] == reference
            assert record['selection_sha256'] == contract['selection_sha256']
            assert record['n_targets'] == len(contract['targets']) == 79
            assert record['cells_real'] == contract['n_real']
            assert record['cells_pred'] == 400
            assert record['exact_target_matched'] and record['primary_calibration_eligible']
            exact_scaled(record, anchors)
            records[task][recipe] = record
            tracked.append(path)
        tracked += [folder/name for name in ['anchors.json', 'anchor_splits.csv',
                    'anchor_summary.csv', 'baseline_mean_response.json',
                    'baseline_mean_response_protocol.json']]
        tasks[task] = dict(reference_sha256=reference, implementation_sha256=implementation,
            headrooms=headrooms, unresolved_labels_absent_from_measured_axis=unresolved,
            local={k: r['overall'] for k, r in records[task].items()},
            correlation=correlations({k: r['overall'] for k, r in records[task].items()}))
        del real, ctrl
        gc.collect()
    macro = {k: float(sum(F(records[t][k]['overall']) for t in TASKS)/len(TASKS)) for k in RECIPES}
    calibration = rt.json.loads(Path('artifacts/calibration.json').read_text())
    for k, value in macro.items():
        close(value, calibration['common_native_macro']['local'][k])
    corr = correlations(macro)
    close(corr['spearman'], calibration['common_native_macro']['spearman'])
    close(corr['pearson'], calibration['common_native_macro']['pearson'])
    pairs = []
    for a, b in combinations(RECIPES, 2):
        difference = macro[a]-macro[b]
        per_task = {}
        for task in TASKS:
            aa, bb = [records[task][k] for k in [a, b]]
            contributions = {metric: (aa['scaled'][metric]-bb['scaled'][metric])/6 for metric in SIX}
            close(sum(contributions.values()), aa['overall']-bb['overall'])
            per_task[task] = dict(overall_difference=aa['overall']-bb['overall'],
                                  metric_contributions=contributions)
        pairs.append(dict(a=a, b=b, official_difference=OFFICIAL[a]-OFFICIAL[b],
            macro_difference=difference, matches_official=difference*(OFFICIAL[a]-OFFICIAL[b]) > 0,
            per_task=per_task))
    mismatches = [p for p in pairs if not p['matches_official']]
    gate = rt.json.loads(Path('artifacts/calibration_gate.json').read_text())
    assert gate['status'] == 'CLOSED' and gate['phase2_allowed'] is False
    assert gate['macro']['all_pairs_rank_consistent'] == (not mismatches)
    # HOLDOUT contracts must equal the pre-CD4-depth192 revision. No HOLDOUT matrix is opened.
    old = rt.json.loads(Path('artifacts/protocol_revisions/CD4_depth192/manifest.json.before').read_text())
    for task, contract in manifest.items():
        if contract['split'] == 'holdout':
            assert contract == old[task]
    assert not Path('artifacts/holdout_events.jsonl').exists()
    assert not [p for p in Path('artifacts/scores').glob('*/*_s42.json') if 'holdout' in str(p).lower()]
    for p, digest in original.items():
        assert rt.sha(p) == digest
    rt.save('artifacts/phase1_calibration_audit.json', dict(
        status='PASS saved identities, official split summaries and exact score arithmetic; calibration FAILS rank gate',
        utc=datetime.now(timezone.utc).isoformat(), verifier_sha256=rt.sha(__file__),
        metadata_members_verified=len(index['entries']), tasks=tasks,
        macro=dict(local=macro, correlation=corr, all_pairs=pairs,
                   mismatched_pairs=[[p['a'], p['b']] for p in mismatches]),
        without_v3b_diagnostic_only=correlations({k: v for k, v in macro.items() if k != 'v3b'}),
        input_hashes={str(p): rt.sha(p) for p in tracked},
        holdout_scored=False, protocols_changed=False, phase2_allowed=False,
        interpretation='All six frozen mean scales are usable. CD4 fidelity headroom is small; its amplification of raw fidelity differences dominates the v4a/v4b reversal. This is descriptive attribution, not proof of a biological cause or permission to change splits, metric weights, targets or contexts. External exact AtlasShift/Sidechain reproduction remains unresolved; no official forecast is justified.'))
    print('PASS audit; rank gate CLOSED:', corr, 'mismatches', [(p['a'], p['b']) for p in mismatches], flush=True)


if __name__ == '__main__':
    main()

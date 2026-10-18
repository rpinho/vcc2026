"""Read-only Phase1 attribution of the completed, frozen v5a SEARCH result.

Replays the SAME count streams, never chooses an amplitude or changes a panel.
Independent per-coordinate arithmetic checks the official MSE numerator and
denominator. Reference-only regime diagnostics and a fixed-bulk correction bound
describe why a zero scaled score need not imply an incomplete correction refund.
"""
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import gc
import hashlib
import math

import runtime as rt
import numpy as np
import pandas as pd
import benchmark
import bounded_score
import official_scoring as scoring
import provenance
import stat_cache
import v5a_loco_baseline as baseline
from audit_phase1_calibration import exact_scaled
from cell_eval2.metrics.delta import (mse_unbiased_capped, distance_unbiased,
                                     PRED_TRACE_CAP_K, _across_pert_budget)

PLAN = Path('artifacts/v5a_mse_diagnostic_plan_round32.json')
OUT = Path('artifacts/v5a_mse_diagnostic_round32.json')
RECEIPTS = Path('artifacts/v5a_completed_audit_round32.json')


def load(path):
    return rt.json.loads(Path(path).read_text())


def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=2e-10, abs_tol=2e-13), (a, b)


def weighted_components(values, labels, genes):
    """Direct centred sums, with independent own-gene omission and row weights."""
    index = {g: j for j, g in enumerate(genes)}
    weights = np.array([1/(len(genes)-int(t in index)) for t in labels])
    centre_ss = spread_ss = 0.0
    for j, gene in enumerate(genes):
        keep = labels != gene
        w, v = weights[keep], values[keep, j]
        if not len(v):
            continue
        centre = float(np.dot(w, v)/sum(w))
        centre_ss += float(sum(w)*centre**2)
        if len(v) >= 2:
            spread_ss += float(np.dot(w, (v-centre)**2))
    return weights, centre_ss, spread_ss


def check_receipts(plan):
    spec = baseline.read_spec()
    result = load('artifacts/v5a_loco_results.json')
    audit = load(baseline.AUDIT)
    assert result['status'].startswith('PASS completed')
    assert result['spec_sha256'] == audit['spec_sha256'] == rt.sha(baseline.SPEC)
    assert result['runner_sha256'] == rt.sha('run_v5a_loco.py')
    assert audit['status'] == 'PASS' and audit['verifier_sha256'] == rt.sha('verify_v5a_loco.py')
    assert result['gate_unchanged_sha256'] == rt.sha('artifacts/calibration_gate.json')
    assert result['external_official_label'] is result['predicted_official_score'] is None
    assert not result['holdout_scored'] and not result['phase2_allowed']
    assert set(result['tasks']) == set(baseline.TASKS)
    ends = [rt.json.loads(line) for line in Path('ledger/steps.jsonl').read_text().splitlines()]
    assert any(e.get('name') == 'score_v5a_loco_round31' and e.get('event') == 'end'
               and e['rc'] == 0 for e in ends)
    index = load('artifacts/frozen_metadata_index.json')
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    assert rt.sha(index['archive']) == index['archive_sha256']
    historical, tasks = {}, {}
    for task in baseline.TASKS:
        anchors = load(f'artifacts/scores/{task}/anchors.json')
        record = result['tasks'][task]
        assert record['spec_sha256'] == rt.sha(baseline.SPEC)
        assert record['audit_sha256'] == rt.sha(baseline.AUDIT)
        assert record['anchors_sha256'] == rt.sha(f'artifacts/scores/{task}/anchors.json')
        assert record['reference_sha256'] == anchors['real_sha256']
        path = Path('artifacts/baseline_scores')/baseline.BASELINE_ID/task/'s20260910.json'
        assert load(path) == record
        exact_scaled(record, anchors)
        identity = provenance.identity(task)
        assert identity == spec['tasks'][task]['historical_implementation_identity']
        for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
            path = f'artifacts/scores/{task}/{recipe}_s42.json'
            old = load(path)
            assert old['implementation_sha256'] == identity
            exact_scaled(old, anchors)
            historical[path] = rt.sha(path)
        d = record['diagnostics']['prediction_only_B_selfclaim']
        close(d['B']/d['full_prediction_jackknife_claim'], d['ratio'])
        tasks[task] = {'raw': record['raw'], 'scaled': record['scaled'],
                       'overall': record['overall'], 'diagnostics': record['diagnostics']}
    close(sum(Fraction(r['overall']) for r in result['tasks'].values())/2,
          result['equal_context_macro'])
    for path, digest in plan['preserve_sha256'].items():
        assert rt.sha(path) == digest, path
    rt.save(RECEIPTS, {'status': 'PASS completed v5a receipts and exact six-score arithmetic',
        'utc': datetime.now(timezone.utc).isoformat(), 'auditor_sha256': rt.sha(__file__),
        'plan_sha256': rt.sha(PLAN), 'result_sha256': rt.sha('artifacts/v5a_loco_results.json'),
        'historical_sha256': historical, 'frozen_metadata_members': len(index['entries']),
        'tasks': tasks, 'equal_context_macro': result['equal_context_macro'],
        'phase2_allowed': False, 'holdout_scored': False})
    print('PASS completed v5a saved-score audit', result['equal_context_macro'], flush=True)
    return spec, result


def diagnose(task, spec, record):
    contract, ctrl, real = benchmark.setup(task)
    assert contract['split'] == 'search'
    assert scoring.matrix_id(real) == record['reference_sha256']
    cfg = scoring.config()
    genes = ctrl.var_names.to_numpy(str)
    rs, ps = {}, {}
    labels = real.obs.target_gene.astype(str).to_numpy()
    for target in sorted(set(labels)):
        rs[target] = bounded_score.stat(real[labels == target].copy(), cfg)
    # Replay unchanged construction; reference is never passed to the predictor.
    info, digest = {}, hashlib.sha256()
    for i, (target, block) in enumerate(baseline.predictions(ctrl, task, spec, info)):
        for value in [block.data, block.indices, block.indptr]:
            digest.update(value.tobytes())
        digest.update(target.encode())
        ps[target] = bounded_score.stat(bounded_score.make(block, target, genes), cfg)
        if i % 20 == 0:
            print(task, 'MSE-only replay', i+1, '/', len(contract['targets']), flush=True)
    assert digest.hexdigest() == record['prediction_sha256']
    for key in ['alpha', 'claim_est_probe', 'B_v4a_profile', 'B_target_profile']:
        close(info[key], record['diagnostics'][key])
    benchmark.release_models()
    pb, pm = bounded_score.assemble(ps)
    rb, rm = bounded_score.assemble(rs)
    targets, pred = pb['bulk_lognorm']
    ri = {t: i for i, t in enumerate(rb['bulk_lognorm'][0])}
    refs = rb['bulk_lognorm'][1][[ri[t] for t in targets]]
    refjk = rm['bulk_lognorm'].jk[[ri[t] for t in targets]]
    refctrl = rb['bulk_lognorm'][1][ri[bounded_score.NT]]
    ctrl_jk = float(rm['bulk_lognorm'].jk[ri[bounded_score.NT]])
    predjk = pm['bulk_lognorm'].jk
    weights, _, B = weighted_components(pred, targets, genes)
    _, ref_centre_ss, ref_B = weighted_components(refs, targets, genes)
    _, error_centre_ss, error_spread_ss = weighted_components(pred-refs, targets, genes)
    gene_index = {g: j for j, g in enumerate(genes)}
    excl = {i: gene_index[t] for i, t in enumerate(targets) if t in gene_index}
    close(B, _across_pert_budget(pred, excl, weights))
    close(ref_B, _across_pert_budget(refs, excl, weights))
    close(B, record['diagnostics']['prediction_only_B_selfclaim']['B'])
    close(float(weights @ predjk),
          record['diagnostics']['prediction_only_B_selfclaim']['full_prediction_jackknife_claim'])
    claim_terms = np.minimum(predjk, PRED_TRACE_CAP_K*refjk)
    claim = float(weights @ claim_terms)
    ratio = min(1.0, B/claim) if claim > 0 else 1.0
    official_n = mse_unbiased_capped(comparator='bulk_lognorm', pred_bulk=pb['bulk_lognorm'],
        real_bulk=rb['bulk_lognorm'], pred_moments=pm['bulk_lognorm'],
        real_moments=rm['bulk_lognorm'], genes=genes, pred_bulk_full=pb['bulk_lognorm'])
    official_d = distance_unbiased(comparator='bulk_lognorm', real_bulk=rb['bulk_lognorm'],
        real_moments=rm['bulk_lognorm'], genes=genes)
    saved = pd.read_csv(Path('artifacts/baseline_scores')/baseline.BASELINE_ID/task/'tidy.csv',
                        float_precision='round_trip')
    rows, nums, dens, plug = {}, [], [], []
    for i, target in enumerate(targets):
        keep = genes != target
        sq = float(np.dot(pred[i, keep]-refs[i, keep], pred[i, keep]-refs[i, keep]))
        ctrl_sq = float(np.dot(refs[i, keep]-refctrl[keep], refs[i, keep]-refctrl[keep]))
        numerator = (sq-ratio*claim_terms[i]-refjk[i])*weights[i]
        denominator = (ctrl_sq-refjk[i]-ctrl_jk)*weights[i]
        close(numerator, official_n[target]); close(denominator, official_d[target])
        for metric, expected in [('expr_mse_unbiased_capped', numerator),
                                 ('expr_distance_unbiased', denominator)]:
            saved_value = saved.loc[(saved.perturbation == target) & (saved.metric == metric), 'value']
            assert len(saved_value) == 1
            close(saved_value.iloc[0], expected)
        rows[target] = {'plugin_error': sq*weights[i], 'prediction_jk': predjk[i]*weights[i],
            'reference_jk': refjk[i]*weights[i], 'capped_prediction_refund': ratio*claim_terms[i]*weights[i],
            'numerator': numerator, 'denominator': denominator}
        nums.append(numerator); dens.append(denominator); plug.append(sq*weights[i])
    count = len(targets)
    numerator, denominator = float(np.mean(nums)), float(np.mean(dens))
    normalized = sum(nums)/sum(dens)
    for key, expected in [('expr_mse_unbiased_capped', numerator),
                          ('expr_distance_unbiased', denominator),
                          ('expr_mse_unbiased_capped_norm', normalized)]:
        close(record['raw'][key], expected)
    close(sum(plug), error_centre_ss+error_spread_ss)
    anchors = load(f'artifacts/scores/{task}/anchors.json')
    key = 'expr_mse_unbiased_capped_norm'
    base, rep = [float(anchors[x][key]) for x in ['baseline', 'replicate']]
    assert denominator > 0 and base > rep
    unbounded_scaled = (normalized-base)/(rep-base)
    close(min(1, max(0, unbounded_scaled)), record['scaled'][key])
    # Fixed observed bulk: ANY emitter's permitted refund <= min(B, k*sum(w*jk_real)).
    # This is a bound on emission-only improvement, never a candidate or new score.
    max_refund = min(B, float(weights @ (PRED_TRACE_CAP_K*refjk)))
    best_possible_num = float(np.mean(plug)) - float(weights @ refjk)/count - max_refund/count
    best_possible_norm = best_possible_num/denominator
    return {'reference_sha256': record['reference_sha256'], 'prediction_sha256': digest.hexdigest(),
        'actual_gene_count': len(genes), 'target_count': count, 'own_gene_unresolved':
            [t for t in targets if t not in gene_index], 'plugin_error_mean': float(np.mean(plug)),
        'reference_jk_subtraction_mean': float(weights @ refjk)/count,
        'prediction_selfclaim_total': float(weights @ predjk), 'prediction_B_total': B,
        'reference_capped_claim_total': claim, 'retained_fraction': ratio,
        'actual_prediction_refund_mean': ratio*claim/count,
        'numerator_mean': numerator, 'denominator_mean': denominator, 'normalized_MSE': normalized,
        'anchor_baseline_normalized_MSE': base, 'anchor_replicate_normalized_MSE': rep,
        'scaled_before_clamp': unbounded_scaled, 'scaled_after_clamp': record['scaled'][key],
        'plugin_error_centre_mean': error_centre_ss/count,
        'plugin_error_centred_residual_mean': error_spread_ss/count,
        'centre_error_fraction': error_centre_ss/sum(plug),
        'reference_B_total': ref_B, 'reference_centre_ss_total': ref_centre_ss,
        'reference_max_claim_total': float(weights @ (PRED_TRACE_CAP_K*refjk)),
        'reference_B_over_max_claim': ref_B/float(weights @ (PRED_TRACE_CAP_K*refjk)),
        'fixed_bulk_emission_only_bound': {'max_permitted_refund_mean': max_refund/count,
            'remaining_refund_headroom_mean': (max_refund-ratio*claim)/count,
            'optimistic_normalized_MSE_lower_bound': best_possible_norm,
            'optimistic_scaled_MSE_upper_bound': min(1, max(0, (best_possible_norm-base)/(rep-base))),
            'scope': 'Fixed observed pseudobulk; algebraic bound across possible within-target scatter, not an emitted candidate'},
        'rows': rows}


def main():
    assert not OUT.exists() and not RECEIPTS.exists(), 'Preserve immutable round receipts'
    plan = load(PLAN)
    assert plan['auditor_sha256'] == rt.sha(__file__)
    spec, result = check_receipts(plan)
    stat_cache.install()
    tasks = {}
    for task in baseline.TASKS:
        tasks[task] = diagnose(task, spec, result['tasks'][task])
        print('PASS MSE decomposition', task, rt.json.dumps({k: v for k, v in tasks[task].items()
                                                           if k != 'rows'}), flush=True)
        gc.collect()
    baseline.read_spec()
    for path, digest in plan['preserve_sha256'].items():
        assert rt.sha(path) == digest, path
    assert not Path('artifacts/holdout_events.jsonl').exists()
    rt.save(OUT, {'status': 'PASS same-count-stream official MSE and independent arithmetic',
        'utc': datetime.now(timezone.utc).isoformat(), 'auditor_sha256': rt.sha(__file__),
        'plan_sha256': rt.sha(PLAN), 'completed_score_audit_sha256': rt.sha(RECEIPTS),
        'tasks': tasks, 'phase2_allowed': False, 'holdout_scored': False,
        'limitations': ['Observed SEARCH attribution only; no causal control-centering ablation.',
            'Reference-derived decomposition is diagnostic and never a predictor input.',
            'No amplitude, emission, panel, anchor, role, historical gate or model selection change.',
            'No external official pairing or forecast; benchmark still fails historical rank gate.']})


if __name__ == '__main__':
    main()

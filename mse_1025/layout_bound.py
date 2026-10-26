"""Analytic limit on layout-only improvements under the official two correction caps.

This is a bound derived from saved official moments, not a replacement scorer or
a candidate prediction. Fixed group gene totals fix both plug-in error and the
whole-panel spread budget, regardless of the arrangement into cells.
"""
import ast
import hashlib
import json
import math
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def derive_layout_bound(decomposition):
    totals = decomposition['totals']
    denominator = totals['denominator_sum']
    assert denominator > 0, 'A positive denominator is required for this bound'
    source = ROOT / 'inputs/validated_delta.py.txt'
    expected = json.loads((ROOT/'inputs/scoring_provenance.json').read_text())
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    assert source_hash == expected['package_files']['metrics/delta.py']
    tree = ast.parse(source.read_text())
    constants = {
        target.id: ast.literal_eval(node.value)
        for node in tree.body if isinstance(node, ast.Assign)
        for target in node.targets
        if isinstance(target, ast.Name) and target.id == 'PRED_TRACE_CAP_K'
    }
    cap_k = constants['PRED_TRACE_CAP_K']
    assert cap_k == 1.0
    real_correction = totals['numerator_real_correction_sum']
    budget = decomposition['validation']['across_target_correction_budget']
    upper_correction = min(cap_k * real_correction, budget)
    lower_numerator = totals['numerator_plugin_sum'] - real_correction - upper_correction
    floor = lower_numerator / denominator
    control = decomposition['control_reference']['raw']['expr_mse_unbiased_capped_norm']
    actual = totals['numerator_effective_prediction_correction_sum']
    assert actual <= upper_correction + 1e-12
    assert floor <= totals['ratio'] + 1e-12
    assert math.isclose(
        totals['ratio'] - floor, (upper_correction - actual) / denominator,
        abs_tol=1e-12,
    )
    return dict(
        status='ANALYTIC LOWER BOUND; not an achieved score',
        scope='Same complete CD4 SEARCH target/gene panel and exact incumbent per-target group gene totals',
        official_prediction_trace_cap_k=cap_k,
        official_delta_sha256=source_hash,
        source='inputs/validated_delta.py.txt; byte-identical to the validated installed scorer',
        plugin_numerator_sum=totals['numerator_plugin_sum'],
        real_correction_sum=real_correction,
        correction_budget=budget,
        correction_upper_bound_sum=upper_correction,
        actual_correction_sum=actual,
        numerator_lower_bound_sum=lower_numerator,
        denominator_sum=denominator,
        raw_mse_lower_bound=floor,
        control_reference_raw_mse=control,
        lower_bound_excess_over_control=floor-control,
        best_possible_skill_against_control=1-floor/control,
        maximum_possible_raw_mse_improvement=totals['ratio']-floor,
        layout_only_can_beat_control=floor < control,
        proof='C_effective <= sum_p w_p min(C_pred_p, k*C_real_p) '
              '<= k*sum_p w_p C_real_p, and C_effective <= B_pred. '
              'Thus N >= plugin_error - C_real - min(k*C_real, B_pred). '
              'Fixed group gene totals fix plugin_error and B_pred; D is reference-only.',
        qualifications=[
            'This deliberately optimistic bound need not be attainable by any cell layout.',
            'Changing group gene totals changes the prediction and requires recomputing this bound.',
            'It bounds raw expression MSE only; no conclusion about the other five metrics or the hidden validation panel.',
        ],
    )


if __name__ == '__main__':
    assert os.getpriority(os.PRIO_PROCESS, 0) >= 10
    assert shutil.disk_usage(ROOT).free >= 33*1024**3, 'Proactive disk hold'
    path = ROOT / 'results/decomposition.json'
    decomposition = json.loads(path.read_text())
    bound = derive_layout_bound(decomposition)
    decomposition['layout']['fixed_totals_lower_bound'] = bound
    path.write_text(json.dumps(decomposition, indent=2, allow_nan=False)+'\n')
    receipt = dict(utc=datetime.now(timezone.utc).isoformat(), **bound)
    (ROOT/'results/layout_bound.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))

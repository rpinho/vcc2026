"""Independent saved-scale, identity and census audit; no model or response search."""
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import math

import runtime as rt
import numpy as np
import pandas as pd


def main():
    task = 'CD4_depth120_search'
    folder = Path('artifacts/scores')/task
    inputs = [folder/name for name in ['invalid_scale.json', 'anchor_splits.csv',
        'anchor_summary.csv', 'baseline_mean_response.json', 'baseline_mean_response_protocol.json']]
    invalid, baseline, protocol = [rt.json.loads((folder/name).read_text()) for name in [
        'invalid_scale.json', 'baseline_mean_response.json', 'baseline_mean_response_protocol.json']]
    setup_audit = rt.json.loads(Path('artifacts/CD4_depth120_activation_audit.json').read_text())
    assert invalid['reference_sha256'] == protocol['reference_sha256'] == setup_audit['reference_sha256']
    assert setup_audit['status'] == 'PASS' and setup_audit['anchor_splits'] == 5
    assert protocol['builder_sha256'] == rt.sha('baseline_protocol.py')
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    assert rt.sha(index['archive']) == index['archive_sha256']
    receipt = rt.json.loads(Path('artifacts/tasks', task+'.json').read_text())
    assert receipt['sha256'] == rt.sha(Path('artifacts/tasks', task+'.h5ad'))
    splits = pd.read_csv(folder/'anchor_splits.csv', float_precision='round_trip')
    summary = pd.read_csv(folder/'anchor_summary.csv', float_precision='round_trip').set_index('metric')
    seeds = np.random.SeedSequence(0).generate_state(5).astype(np.int64)
    for metric, rows in splits.groupby('metric'):
        rows = rows.sort_values('split_index')
        assert rows.split_index.tolist() == list(range(5))
        assert np.array_equal(rows.seed, seeds)
        values = rows.value.to_numpy()
        for key, value in [('replicate', values.mean()), ('replicate_sd', values.std(ddof=0)),
            ('replicate_min', values.min()), ('replicate_max', values.max())]:
            assert math.isclose(float(summary.loc[metric, key]), float(value), rel_tol=1e-13, abs_tol=1e-14), (metric, key, summary.loc[metric, key], value)
    from official_scoring import SIX
    headrooms = {}
    for metric in SIX:
        b = float(baseline[metric])
        r = float(summary.loc[metric, 'replicate'])
        delta = Fraction(r)-Fraction(b)
        positive = metric not in ['expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae']
        headrooms[metric] = dict(baseline=b, replicate=r,
            exact_binary64_headroom=str(delta), usable_direction=bool(delta > 0 if positive else delta < 0),
            split_values=splits.loc[splits.metric == metric, 'value'].tolist())
    failed = [key for key, entry in headrooms.items() if not entry['usable_direction']]
    assert failed == ['de_wilcoxon_direction_fidelity_yield_raw']
    assert failed[0] in invalid['error'] and 'no usable replicate scale' in invalid['error']
    assert not (folder/'anchors.json').exists() and not list(folder.glob('*_s42.json'))
    # Independently recalculate count-only availability, including a depth-400
    # alternative's feasibility. This does not select or register a new task.
    census = pd.read_parquet('artifacts/CD4_cell_census.parquet', columns=['target_gene', 'batch'])
    table = census.groupby(['target_gene', 'batch'], observed=True).size().unstack(fill_value=0)
    controls = table.loc['non-targeting']
    table = table.drop(index='non-targeting')
    contract = rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())[task]
    fixed = table.loc[contract['targets']]
    assert len(fixed) == 300 and int(fixed.min(axis=0).sum()) == 120
    upper = int(np.sort(table.to_numpy(), axis=0)[-300].sum())
    assert upper == 391
    gwps = pd.read_parquet('artifacts/GWPS_raw_obs.parquet', columns=['gene'])
    gwps_n = gwps.groupby('gene', observed=True).size()
    available = table.index.intersection(gwps_n[gwps_n >= 200].index)
    # Balanced 400-cell quota is 16 per lane plus one extra in 16 lanes.
    # For this descriptive feasibility bound, require 17 in every lane (408).
    capacity = []
    for q in [5, 8, 10, 17]:
        eligible = (table.loc[available] >= q).all(axis=1)
        capacity.append(dict(per_lane_quota=q, real_cells=24*q,
            matched_GWPS200_targets=int(eligible.sum())))
    rt.save('artifacts/CD4_depth120_scale_audit.json', dict(
        status='PASS independent invalid-scale and count-capacity audit',
        utc=datetime.now(timezone.utc).isoformat(), verifier_sha256=rt.sha(__file__),
        task=task, reference_sha256=invalid['reference_sha256'],
        metadata_members_verified=len(index['entries']), headrooms=headrooms,
        input_hashes={str(path): rt.sha(path) for path in inputs},
        official_five_seeds=seeds.tolist(),
        all_cd4_300_target_common_quota_upper_bound=upper,
        fixed_300_target_maximum_common_depth=120,
        count_only_alternative_capacity=capacity,
        controls_per_lane_min=int(controls.min()),
        new_task_registered=False, holdout_scored=False,
        conclusion='Frozen 300-target depth repair exhausted. Exclude the whole task; gate CLOSED. Further work requires a separately frozen response-independent design, not different splits or thresholds.',
        limitation='Two of five fidelity splits exceed baseline; selecting those would be biased. More depth may not fix the scale. Descriptive alternative capacity is not a task selection or evidence of predictive validity.'))
    print('PASS CD4 depth120 invalid scale; fixed depth exhausted; census alternatives', capacity, flush=True)


if __name__ == '__main__':
    main()

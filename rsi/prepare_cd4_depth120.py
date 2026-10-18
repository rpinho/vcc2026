"""Audit the failed CD4 scale and freeze a count-only, incremental depth repair.

No response matrix is read except for hashing the already materialized base task.
No model score, HOLDOUT outcome, DE gene list or response magnitude selects rows.
"""
import runtime as rt
from pathlib import Path
from fractions import Fraction
from itertools import combinations
from datetime import datetime, timezone
import copy
import numpy as np
import pandas as pd

NT = 'non-targeting'
BASE = 'CD4_matched_search'
NAME = 'CD4_depth120_search'
FOLDER = Path('artifacts/tasks')


def upper_bound(counts, n):
    """Any n-target common quota is bounded by each lane's nth largest count."""
    return int(np.sort(counts, axis=0)[-n].sum())


def main():
    metadata = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, entry in metadata['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    # Exhaustive independent fixture: the integer bound never understates capacity.
    example = np.array([[1, 7], [3, 2], [8, 0], [2, 4]], dtype=np.int64)
    for n in range(1, 5):
        exact = max(sum(min(example[i, b] for i in rows) for b in range(2))
                    for rows in combinations(range(4), n))
        assert upper_bound(example, n) >= exact

    scores = Path('artifacts/scores') / BASE
    invalid = rt.json.loads((scores/'invalid_scale.json').read_text())
    preflight = rt.json.loads(Path('artifacts/calibration_restore_preflight.json').read_text())
    assert invalid['reference_sha256'] == preflight['tasks'][BASE]['reference_sha256']
    baseline = rt.json.loads((scores/'baseline_mean_response.json').read_text())
    splits = pd.read_csv(scores/'anchor_splits.csv')
    summary = pd.read_csv(scores/'anchor_summary.csv').set_index('metric')
    seeds = np.random.SeedSequence(0).generate_state(3).tolist()
    headrooms = {}
    from verify_results import SIX
    for metric in SIX:
        part = splits[splits.metric == metric].sort_values('split_index')
        assert part.split_index.tolist() == [0, 1, 2]
        assert part.seed.tolist() == seeds
        values = part.value.to_numpy(float)
        anchor = float(summary.loc[metric, 'replicate'])
        assert abs(float(np.mean(values)) - anchor) < 1e-14
        assert abs(float(np.std(values)) - float(summary.loc[metric, 'replicate_sd'])) < 1e-14
        delta = Fraction(anchor) - Fraction(float(baseline[metric]))
        direction = -1 if metric in ['expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae'] else 1
        headrooms[metric] = dict(baseline=baseline[metric], replicate=anchor,
                                exact_binary64_headroom=str(delta),
                                usable_direction=bool(direction * delta > 0),
                                split_values=values.tolist())
    fidelity = 'de_wilcoxon_direction_fidelity_yield_raw'
    assert not headrooms[fidelity]['usable_direction']
    assert all(v < baseline[fidelity] for v in headrooms[fidelity]['split_values'])
    assert not (scores/'anchors.json').exists(), 'Invalid task must not have usable anchors'

    census = pd.read_parquet('artifacts/CD4_cell_census.parquet')
    table = census.groupby(['target_gene', 'batch'], observed=True).size().unstack(fill_value=0)
    table.columns = table.columns.astype(str)
    manifest = rt.json.loads((FOLDER/'manifest.json').read_text())
    base = manifest[BASE]
    targets, batches = sorted(base['targets']), sorted(table.columns)
    capacity = table.loc[targets, batches].min(axis=0).to_numpy(np.int64)
    assert np.array_equal(capacity, np.full(24, 5)), capacity
    all_counts = table.drop(index=NT).to_numpy(np.int64)
    global_bound = upper_bound(all_counts, 300)
    assert global_bound < 400
    # Keep targets, controls and every acquired response row. Only add one cell/lane.
    old_path = FOLDER/f'{BASE}_selection.parquet'
    old = pd.read_parquet(old_path)
    assert rt.sha(old_path) == base['selection_sha256']
    inherited = set(old.source_row)
    pool = census[census.target_gene.isin(targets) & ~census.source_row.isin(inherited)]
    groups = pool.groupby(['target_gene', 'batch'], observed=True)
    rng = np.random.default_rng(20261017120)
    additions = []
    for target in targets:
        for batch in batches:
            group = groups.get_group((target, batch)).sort_values('source_row')
            row = group.iloc[[int(rng.integers(len(group)))]].copy()
            row['evaluation_role'] = 'perturbation'
            additions.append(row)
    extra = pd.concat(additions, ignore_index=True)[old.columns]
    selection = pd.concat([old, extra], ignore_index=True).sort_values('source_row').reset_index(drop=True)
    assert len(extra) == 7200 and selection.source_row.is_unique
    counts = selection[selection.evaluation_role == 'perturbation'].groupby(
        ['target_gene', 'batch'], observed=True).size()
    assert len(counts) == 7200 and (counts == 5).all()
    for role in ['input_control', 'reference_control']:
        pd.testing.assert_frame_equal(
            selection[selection.evaluation_role == role].sort_values('source_row').reset_index(drop=True),
            old[old.evaluation_role == role].sort_values('source_row').reset_index(drop=True))
    base_receipt = rt.json.loads((FOLDER/f'{BASE}.json').read_text())
    assert rt.sha(FOLDER/f'{BASE}.h5ad') == base_receipt['sha256']
    journal = Path('artifacts/CD4_destination_charges.jsonl')
    charges = [rt.json.loads(line) for line in journal.read_text().splitlines()]
    spent = sum(row['bytes'] for row in charges if row['event'] == 'charge')
    remaining = 2_600_000_000 - spent
    assert remaining > 251_000_000
    plan = copy.deepcopy(base)
    plan.update(n_real=120, source='CD4_raw_maximum_frozen_target_depth',
                selection_policy='Same 300 frozen targets; maximal common 24-lane quota from census = 5. Keep all old rows and controls; add one count-selected cell per target/lane.',
                parent_task=BASE, parent_task_sha256=base_receipt['sha256'])
    plan['audit'].update(target_batch_quotas={b: 5 for b in batches},
                         control_ratio=32/5, input_control_ratio=16,
                         reference_control_ratio=32/5)
    plan['control_revision'] = 'Exactly the original disjoint input/reference controls, unchanged; control mixtures still equal each target mixture.'
    plan.pop('budget_fallback', None)
    plan.pop('budget_fallback_used', None)
    output = FOLDER/f'{NAME}_selection.parquet'
    increment = FOLDER/f'{NAME}_increment.parquet'
    blueprint = FOLDER/f'{NAME}_blueprint.json'
    if blueprint.exists():
        previous = rt.json.loads(blueprint.read_text())
        assert rt.sha(output) == previous['selection_sha256']
        pd.testing.assert_frame_equal(pd.read_parquet(output), selection)
        pd.testing.assert_frame_equal(pd.read_parquet(increment), extra)
    else:
        selection.to_parquet(output, index=False)
        extra.to_parquet(increment, index=False)
        plan.update(selection_sha256=rt.sha(output), increment_sha256=rt.sha(increment),
                    original_charge_journal_sha256=rt.sha(journal),
                    original_charged_bytes=spent, increment_ceiling_bytes=remaining,
                    increment_seed=20261017120,
                    qualification='Depth repair frozen before new anchors/models. Valid scale and calibration ranks remain UNVERIFIED. No fallback quota or target reselection.')
        rt.save(blueprint, plan)
    result = dict(status='PASS independent invalid-scale and response-blind depth-design audit',
                  utc=datetime.now(timezone.utc).isoformat(), round=13,
                  metadata_members_verified=len(metadata['entries']),
                  reference_sha256=invalid['reference_sha256'], headrooms=headrooms,
                  input_hashes={str(p): rt.sha(p) for p in [scores/'invalid_scale.json', scores/'anchor_splits.csv', scores/'anchor_summary.csv', scores/'baseline_mean_response.json', Path('artifacts/CD4_cell_census.parquet')]},
                  all_cd4_300_target_common_quota_upper_bound=global_bound,
                  upper_bound_fixture='All subsets of a 4-target, 2-lane synthetic census verified independently',
                  fixed_300_target_maximum_common_depth=int(capacity.sum()),
                  next_task=NAME, additional_cells=7200, original_charged_bytes=spent,
                  increment_ceiling_bytes=remaining, blueprint_sha256=rt.sha(blueprint),
                  limitation='Fidelity failure cause is not established causally. More cells may still leave an invalid scale. Integer bound assumes one common lane mixture and no repeated cells.',
                  holdout_scored=False)
    rt.save('artifacts/CD4_scale_design_audit.json', result)
    print(result['status'], 'all-CD4 300-target bound', global_bound,
          'fixed-target maximum', capacity.sum(), 'increment ceiling', remaining, flush=True)


if __name__ == '__main__':
    main()

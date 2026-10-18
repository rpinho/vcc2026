"""Freeze a count-only, deeper CD4 panel without reading expression or scores.

The target count matches H1 SEARCH. Maximize a uniform quota over all 24 lanes
subject to at least 79 targets with >=200 cells in both prepared GWPS banks.
The resulting small-panel PDS remains a stated limitation, never a 300-panel proxy.
"""
from datetime import datetime, timezone
import hashlib
from pathlib import Path

import runtime as rt
import numpy as np
import pandas as pd

NAME = 'CD4_depth192_search'
BASE = 'CD4_depth120_search'
FOLDER = Path('artifacts/tasks')
SALT = 'VCC-RSI-CD4-matched-depth-20261017-v1:'
SEED = 20261017


def main():
    revision_path = Path('artifacts/benchmark_depth192_revision.json')
    assert not revision_path.exists(), 'Design already frozen; do not reselect'
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    census_path = Path('artifacts/CD4_cell_census.parquet')
    census = pd.read_parquet(census_path)
    assert census.source_row.is_unique and not census.low_quality.any()
    table = pd.crosstab(census.target_gene, census.batch).sort_index(axis=1)
    batches = table.columns.astype(str).tolist()
    assert len(batches) == 24
    table = table.drop('non-targeting')
    bank_paths = ['artifacts/cd4_banks/K562_GWPS_raw.npz',
                  'artifacts/cd4_sources/K562_GWPS_bulk.npz']
    available = set(table.index)
    for path in bank_paths:
        with np.load(path) as bank:
            available &= set(bank['targets'][bank['n_cells'] >= 200].astype(str))
    counts = table.loc[sorted(available)]
    capacity = {q: int((counts >= q).all(axis=1).sum())
                for q in range(1, int(counts.min(axis=1).max())+2)}
    quota = max(q for q, n in capacity.items() if n >= 79)
    assert quota == 8, 'Unexpected census; investigate without looking at expression'
    eligible = counts.index[(counts >= quota).all(axis=1)].astype(str).tolist()
    targets = sorted(eligible, key=lambda t: hashlib.sha256((SALT+t).encode()).hexdigest())[:79]
    assert len(eligible) == 84 and capacity[9] == 52
    base_selection = pd.read_parquet(FOLDER/f'{BASE}_selection.parquet')
    inherited = base_selection.loc[(base_selection.evaluation_role != 'perturbation') |
                                  base_selection.target_gene.isin(targets)].copy()
    rng = np.random.default_rng(SEED)
    pieces = [inherited]
    for target in targets:
        for batch in batches:
            old = inherited.loc[(inherited.target_gene == target) & (inherited.batch == batch)]
            needed = quota-len(old)
            assert needed >= 0
            pool = census.loc[(census.target_gene == target) & (census.batch == batch) &
                              ~census.source_row.isin(old.source_row)]
            chosen = pool.iloc[rng.choice(len(pool), needed, replace=False)].copy()
            chosen['evaluation_role'] = 'perturbation'
            pieces.append(chosen)
    selection = pd.concat(pieces, ignore_index=True).sort_values('source_row').reset_index(drop=True)
    increment = selection.loc[~selection.source_row.isin(inherited.source_row)].copy()
    assert selection.source_row.is_unique and len(selection) == 79*192+1920+768
    selection_path = FOLDER/f'{NAME}_selection.parquet'
    increment_path = FOLDER/f'{NAME}_increment.parquet'
    selection.to_parquet(selection_path, index=False)
    increment.to_parquet(increment_path, index=False)
    manifest = rt.json.loads((FOLDER/'manifest.json').read_text())
    base_plan = manifest[BASE]
    plan = dict(base_plan)
    plan.update(context='CD4', split='search', targets=targets, cells_per_target=192, n_real=192,
        source='CD4_Stim8hr_matched_depth192', replicate_splits=5,
        selection_sha256=rt.sha(selection_path), increment_sha256=rt.sha(increment_path),
        parent_task=BASE, parent_task_sha256=rt.sha(FOLDER/f'{BASE}.h5ad'),
        parent_selection_sha256=rt.sha(FOLDER/f'{BASE}_selection.parquet'),
        census_sha256=rt.sha(census_path), eligible_for_primary_calibration=True,
        qualification='UNVERIFIED count-only design; activation and six-scale audit required',
        panel_size=79, per_lane_quota=8, target_order_salt=SALT, selection_seed=SEED,
        additional_cells=len(increment), inherited_cells=len(inherited),
        increment_ceiling_bytes=1_200_000_000,
        previous_charge_journals={p: rt.sha(p) for p in [
            'artifacts/CD4_destination_charges.jsonl', 'artifacts/CD4_depth120_charges.jsonl']})
    # Do not inherit supersession/failed-scale declarations from the invalid parent.
    for key in ['superseded_by', 'scale_status', 'invalid_scale', 'anchor_policy_revision',
                'original_charge_journal_sha256', 'original_charged_bytes', 'increment_seed']:
        plan.pop(key, None)
    plan['matching'] = dict(exact_target_batch_matching=True,
        target_batch_quotas={b: 8 for b in batches}, input_controls_per_batch=80,
        reference_controls_per_batch=32, target_to_input_ratio='1/10', target_to_reference_ratio='1/4')
    plan['audit'] = dict(target_batch_quotas={b: 8 for b in batches},
        control_ratio=4, control_batch_quotas={b: 32 for b in batches},
        per_target_control_mixture_L1={t: 0.0 for t in targets},
        input_reference_disjoint=True, input_control_ratio=10, reference_control_ratio=4,
        input_control_batch_quotas={b: 80 for b in batches})
    plan['selection_policy'] = 'Count-only maximum uniform quota supporting79 covered targets; SHA256 first79; retain available parent counts and sample missing cells once.'
    blueprint_path = FOLDER/f'{NAME}_blueprint.json'
    rt.save(blueprint_path, plan)
    revision = dict(frozen_utc=datetime.now(timezone.utc).isoformat(), round=17,
        status='FROZEN count-only design; expressions/scales/ranks UNVERIFIED',
        reason='300-target CD4 maximum common depth120 failed its frozen scale. '
               'Match H1 SEARCH panel size79 and maximize uniform actual depth over all24 lanes '
               'with >=200-cell evidence in both existing GWPS banks, without new source preparation.',
        primary_search=['H1_rich_search', NAME], anchors={'H1_rich_search': 5, NAME: 5},
        task=NAME, parent=BASE, blueprint_sha256=rt.sha(blueprint_path),
        census_sha256=rt.sha(census_path), bank_sha256={p: rt.sha(p) for p in bank_paths},
        capacity_by_quota=capacity, eligible_targets=sorted(eligible),
        targets=targets, target_order_salt=SALT, selection_seed=SEED,
        count_objective='Maximum uniform all-lane quota with >=79 source-covered targets; SHA256 first79',
        controls='Exactly the existing1920 input and768 reference controls, no HOLDOUT controls',
        policy='If any of six scales fails, exclude this whole task. No favourable splits, '
               'alternative targets, dropped lanes, threshold changes, scale clipping or smaller-panel retries. '
               'No recipe score until all six five-split scales pass. Gate stays CLOSED until two-context '
               'all-pair ranks and external baseline requirements pass.',
        activation_requirement='Independent materialization audit PASS; archive manifest, roles and '
               'anchor policy before activation. Historical frozen metadata/archive remain unchanged.',
        limitations='79 targets per context gives small-panel PDS and a count-enriched gene panel; '
                     '192 real versus400 emitted cells remains unequal-depth. No guarantee of valid scales '
                     'or final300-target unseen-cell-line generalization. Previously observed aggregate '
                     'CD4 failures motivated depth change; no target-level responses/DE/model scores used.',
        holdout_scored=False, expression_payloads_read=False, task_registered=False)
    rt.save(revision_path, revision)
    print('FROZEN', NAME, 'eligible', len(eligible), 'targets', len(targets),
          'quota', quota, 'inherited', len(inherited), 'new', len(increment), flush=True)


if __name__ == '__main__':
    main()

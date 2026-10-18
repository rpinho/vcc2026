"""Independent metadata design verification; optional bounded live count audit."""
from collections import Counter
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import sys

import runtime as rt
import numpy as np
import pandas as pd

NAME = 'CD4_depth192_search'
BASE = 'CD4_depth120_search'
FOLDER = Path('artifacts/tasks')


def main(materialized=False):
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    assert rt.sha(index['archive']) == index['archive_sha256']
    revision_path = Path('artifacts/benchmark_depth192_revision.json')
    revision = rt.json.loads(revision_path.read_text())
    blueprint_path = FOLDER/f'{NAME}_blueprint.json'
    plan = rt.json.loads(blueprint_path.read_text())
    assert revision['blueprint_sha256'] == rt.sha(blueprint_path)
    assert plan['n_real'] == plan['cells_per_target'] == 192
    assert plan['n_pred'] == 400 and plan['replicate_splits'] == 5
    census_path = Path('artifacts/CD4_cell_census.parquet')
    assert rt.sha(census_path) == revision['census_sha256'] == plan['census_sha256']
    census = pd.read_parquet(census_path)
    assert census.source_row.is_unique and not census.low_quality.any()
    row_groups = {key: value for key, value in census.groupby(['target_gene', 'batch'], observed=True, sort=False)}
    available = None
    for path, identity in revision['bank_sha256'].items():
        assert rt.sha(path) == identity
        with np.load(path) as bank:
            targets = {str(t) for t, n in zip(bank['targets'], bank['n_cells']) if n >= 200}
        available = targets if available is None else available & targets
    # Counter and scalar minima independently check pandas' design calculation.
    counter = Counter(zip(census.target_gene.astype(str), census.batch.astype(str)))
    batches = sorted(census.batch.astype(str).unique())
    minima = {t: min(counter[t, b] for b in batches) for t in available}
    eligible = sorted(t for t, n in minima.items() if n >= 8)
    assert len(batches) == 24 and len(eligible) == 84
    assert sum(n >= 9 for n in minima.values()) == 52
    for q, n in revision['capacity_by_quota'].items():
        assert sum(v >= int(q) for v in minima.values()) == n
    assert max(q for q in range(1, max(minima.values())+1)
               if sum(n >= q for n in minima.values()) >= 79) == 8
    targets = sorted(eligible, key=lambda t: hashlib.sha256(
        (revision['target_order_salt']+t).encode()).hexdigest())[:79]
    assert targets == revision['targets'] == plan['targets'] and len(targets) == 79
    assert sorted(eligible) == revision['eligible_targets']
    selection_path = FOLDER/f'{NAME}_selection.parquet'
    increment_path = FOLDER/f'{NAME}_increment.parquet'
    assert rt.sha(selection_path) == plan['selection_sha256']
    assert rt.sha(increment_path) == plan['increment_sha256']
    selection = pd.read_parquet(selection_path).sort_values('source_row').reset_index(drop=True)
    increment = pd.read_parquet(increment_path)
    base_path = FOLDER/f'{BASE}_selection.parquet'
    assert rt.sha(base_path) == plan['parent_selection_sha256']
    base = pd.read_parquet(base_path)
    inherited = base[(base.evaluation_role != 'perturbation') | base.target_gene.isin(targets)]
    assert selection.source_row.is_unique and len(selection) == 17856
    assert len(inherited) == plan['inherited_cells']
    assert len(increment) == plan['additional_cells']
    assert set(selection.source_row) == set(inherited.source_row) | set(increment.source_row)
    assert set(inherited.source_row).isdisjoint(increment.source_row)
    selected = selection.set_index('source_row')
    full_census = census.set_index('source_row')
    for col in ['target_gene', 'batch', 'guide_id']:
        assert selected[col].astype(str).equals(full_census.loc[selected.index, col].astype(str))
    current_counts = Counter(zip(selection.target_gene, selection.batch, selection.evaluation_role))
    for t in targets:
        assert [current_counts[t, b, 'perturbation'] for b in batches] == [8]*24
    assert sorted(selection.loc[selection.evaluation_role == 'perturbation', 'target_gene'].unique()) == sorted(targets)
    for role, n in [('input_control', 80), ('reference_control', 32)]:
        assert set(selection.loc[selection.evaluation_role == role, 'source_row']) == set(base.loc[base.evaluation_role == role, 'source_row'])
        assert [current_counts['non-targeting', b, role] for b in batches] == [n]*24
    # A separately replayed RNG checks all newly selected row identities/order.
    rng = np.random.default_rng(revision['selection_seed'])
    replay = []
    for t in targets:
        for b in batches:
            old_rows = set(inherited.loc[(inherited.target_gene == t) & (inherited.batch == b), 'source_row'])
            group = row_groups[t, b]
            candidates = group.loc[~group.source_row.isin(old_rows)]
            replay.extend(candidates.iloc[rng.choice(len(candidates), 8-len(old_rows), replace=False)].source_row)
    assert sorted(replay) == sorted(increment.source_row)
    for path, expected in plan['previous_charge_journals'].items():
        assert rt.sha(path) == expected
    assert rt.sha(FOLDER/f'{BASE}.h5ad') == plan['parent_task_sha256']
    audit = dict(status='PASS independent count-only depth192 design audit',
        utc=datetime.now(timezone.utc).isoformat(), verifier_sha256=rt.sha(__file__),
        blueprint_sha256=rt.sha(blueprint_path), revision_sha256=rt.sha(revision_path),
        frozen_metadata_members_verified=len(index['entries']), task=NAME,
        targets=79, real_cells_per_target=192, inherited_cells=len(inherited), additional_cells=len(increment),
        input_controls=1920, reference_controls=768, eligible_targets=84, next_quota_capacity=52,
        holdout_scored=False, expression_payloads_read=False,
        limitation='Count selection only; new scales/ranks and external baseline reproduction unresolved.')
    if not materialized:
        rt.save('artifacts/CD4_depth192_design_audit.json', audit)
        print(audit['status'], 'new cells', len(increment), flush=True)
        return
    import anndata as ad
    path = FOLDER/f'{NAME}.h5ad'
    receipt = rt.json.loads(path.with_suffix('.json').read_text())
    assert rt.sha(path) == receipt['sha256']
    assert receipt['blueprint_sha256'] == rt.sha(blueprint_path)
    assert receipt['selection_sha256'] == rt.sha(selection_path)
    parent = ad.read_h5ad(FOLDER/f'{BASE}.h5ad', backed='r')
    child = ad.read_h5ad(path, backed='r')
    try:
        assert child.shape == (17856, 17772)
        assert child.var_names.tolist() == parent.var_names.tolist()
        official = pd.read_csv(rt.ROOT.parent/'data/gene_names.csv').iloc[:, 0].astype(str)
        assert [g for g in official if g in set(child.var_names)] == child.var_names.tolist()
        for col in ['source_row', 'target_gene', 'batch', 'guide_id', 'evaluation_role']:
            assert child.obs[col].astype(str).tolist() == selection[col].astype(str).tolist()
        parent_indices = pd.Index(parent.obs.source_row).get_indexer(inherited.source_row)
        child_indices = pd.Index(child.obs.source_row).get_indexer(inherited.source_row)
        assert np.all(parent_indices >= 0) and np.all(child_indices >= 0)
        for lo in range(0, len(inherited), 256):
            difference = child.X[child_indices[lo:lo+256]]-parent.X[parent_indices[lo:lo+256]]
            difference.eliminate_zeros()
            assert difference.nnz == 0
        for lo in range(0, child.n_obs, 256):
            block = child.X[lo:lo+256]
            assert np.isfinite(block.data).all() and np.all(block.data >= 0)
            assert np.equal(block.data, np.floor(block.data)).all()
            totals = np.asarray(block.sum(1)).ravel()
            assert np.all(totals > 0)
            assert np.all(totals <= full_census.loc[child.obs.source_row.iloc[lo:lo+256], 'total_counts'].to_numpy())
    finally:
        child.file.close()
        parent.file.close()
    charged = 0
    journal = Path('artifacts/CD4_depth192_charges.jsonl')
    for line in journal.read_text().splitlines():
        row = rt.json.loads(line)
        if row['event'] == 'charge':
            assert row['bytes'] > 0
            charged += row['bytes']
            assert row['total'] == charged
    assert charged == receipt['bytes'] <= plan['increment_ceiling_bytes']
    extension_path = FOLDER/f'{NAME}_manifest_extension.json'
    assert rt.json.loads(extension_path.read_text()) == {NAME: plan}
    fixture = rt.json.loads(Path('artifacts/cd4_depth192_fixture.json').read_text())
    assert fixture['status'] == 'PASS' and fixture['reader_sha256'] == rt.sha('acquire_cd4_depth192.py')
    audit.update(status='PASS independent depth192 materialization audit',
        expression_payloads_read=True, sha256=receipt['sha256'], charged_bytes=charged,
        manifest_extension_sha256=rt.sha(extension_path), readouts=17772)
    rt.save('artifacts/CD4_depth192_materialization_audit.json', audit)
    print(audit['status'], flush=True)


if __name__ == '__main__':
    main('--materialized' in sys.argv)

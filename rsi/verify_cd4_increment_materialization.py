"""Independently audit the live depth repair, without computing any response score."""
import runtime as rt
from pathlib import Path
from datetime import datetime, timezone
import anndata as ad
import numpy as np
import pandas as pd

NAME = 'CD4_depth120_search'
BASE = 'CD4_matched_search'
FOLDER = Path('artifacts/tasks')


def main():
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    blueprint_path = FOLDER/f'{NAME}_blueprint.json'
    blueprint = rt.json.loads(blueprint_path.read_text())
    design = rt.json.loads(Path('artifacts/CD4_scale_design_audit.json').read_text())
    assert design['blueprint_sha256'] == rt.sha(blueprint_path)
    receipt = rt.json.loads((FOLDER/f'{NAME}.json').read_text())
    assert receipt['sha256'] == rt.sha(FOLDER/f'{NAME}.h5ad')
    assert receipt['blueprint_sha256'] == rt.sha(blueprint_path)
    assert receipt['selection_sha256'] == blueprint['selection_sha256']
    assert rt.sha(FOLDER/f'{BASE}.h5ad') == blueprint['parent_task_sha256']
    base = ad.read_h5ad(FOLDER/f'{BASE}.h5ad', backed='r')
    result = ad.read_h5ad(FOLDER/f'{NAME}.h5ad', backed='r')
    try:
        selection_path = FOLDER/f'{NAME}_selection.parquet'
        assert rt.sha(selection_path) == blueprint['selection_sha256']
        selection = pd.read_parquet(selection_path).sort_values('source_row').reset_index(drop=True)
        assert selection.source_row.is_unique and len(selection) == result.n_obs == 38688
        columns = ['source_row', 'target_gene', 'batch', 'guide_id', 'evaluation_role']
        for column in columns:
            assert np.array_equal(result.obs[column].astype(str).to_numpy(), selection[column].astype(str).to_numpy()), column
        assert np.array_equal(base.var_names.to_numpy(str), result.var_names.to_numpy(str))
        official = pd.read_csv(rt.ROOT.parent/'data/gene_names.csv').iloc[:, 0].astype(str).tolist()
        measured = set(result.var_names)
        assert [gene for gene in official if gene in measured] == result.var_names.tolist()
        assert result.n_vars == 17772
        observations = result.obs.reset_index(drop=True)
        old_observations = base.obs.reset_index(drop=True)
        inherited_order = pd.Index(observations.source_row).get_indexer(old_observations.source_row)
        assert np.all(inherited_order >= 0) and len(np.unique(inherited_order)) == len(inherited_order)
        for column in columns:
            assert np.array_equal(observations.iloc[inherited_order][column].astype(str).to_numpy(), old_observations[column].astype(str).to_numpy())
        # Bounded independent exact comparison of every inherited response/control count.
        for lo in range(0, base.n_obs, 256):
            high = min(lo+256, base.n_obs)
            difference = result.X[inherited_order[lo:high]] - base.X[lo:high]
            difference.eliminate_zeros()
            assert difference.nnz == 0
        for lo in range(0, result.n_obs, 256):
            matrix = result.X[lo:lo+256]
            assert np.isfinite(matrix.data).all() and np.all(matrix.data >= 0)
            assert np.equal(matrix.data, np.floor(matrix.data)).all()
            totals = np.asarray(matrix.sum(axis=1)).ravel()
            assert np.all((totals > 0) & (totals <= 1_000_000))
        perturbations = observations[observations.evaluation_role.astype(str) == 'perturbation']
        targets = sorted(blueprint['targets'])
        quotas = perturbations.groupby(['target_gene', 'batch'], observed=True).size()
        assert len(targets) == 300 and len(quotas) == 7200 and (quotas == 5).all()
        assert sorted(perturbations.target_gene.astype(str).unique()) == targets
        assert (perturbations.groupby('target_gene', observed=True).size() == 120).all()
        batches = sorted(perturbations.batch.astype(str).unique())
        for role, n_per_lane in [('input_control', 80), ('reference_control', 32)]:
            current = observations[observations.evaluation_role.astype(str) == role]
            previous = old_observations[old_observations.evaluation_role.astype(str) == role]
            assert set(current.source_row) == set(previous.source_row)
            assert (current.groupby('batch', observed=True).size().reindex(batches, fill_value=0) == n_per_lane).all()
        added = observations[~observations.source_row.isin(old_observations.source_row)]
        assert len(added) == 7200 and (added.evaluation_role.astype(str) == 'perturbation').all()
        expected_extra = pd.read_parquet(FOLDER/f'{NAME}_increment.parquet')
        assert set(added.source_row) == set(expected_extra.source_row)
        census = pd.read_parquet('artifacts/CD4_cell_census.parquet').set_index('source_row')
        for column in ['target_gene', 'batch', 'guide_id']:
            assert np.array_equal(census.loc[added.source_row, column].astype(str).to_numpy(), added[column].astype(str).to_numpy())
    finally:
        base.file.close()
        result.file.close()
    journals = ['artifacts/CD4_destination_charges.jsonl', 'artifacts/CD4_depth120_charges.jsonl']
    totals = []
    for journal in journals:
        total = 0
        for line in Path(journal).read_text().splitlines():
            row = rt.json.loads(line)
            if row['event'] == 'charge':
                assert row['bytes'] > 0
                total += row['bytes']
                assert total == row['total']
        totals.append(total)
    assert totals[0] == blueprint['original_charged_bytes']
    assert totals[1] == receipt['bytes'] and sum(totals) == receipt['cumulative_cd4_bytes'] <= 2_600_000_000
    extension_path = FOLDER/f'{NAME}_manifest_extension.json'
    extension = rt.json.loads(extension_path.read_text())
    assert extension == {NAME: blueprint}
    fixture = rt.json.loads(Path('artifacts/cd4_increment_fixture.json').read_text())
    assert fixture['reader_sha256'] == rt.sha('acquire_cd4_increment.py') and fixture['status'] == 'PASS'
    output = dict(status='PASS independent count/selection/control/inheritance/budget audit',
                  utc=datetime.now(timezone.utc).isoformat(), task=NAME,
                  sha256=receipt['sha256'], blueprint_sha256=rt.sha(blueprint_path),
                  frozen_metadata_members_verified=len(index['entries']),
                  targets=300, cells_per_target=120, inherited_cells=31488, additional_cells=7200,
                  input_controls=1920, reference_controls=768, readouts=17772,
                  cumulative_cd4_bytes=sum(totals), ceiling_bytes=2_600_000_000,
                  manifest_extension_sha256=rt.sha(extension_path),
                  holdout_scored=False,
                  limitation='Counts and provenance only. No new baseline, replicate anchor, recipe score or scale validity was evaluated. Frozen global manifest/roles remain unchanged pending explicit design revision.')
    rt.save('artifacts/CD4_depth120_materialization_audit.json', output)
    print(output['status'], 'CD4 charged', sum(totals), flush=True)


if __name__ == '__main__':
    main()

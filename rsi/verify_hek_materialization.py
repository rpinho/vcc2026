"""Read-only frozen HEK task/count/transfer audit; never score HOLDOUT."""
import argparse
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

NAME = 'HEK293T_holdout'


def verify(projected=False):
    folder = Path('artifacts/tasks')
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for file, expected in index['entries'].items():
        assert Path(file).stat().st_size == expected['bytes'] and rt.sha(file) == expected['sha256'], file
    manifest = rt.json.loads((folder/'manifest.json').read_text())[NAME]
    assert manifest['split'] == 'holdout' and manifest['context'] == 'HEK293T'
    selection = folder/(NAME+'_selection.parquet')
    selected = pd.read_parquet(selection).sort_values('cell_integer_id').reset_index(drop=True)
    assert selected.cell_integer_id.is_unique and selected.source_row.is_unique
    assert selected.pass_guide_filter.all()
    roles = selected.evaluation_role.astype(str)
    assert roles.value_counts().to_dict() == dict(perturbation=19200, input_control=2000, reference_control=2000)
    response = selected.loc[roles == 'perturbation']
    assert len(manifest['targets']) == 300 and manifest['n_real'] == 64
    assert set(response.target_gene.astype(str)) == set(manifest['targets'])
    assert (response.target_gene.astype(str).value_counts() == 64).all()
    assert set(selected.loc[roles != 'perturbation', 'target_gene'].astype(str)) == {'non-targeting'}
    entries = [rt.json.loads(line) for line in Path('artifacts/HEK293T_charges.jsonl').read_text().splitlines()]
    bindings, charged = {}, 0
    for entry in entries:
        if entry['event'] == 'bind':
            assert entry['key'] not in bindings
            bindings[entry['key']] = entry['value']
        else:
            assert entry['event'] == 'charge' and entry['bytes'] > 0
            charged += entry['bytes']
            assert charged == entry['total'] <= 3_000_000_000
    assert bindings['limit_bytes'] == 3_000_000_000
    frozen = bindings['frozen_inputs']
    migration = bindings['metadata_decode_fix_v1']
    assert migration['prior_inputs'] == frozen
    current = migration['inputs']
    for key, file in [('selection', selection), ('genes', Path('artifacts/HEK293T_genes.parquet')),
                      ('offsets', Path('artifacts/HEK293T_offset_mapping.npz')),
                      ('implementation', Path('acquire_orion.py')), ('bridge', Path('durable_public_store.py'))]:
        assert current[key] == rt.sha(file)
        if key not in ['implementation', 'bridge']:
            assert frozen[key] == current[key]
    assert current['version'] == frozen['version'] == 908
    assert bindings['repository_revision'] == '598aa5442a4ee9a8f1a383a843ab5ce142904d04'
    files = {key: value for key, value in bindings.items() if key.startswith('file:')}
    assert files and all(value['etag'] and value['size'] > 0 for value in files.values())
    steps = [rt.json.loads(line) for line in Path('ledger/steps.jsonl').read_text().splitlines()]
    last = [entry for entry in steps if entry['name'] == 'acquire_orion_HEK293T'][-1]
    assert last['event'] == 'end' and last['rc'] == 0
    assert not Path('artifacts/HEK293T_payload_checkpoint').exists()
    fresh = rt.json.loads((folder/(NAME+'.json')).read_text())
    old = rt.json.loads(Path('artifacts/previous_job_downloads/cold_restore_receipts/tasks', NAME+'.json').read_text())
    transfer = rt.json.loads(Path('artifacts/HEK293T_raw_transfer.json').read_text())
    assert transfer['status'] == 'complete' and transfer['cells_done'] == len(selected)
    assert transfer['bytes'] == fresh['bytes'] == charged
    assert transfer['version'] == fresh['version'] == current['version'] == old['version']
    assert transfer['repository_revision'] == fresh['repository_revision'] == bindings['repository_revision']
    assert fresh['selection_sha256'] == current['selection'] == old['selection_sha256'] == manifest['selection_sha256']
    assert fresh['source'] == old['source']
    actual = rt.sha(folder/(NAME+'.h5ad'))
    archived = old['sha256'] if projected else old['unprojected_sha256']
    assert actual == fresh['sha256'] == archived, 'Changed bytes: investigate before using task'
    genes = pd.read_parquet('artifacts/HEK293T_genes.parquet').gene_id.astype(str).drop_duplicates().to_numpy()
    if projected:
        official = pd.read_csv(rt.ROOT.parent/'data/gene_names.csv').iloc[:, 0].astype(str)
        genes = official[official.isin(genes)].to_numpy()
        assert len(genes) == 18106 and fresh['unprojected_sha256'] == old['unprojected_sha256']
        assert fresh['readout_policy'] == 'measured official intersection'
    else:
        assert len(genes) == 38606 and 'readout_policy' not in fresh
    a = ad.read_h5ad(folder/(NAME+'.h5ad'), backed='r')
    total = 0
    try:
        assert a.n_obs == len(selected) and a.var_names.is_unique
        np.testing.assert_array_equal(a.var_names.to_numpy(), genes)
        np.testing.assert_array_equal(a.obs_names.to_numpy(), [f'HEK293T_{i}' for i in selected.cell_integer_id])
        for column in ['target_gene', 'batch', 'evaluation_role', 'source_row', 'cell_integer_id']:
            np.testing.assert_array_equal(a.obs[column].to_numpy(), selected[column].to_numpy())
        for lo in range(0, a.n_obs, 256):
            x = a.X[lo:lo+256]
            assert sparse.issparse(x) and np.issubdtype(x.dtype, np.integer)
            assert (x.data >= 0).all()
            totals = np.asarray(x.sum(axis=1)).ravel()
            original = selected.total_counts.iloc[lo:lo+256].to_numpy()
            assert (totals > 0).all()
            if projected:
                assert (totals <= original).all()
            else:
                np.testing.assert_array_equal(totals, original)
            total += int(totals.sum())
        shape = list(a.shape)
    finally:
        a.file.close()
    return dict(status='PASS', utc=datetime.now(timezone.utc).isoformat(), projected=projected,
        shape=shape, total_umis=total, sha256=actual, archived_sha256=archived,
        byte_identical_to_archived=True, frozen_metadata_members_verified=len(index['entries']),
        arms=roles.value_counts().to_dict(), charged_bytes=charged, ceiling_bytes=3_000_000_000,
        repository_revision=bindings['repository_revision'], immutable_files_bound=len(files),
        verifier_sha256=rt.sha(__file__), holdout_scored=False,
        limitation='Materialization identity/count audit only; pooled batch mismatch remains as frozen')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--projected', action='store_true')
    args = parser.parse_args()
    receipt = 'artifacts/HEK293T_materialization_' + ('projected' if args.projected else 'raw') + '_audit.json'
    try:
        result = verify(args.projected)
    except Exception as error:
        rt.save(receipt, dict(status='FAIL', error=f'{type(error).__name__}: {error}', holdout_scored=False))
        raise
    rt.save(receipt, result)
    print(result['status'], receipt, flush=True)

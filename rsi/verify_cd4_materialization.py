"""Audit restored CD4 counts and frozen design without evaluating a model.

--preflight checks metadata while acquisition is running. Full verification requires
the acquisition's successful step receipt. A changed task hash fails closed: the
source var table was not separately archived, so full measured-axis equality then
needs investigation rather than assuming an arbitrary official-gene subset is valid.
"""
import argparse
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

NAME = 'CD4_matched_search'
COLUMNS = ['target_gene', 'batch', 'guide_id', 'evaluation_role', 'source_row']


def audit_selection(selected, blueprint):
    assert selected.source_row.is_unique, 'Duplicate source cells'
    assert selected.source_row.is_monotonic_increasing, 'Frozen source order changed'
    roles = selected.evaluation_role.astype(str)
    assert set(roles) == {'perturbation', 'input_control', 'reference_control'}
    lanes = sorted(blueprint['audit']['target_batch_quotas'])
    response = selected.loc[roles == 'perturbation']
    assert set(response.target_gene.astype(str)) == set(blueprint['targets'])
    for target, group in response.groupby('target_gene', observed=True):
        assert len(group) == blueprint['n_real'], target
        np.testing.assert_array_equal(group.batch.astype(str).value_counts().reindex(lanes, fill_value=0),
            [blueprint['audit']['target_batch_quotas'][lane] for lane in lanes])
    for role, key in [('input_control', 'input_control_batch_quotas'),
                      ('reference_control', 'control_batch_quotas')]:
        rows = selected.loc[roles == role]
        assert set(rows.target_gene.astype(str)) == {'non-targeting'}
        expected = [blueprint['audit'][key][lane] for lane in lanes]
        assert len(rows) == sum(expected)
        np.testing.assert_array_equal(rows.batch.astype(str).value_counts().reindex(lanes, fill_value=0), expected)
    # Unique source rows above imply disjoint response/input/reference arms.
    return {role: int((roles == role).sum()) for role in sorted(set(roles))}


def audit_matrix(a, selected, expected_genes):
    assert a.n_obs == len(selected) and a.var_names.is_unique
    for column in COLUMNS:
        np.testing.assert_array_equal(a.obs[column].to_numpy(), selected[column].to_numpy())
    np.testing.assert_array_equal(a.var_names.to_numpy(), expected_genes)
    total_umis = 0
    for start in range(0, a.n_obs, 256):
        x = a.X[start:start+256]
        assert sparse.issparse(x), 'Expected sparse raw counts'
        values = x.data
        assert np.issubdtype(values.dtype, np.integer), 'Non-integer count dtype'
        assert np.all(values >= 0), 'Negative counts'
        totals = np.asarray(x.sum(axis=1)).ravel()
        assert np.all((totals > 0) & (totals <= 1_000_000)), 'Invalid cell depth'
        total_umis += int(totals.sum())
    return total_umis


def verify(preflight=False):
    folder = Path('artifacts/tasks')
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for file, expected in index['entries'].items():
        assert Path(file).stat().st_size == expected['bytes'] and rt.sha(file) == expected['sha256'], file
    blueprint = rt.json.loads((folder/'cd4_wide_blueprint.json').read_text())
    selection_path = folder/(NAME+'_selection.parquet')
    selected = pd.read_parquet(selection_path)
    selection_hash = rt.sha(selection_path)
    old = rt.json.loads(Path('artifacts/previous_job_downloads/cold_restore_receipts/tasks', NAME+'.json').read_text())
    assert selection_hash == blueprint['selection_sha256'] == old['selection_sha256']
    assert len(blueprint['targets']) == 300 and blueprint['n_real'] == 96
    quotas = blueprint['audit']['target_batch_quotas']
    assert len(quotas) == 24 and set(quotas.values()) == {4}
    arms = audit_selection(selected, blueprint)
    assert arms == dict(input_control=1920, perturbation=28800, reference_control=768)
    census = pd.read_parquet('artifacts/CD4_cell_census.parquet').set_index('source_row')
    assert census.index.is_unique
    rows = census.loc[selected.source_row]
    for column in ['target_gene', 'batch', 'guide_id', 'lane_id', 'guide_type', 'total_counts', 'low_quality']:
        np.testing.assert_array_equal(selected[column].to_numpy(), rows[column].to_numpy())
    assert not rows.low_quality.astype(bool).any() and (rows.total_counts > 0).all()
    del census, rows
    entries = [rt.json.loads(line) for line in Path('artifacts/CD4_destination_charges.jsonl').read_text().splitlines()]
    bindings = {e['key']: e['value'] for e in entries if e['event'] == 'bind'}
    frozen = bindings['frozen_inputs']
    assert frozen['blueprint'] == rt.sha(folder/'cd4_wide_blueprint.json')
    assert frozen['selection'] == selection_hash
    assert frozen['census'] == rt.sha('artifacts/CD4_cell_census.parquet') == old['source_census_sha256']
    assert frozen['implementation'] == rt.sha('acquire_cd4_destination.py')
    source = bindings['source_identity']
    feasibility = rt.json.loads(Path('artifacts/CD4_cell_feasibility.json').read_text())
    assert source['url'] == frozen['source'] == old['source'] == feasibility['source']
    assert source['etag'] == old['etag'] == feasibility['etag']
    assert source['size'] == feasibility['source_size']
    assert bindings['limit_bytes'] == 2_600_000_000
    plan = bindings['selection_and_plan']
    assert plan['inputs'] == frozen and plan['source'] == source
    import hashlib
    assert plan['rows_sha256'] == hashlib.sha256(selected.source_row.to_numpy(np.int64).tobytes()).hexdigest()
    result = dict(status='PASS metadata preflight; live counts UNVERIFIED',
        utc=datetime.now(timezone.utc).isoformat(), frozen_metadata_members_verified=len(index['entries']),
        arms=arms, lanes=24, targets=300, selection_sha256=selection_hash,
        source=source, verifier_sha256=rt.sha(__file__), holdout_scored=False)
    if preflight:
        return result
    step_rows = [rt.json.loads(line) for line in Path('ledger/steps.jsonl').read_text().splitlines()]
    last = [e for e in step_rows if e['name'] == 'acquire_cd4_destination'][-1]
    assert last['event'] == 'end' and last['rc'] == 0, 'Acquisition must end successfully before audit'
    assert not Path('artifacts/CD4_destination_payload_checkpoint').exists()
    fresh = rt.json.loads((folder/(NAME+'.json')).read_text())
    actual_hash = rt.sha(folder/(NAME+'.h5ad'))
    assert actual_hash == fresh['sha256'], 'Fresh task differs from its receipt'
    assert actual_hash == old['sha256'], 'Changed bytes: investigate measured axis/count identity before calibration'
    assert fresh['selection_sha256'] == selection_hash
    assert fresh['source_census_sha256'] == frozen['census']
    assert fresh['source'] == source['url'] and fresh['etag'] == source['etag'] and fresh['source_size'] == source['size']
    charged = sum(e['bytes'] for e in entries if e['event'] == 'charge')
    transfer = rt.json.loads(Path('artifacts/CD4_destination_transfer.json').read_text())
    assert charged == fresh['bytes'] == transfer['bytes'] <= bindings['limit_bytes']
    assert transfer['etag'] == source['etag'] and transfer['source_size'] == source['size']
    a = ad.read_h5ad(folder/(NAME+'.h5ad'), backed='r')
    try:
        official = pd.read_csv(rt.ROOT.parent/'data/gene_names.csv').iloc[:, 0].astype(str)
        expected_genes = official[official.isin(a.var_names)].to_numpy()
        assert len(expected_genes) == 17772, 'Archived measured official intersection size'
        total = audit_matrix(a, selected, expected_genes)
        result.update(status='PASS', shape=list(a.shape), total_umis=total)
    finally:
        a.file.close()
    result.update(sha256=actual_hash, archived_sha256=old['sha256'], byte_identical_to_archived=True,
        charged_bytes=charged, ceiling_bytes=bindings['limit_bytes'],
        axis_evidence='Official order and 17,772 readouts checked directly; complete measured-axis identity from byte-exact archive match',
        limitation='Materialization identity/count audit only; no model or calibration claim')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--preflight', action='store_true')
    args = parser.parse_args()
    receipt = 'artifacts/CD4_materialization_preflight.json' if args.preflight else 'artifacts/CD4_materialization_audit.json'
    try:
        result = verify(args.preflight)
    except Exception as error:
        rt.save(receipt, dict(status='FAIL', utc=datetime.now(timezone.utc).isoformat(),
            error=f'{type(error).__name__}: {error}', verifier_sha256=rt.sha(__file__), holdout_scored=False))
        raise
    rt.save(receipt, result)
    print(result['status'], receipt, flush=True)

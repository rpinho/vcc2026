"""Independent frozen-row/axis/count audit; never invokes a HOLDOUT scorer."""
import runtime as rt
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

folder = Path('artifacts/tasks')
archived = Path('artifacts/previous_job_downloads/cold_restore_receipts/tasks')
blueprint = rt.json.loads((folder/'h1_rich_blueprint.json').read_text())
index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
for file, expected in index['entries'].items():
    assert Path(file).stat().st_size == expected['bytes'] and rt.sha(file) == expected['sha256'], file
official = pd.read_csv(rt.ROOT.parent/'data/gene_names.csv').iloc[:, 0].astype(str)
measured = set(pd.read_parquet('artifacts/H1_full_var.parquet').index.astype(str))
expected_genes = official[official.isin(measured)].to_numpy()
tasks, controls = {}, []
for name in ['H1_rich_search', 'H1_rich_holdout']:
    m = blueprint[name]
    selection_path = folder/(name+'_selection.parquet')
    selected = pd.read_parquet(selection_path)
    assert selected.source_row.is_unique
    old = rt.json.loads((archived/(name+'.json')).read_text())
    fresh = rt.json.loads((folder/(name+'.json')).read_text())
    actual_hash = rt.sha(folder/(name+'.h5ad'))
    assert fresh['sha256'] == actual_hash
    assert fresh['selected_rows_sha256'] == rt.sha(selection_path) == old['selected_rows_sha256']
    assert fresh['source_generation'] == old['source_generation'] == '1765904883947296'
    a = ad.read_h5ad(folder/(name+'.h5ad'), backed='r')
    try:
        for column in ['target_gene', 'batch', 'guide_id', 'evaluation_role', 'source_row']:
            np.testing.assert_array_equal(a.obs[column].to_numpy(), selected[column].to_numpy())
        np.testing.assert_array_equal(a.var_names.to_numpy(), expected_genes)
        roles = selected.evaluation_role.astype(str)
        response = selected.loc[roles=='perturbation']
        assert set(response.target_gene.astype(str)) == set(m['targets'])
        assert (response.groupby('target_gene', observed=True).size() == 400).all()
        batches = sorted(selected.batch.astype(str).unique())
        assert len(batches)==48
        quotas = m['audit']['target_batch_quotas']
        for _, group in response.groupby('target_gene', observed=True):
            np.testing.assert_array_equal(group.batch.astype(str).value_counts().reindex(batches, fill_value=0), [quotas[b] for b in batches])
        for role in ['input_control', 'reference_control']:
            rows = selected.loc[roles==role]
            assert len(rows)==3200 and set(rows.target_gene.astype(str))=={'non-targeting'}
            np.testing.assert_allclose(rows.batch.astype(str).value_counts().reindex(batches, fill_value=0).to_numpy()/3200,
                np.array([quotas[b] for b in batches])/400, atol=0, rtol=0)
            controls.append(set(rows.source_row))
        total_umis = 0
        for start in range(0, a.n_obs, 256):
            x = a.X[start:start+256]
            assert np.issubdtype(x.dtype, np.integer) and np.all(x.data>=0)
            totals = np.asarray(x.sum(axis=1)).ravel()
            assert np.all(totals>0)
            total_umis += int(totals.sum())
        tasks[name] = dict(shape=list(a.shape), targets=len(m['targets']), real_cells_per_target=400,
            input_controls=3200, reference_controls=3200, batches=48, sha256=actual_hash,
            archived_sha256=old['sha256'], byte_identical_to_archived=actual_hash==old['sha256'],
            selection_sha256=rt.sha(selection_path), total_umis=total_umis,
            checks='Frozen observation order, official measured axis, integer nonnegative real-cell counts, exact batch quotas')
    finally:
        a.file.close()
assert all(a.isdisjoint(b) for i,a in enumerate(controls) for b in controls[i+1:])
assert set(blueprint['H1_rich_search']['targets']).isdisjoint(blueprint['H1_rich_holdout']['targets'])
charges = [rt.json.loads(line) for line in Path('artifacts/H1_rich_charges.jsonl').read_text().splitlines()]
bindings = {row['key']: row['value'] for row in charges if row['event']=='bind'}
total = sum(row['bytes'] for row in charges if row['event']=='charge')
receipt = rt.json.loads(Path('artifacts/H1_rich_acquisition.json').read_text())
assert receipt['bytes']==total<=5_500_000_000
assert receipt['source_etag']==bindings['source_identity']['etag']=='"-CKCuwqrMwpEDEAE="'
assert receipt['source_size']==bindings['source_identity']['size']==15_482_497_461
assert not Path('artifacts/H1_rich_payload_checkpoint').exists()
rt.save('artifacts/live_restore_round6_h1.json', dict(status='PASS', utc=datetime.now(timezone.utc).isoformat(),
    round=6, frozen_metadata_members_verified=len(index['entries']), tasks=tasks,
    all_four_control_arms_disjoint=True, charged_bytes=total, ceiling_bytes=5_500_000_000,
    source_etag=receipt['source_etag'], source_size=receipt['source_size'], holdout_scored=False,
    verifier_sha256=rt.sha(__file__), limitation='Identity and count audit only; no calibration or biological score claim'))
print('PASS independent H1 rich materialization audit', tasks)

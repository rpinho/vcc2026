"""Count-only feasibility certificates for existing, unscored native contexts.

Reads observation metadata, never expression matrices. This selects no panel and
changes no protocol. Coordinate bounds apply even if some batches are omitted.
"""
from collections import Counter
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
import hashlib
import json

import h5py
import numpy as np
import pandas as pd

import runtime as rt


def read_column(obs, name):
    column = obs[name]
    if isinstance(column, h5py.Group):
        categories = column['categories'].asstr()[:]
        codes = column['codes'][:]
    elif '__categories' in obs and name in obs['__categories']:
        categories = obs['__categories'][name].asstr()[:]
        codes = column[:]
    else:
        return column.asstr()[:] if column.dtype.kind in 'OS' else column[:]
    assert np.all(codes >= 0), (name, 'missing categorical label')
    return categories[codes]


def upper_bound(counts, panel_size):
    if len(counts) < panel_size:
        return None
    return np.sort(counts, axis=0)[-panel_size].astype(np.int64)


def check_bound_logic():
    # Exhaust every subset in a small asymmetric example, including zero quotas.
    counts = np.array([[8, 0, 4], [6, 4, 3], [2, 7, 5], [0, 8, 2], [3, 3, 3]])
    cases = 0
    for k in range(1, len(counts)+1):
        ceiling = upper_bound(counts, k)
        for subset in combinations(range(len(counts)), k):
            quota = counts[list(subset)].min(axis=0)
            assert np.all(quota <= ceiling)
            assert int(quota.sum()) <= int(ceiling.sum())
            cases += 1
    assert upper_bound(counts, 6) is None
    return cases


def inspect(context, frame, path, source_sets):
    frame = frame.astype(str)
    assert frame.row_id.is_unique
    fingerprint = hashlib.sha256(frame.to_csv(index=False).encode()).hexdigest()
    # Separate implementations: pandas cross-tabulation and observation counters.
    table = pd.crosstab(frame.target, frame.batch).sort_index().sort_index(axis=1)
    counter = Counter(zip(frame.target, frame.batch))
    assert int(table.to_numpy().sum()) == len(frame)
    for target, row in table.iterrows():
        assert row.tolist() == [counter[(target, batch)] for batch in table.columns]
    controls = table.loc['non-targeting']
    table = table.drop(index='non-targeting')
    target_path = Path('artifacts') / (context+'_count_only_capacity.parquet')
    table.to_parquet(target_path)
    results = {}
    for scope, available in dict(all_observed=set(table.index), **source_sets).items():
        tab = table.loc[sorted(set(table.index) & available)]
        counts = tab.to_numpy(np.int64)
        bounds = {}
        for k in [20, 79, 100, 200, 300]:
            ceiling = upper_bound(counts, k)
            if ceiling is None:
                bounds[str(k)] = dict(status='insufficient targets')
                continue
            # Independent scalar kth-largest verification, with no numpy sort.
            independent = [sorted((counter[(t, b)] for t in tab.index), reverse=True)[k-1]
                           for b in table.columns]
            assert ceiling.tolist() == independent
            bounds[str(k)] = dict(
                unlimited_control_depth_upper_bound=int(ceiling.sum()),
                per_batch_upper_bounds=dict(zip(table.columns, map(int, ceiling))),
                impossible_at_400_cells=bool(ceiling.sum() < 400))
        totals = counts.sum(axis=1)
        results[scope] = dict(
            candidate_targets=len(tab), actual_cells_ge400_targets=int((totals >= 400).sum()),
            actual_cells_ge800_targets=int((totals >= 800).sum()),
            common_mixture_bounds=bounds,
            uniform_quota_capacities=[dict(quota_per_batch=q, cells=q*len(table.columns),
                targets=int((counts >= q).all(axis=1).sum())) for q in [1, 2, 4, 8, 10]])
    stat = path.stat()
    return dict(
        path=str(path), raw_file_size=stat.st_size, raw_file_mtime_ns=stat.st_mtime_ns,
        full_raw_file_sha256='NOT COMPUTED; expression payload not read',
        selected_metadata_csv_sha256=fingerprint, selected_metadata_columns=list(frame.columns),
        observations=len(frame), perturbation_targets=len(table), batches=len(table.columns),
        controls_per_batch={str(k): int(v) for k, v in controls.items()},
        count_table_path=str(target_path), count_table_sha256=rt.sha(target_path), scopes=results)


def main():
    fixtures = check_bound_logic()
    index = json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path, receipt in index['entries'].items():
        assert rt.sha(path) == receipt['sha256'], path
    assert rt.sha(index['archive']) == index['archive_sha256']
    # This census covers the full raw GWPS metadata, beyond already built banks.
    source_path = Path('artifacts/GWPS_raw_obs.parquet')
    gwps = pd.read_parquet(source_path, columns=['gene'])
    n = gwps.groupby('gene', observed=True).size()
    source_sets = {'full_raw_GWPS_ge200': set(n[n >= 200].index.astype(str))}
    bank_paths = {
        'prepared_base_GWPS': Path('artifacts/sources/K562_GWPS_bulk.npz'),
        'prepared_CD4_GWPS_raw': Path('artifacts/cd4_banks/K562_GWPS_raw.npz'),
        'historical_v4_GWPS_raw': rt.ROOT.parent/'public_data/model_v4_signature_banks/K562_GWPS_raw.npz'}
    for name, path in bank_paths.items():
        with np.load(path) as bank:
            source_sets[name] = set(bank['targets'].astype(str))
    output = {}
    for context in ['HepG2', 'Jurkat']:
        path = rt.ROOT.parent/'public_data'/(context+'.h5ad')
        with h5py.File(path) as f:
            obs = f['obs']
            frame = pd.DataFrame(dict(row_id=read_column(obs, str(obs.attrs['_index'])),
                target=read_column(obs, 'gene'), batch=read_column(obs, 'gem_group')))
        output[context] = inspect(context, frame, path, source_sets)
    path = rt.ROOT.parent/'model_v4/artifacts/HCT116_cells.parquet'
    obs = pd.read_parquet(path, columns=['cell_integer_id', 'gene_target', 'sample', 'pass_guide_filter'])
    obs = obs.loc[obs.pass_guide_filter].copy()
    frame = obs.rename(columns=dict(cell_integer_id='row_id', gene_target='target', sample='batch'))
    frame['target'] = frame.target.replace({'Non-Targeting': 'non-targeting'})
    output['HCT116_cached_subset'] = inspect('HCT116_cached_subset', frame[['row_id', 'target', 'batch']], path, source_sets)
    for name, record in output.items():
        all_targets = record['scopes']['all_observed']
        print(name, '400-cell targets', all_targets['actual_cells_ge400_targets'],
              '300-target common-depth upper bound', all_targets['common_mixture_bounds']['300'], flush=True)
    rt.save('artifacts/native_capacity_audit.json', dict(
        status='PASS independently checked metadata counts and exact integer bounds',
        utc=datetime.now(timezone.utc).isoformat(), verifier_sha256=rt.sha(__file__),
        frozen_metadata_members_verified=len(index['entries']), exhaustive_subset_fixture_cases=fixtures,
        gwps_census_sha256=rt.sha(source_path), gwps_ge200_targets=len(source_sets['full_raw_GWPS_ge200']),
        bank_hashes={name: rt.sha(path) for name, path in bank_paths.items()}, contexts=output,
        proof='For any k-target subset S and common nonnegative integer batch quotas q, '
              'q[b] <= min_{t in S} counts[t,b] <= kth-largest_t counts[t,b]. '
              'Sum these coordinate ceilings. Zero quotas permit arbitrary omitted batches. '
              'Ignoring control supply only enlarges the bound. Separately, fewer than k targets '
              'with total cells >=400 excludes even target-specific matching at depth400.',
        expression_payloads_read=False, holdout_scored=False, new_task_registered=False,
        limitation='Availability certificates concern existing cached subsets, not full public releases. '
                   'A depth ceiling below400 does not prove a lower-depth official scale invalid. '
                   'Full raw GWPS metadata availability does not mean a materialized bank exists. '
                   'Raw H5AD byte identities are not asserted: only the selected metadata are fingerprinted.'))
    print('PASS native capacity audit; no panel selected, expressions read or protocol changed', flush=True)


if __name__ == '__main__':
    main()

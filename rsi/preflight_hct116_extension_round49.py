"""Freeze a metadata-only source extension; no expression/network/model access."""
import argparse
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt
import lance
import numpy as np
import pandas as pd

ROOT = Path('artifacts/HCT116_extension_preflight_round49')
PLAN = ROOT/'plan_retry.json'
META = Path('artifacts/HCT116_full_cells_metadata_round48')
COLS = ['cell_integer_id', 'gene_target', 'sample', 'pass_guide_filter',
        'cell_start_index', 'n_genes_by_counts', 'total_counts']
BANK = Path('artifacts/cd4_banks/HCT116.npz')
OFFSETS = Path('artifacts/HCT116_offset_mapping.npz')
FRAGS = Path('../model_v4/artifacts/HCT116_fragments.json')


def load(p):
    return rt.json.loads(Path(p).read_text())


def save_new(p, obj):
    data = rt.json.dumps(obj, indent=2, allow_nan=False)+'\n'
    with Path(p).open('x') as f:
        f.write(data)


def freeze():
    ROOT.mkdir(exist_ok=True)
    assert not PLAN.exists() and not (ROOT/'preflight.json').exists()
    paths = [Path(__file__), Path('audit_hct116_extension_round49.py'),
        META/'receipt.json', META/'audit_retry.json', META/'audit_plan_retry.json',
        Path('artifacts/tasks/manifest.json'), Path('artifacts/benchmark_roles.json'),
        Path('artifacts/calibration_gate.json'), Path('artifacts/download_allocations.json'),
        Path('artifacts/frozen_metadata_index.json'), BANK, OFFSETS, FRAGS,
        Path('../model_v4/artifacts/HCT116_fragment_boundaries.npz'),
        Path('../model_v4/artifacts/HCT116_offset_verification.json'),
        Path('../model_v4/artifacts/HCT116_cells.parquet'),
        Path('../model_v4/artifacts/HCT116_genes.parquet'),
        Path('../model_v4/prepare_orion.py'), Path('../model_v4/correct_orion_controls.py'),
        Path('vendor/v4_model.py'), Path('vendor/atlas_model.py'), Path('recipes.py')]
    paths += [Path(e['path']) for e in load(META/'receipt.json')['files'].values()]
    paths += sorted(Path('artifacts/scores').glob('*/*_s42.json'))
    save_new(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        scope='Existing CD4_depth192_search79 targets only. Inspect all their full HCT116 metadata, keep all guide-pass cells in historical present expression blocks. Append only targets with original bank n<20 and full filtered n>=20. HMGXB4 stays on the destination panel and abstains below20; do not acquire its two cells. Controls: original first250 by cell ID per sample BEFORE block exclusion; no replacement. No source activation, expression reads or biological scoring.',
        columns=COLS, minimum_source_cells=20, controls_per_sample=250,
        raw_moment_contract='Reuse original prepare_orion/correct_orion_controls: count sums and mapped library depths; per-cell proportions divided by FULL metadata total_counts; squared proportions; sample count weights for CPM and library weights for bulk; original conservative mean variances and float32 publication. Exact original readout axis and old bank rows retained. New namespace only.',
        centering_contract='Any future union bank retains437 historical target rows plus eligible appended targets; raw old rows byte-identical. Source-common subtraction may change when its eligible population grows. That is a new recipe/source identity and requires a frozen separate Phase1 revision, never historical score reassignment. Do not recenter over just the destination79.',
        identity_limit='New cell metadata is pinned at598aa5442a4ee9a8f1a383a843ab5ce142904d04; historical expression version682 offset/fragment map is checked for internal compatibility only, NOT freshly expression-certified. File-local logical record ranges are exact conditional on this map. Compressed HTTP byte ranges/cost require pinned expression manifest/file metadata; bytes per record is not a transfer estimate.',
        io_policy='No public I/O or reallocation in this preflight. No expression payload until fresh expression identity, exact bounded request plan and separately preallocated durable reader/checkpoint/failure audit. Keep all historical banks and task/control/source/seed/score/HOLDOUT contracts.',
        sha256={str(p):rt.sha(p) for p in paths}))
    print('FROZEN', PLAN, flush=True)


def main():
    assert not (ROOT/'preflight.json').exists(), 'Preserve completed evidence'
    plan = load(PLAN)
    for p, h in plan['sha256'].items():
        assert rt.sha(p) == h, p
    assert load(META/'audit_retry.json')['status'].startswith('PASS')
    manifest = load('artifacts/tasks/manifest.json')
    panel = manifest['CD4_depth192_search']['targets']
    assert len(panel) == len(set(panel)) == 79
    ds = lance.dataset(str((META/'release/cells.lance').resolve()))
    cells = ds.to_table(columns=COLS).to_pandas()
    ids = cells.cell_integer_id.to_numpy(np.int64)
    nnz = cells.n_genes_by_counts.to_numpy(np.int64)
    starts = cells.cell_start_index.to_numpy(np.int64)
    assert np.array_equal(ids, np.arange(len(cells)))
    assert np.all(nnz > 0) and starts[0] == 0
    assert np.array_equal(starts[1:], starts[:-1]+nnz[:-1])
    assert np.isfinite(cells.total_counts).all() and (cells.total_counts > 0).all()
    with np.load(OFFSETS, allow_pickle=False) as d:
        off = {k:d[k] for k in d.files}
    assert int(off['version']) == 682
    np.testing.assert_array_equal(starts[::5000], off['original_start'])
    block_lengths = np.add.reduceat(nnz, np.arange(0,len(cells),5000))
    present_blocks = off['logical_start'] >= 0
    np.testing.assert_array_equal(block_lengths[present_blocks], off['block_length'][present_blocks])
    assert np.flatnonzero(~present_blocks).tolist() == [374]
    block = ids//5000
    present = present_blocks[block]
    eligible = cells.loc[cells.pass_guide_filter & present].copy()
    with np.load(BANK, allow_pickle=False) as d:
        old = dict(zip(d['targets'].astype(str), map(int,d['n_cells'])))
        genes = d['genes'].astype(str)
    assert len(old) == 437
    original = pd.read_parquet('../model_v4/artifacts/HCT116_cells.parquet', columns=COLS)
    replay = cells.iloc[original.cell_integer_id.to_numpy(np.int64)].reset_index(drop=True)
    pd.testing.assert_frame_equal(replay[COLS].astype({'gene_target':str,'sample':str}),
        original[COLS].reset_index(drop=True).astype({'gene_target':str,'sample':str}), check_dtype=False)
    all_counts = eligible.gene_target.value_counts()
    for t,n in old.items():
        assert all_counts.get(t,0) == n, (t,n,all_counts.get(t,0))
    raw_counts = cells.loc[cells.pass_guide_filter].gene_target.value_counts()
    target_rows = []
    for t in panel:
        n = int(all_counts.get(t,0)); prior = old.get(t,0)
        action = 'retain_existing' if prior >= 20 else ('append' if n >= 20 else 'abstain_below20')
        selected = eligible.loc[eligible.gene_target == t]
        target_rows.append(dict(target=t, old_bank_cells=prior, full_guide_cells=int(raw_counts.get(t,0)),
            block_present_cells=n, source_samples=int(selected['sample'].nunique()),
            expression_records=int(selected.n_genes_by_counts.sum()), action=action))
    table = pd.DataFrame(target_rows)
    additions = table.loc[table.action == 'append','target'].tolist()
    assert not set(additions) & set(manifest['H1_holdout']['targets'])
    ntc_before = cells.loc[cells.pass_guide_filter & (cells.gene_target == 'Non-Targeting')].sort_values(
        ['sample','cell_integer_id']).groupby('sample',observed=True).head(250)
    ntc = ntc_before.loc[present[ntc_before.index]].copy()
    old_ntc = original.loc[original.pass_guide_filter & (original.gene_target == 'Non-Targeting')].sort_values(
        ['sample','cell_integer_id']).groupby('sample',observed=True).head(250)
    old_ntc = old_ntc.loc[present_blocks[old_ntc.cell_integer_id.to_numpy(np.int64)//5000]]
    assert sorted(ntc.cell_integer_id) == sorted(old_ntc.cell_integer_id)
    pert = eligible.loc[eligible.gene_target.isin(additions)].copy()
    assert set(pert['sample']) <= set(ntc['sample'])
    pert['source_role']='append_target'; ntc['source_role']='source_control'
    selected = pd.concat([pert,ntc],ignore_index=True).sort_values('cell_integer_id').reset_index(drop=True)
    assert selected.cell_integer_id.is_unique
    assert (table.action == 'abstain_below20').sum() == 1
    assert table.loc[table.action == 'abstain_below20','target'].tolist() == ['HMGXB4']
    assert int(all_counts['HMGXB4']) == 2
    selected.to_parquet(ROOT/'selection.parquet',index=False)
    table.to_csv(ROOT/'targets.tsv',sep='\t',index=False)
    pd.crosstab(pert.gene_target,pert['sample']).reindex(index=additions,
        columns=sorted(ntc['sample'].unique()),fill_value=0).to_csv(ROOT/'target_sample_counts.tsv',sep='\t')
    controls = ntc.groupby('sample',observed=True).agg(cells=('cell_integer_id','size'),
        expression_records=('n_genes_by_counts','sum'), metadata_total_counts=('total_counts','sum'))
    controls.to_csv(ROOT/'controls.tsv',sep='\t')

    frags = load(FRAGS)
    with np.load('../model_v4/artifacts/HCT116_fragment_boundaries.npz') as d:
        np.testing.assert_array_equal([f['id'] for f in frags],d['fragment_ids'])
        np.testing.assert_array_equal([f['physical_rows'] for f in frags],d['physical_rows'])
    boundaries = np.r_[0,np.cumsum([f['physical_rows'] for f in frags],dtype=np.int64)]
    assert boundaries[-1] == off['block_length'].sum()
    ranges = []
    for row in selected.itertuples(index=False):
        b = int(row.cell_integer_id)//5000
        start = int(row.cell_start_index-off['original_start'][b]+off['logical_start'][b])
        stop = start+int(row.n_genes_by_counts)
        assert off['logical_start'][b] <= start < stop <= off['logical_start'][b]+off['block_length'][b]
        while start < stop:
            fi = int(np.searchsorted(boundaries,start,side='right')-1)
            end = min(stop,int(boundaries[fi+1]))
            ranges.append(dict(cell_integer_id=int(row.cell_integer_id),fragment_id=int(frags[fi]['id']),
                path=frags[fi]['path'], record_start=int(start-boundaries[fi]),record_stop=int(end-boundaries[fi])))
            start=end
    ranges = pd.DataFrame(ranges).sort_values(['fragment_id','record_start']).reset_index(drop=True)
    for _,rows in ranges.groupby('fragment_id',sort=False):
        assert (rows.record_start.to_numpy()[1:] >= rows.record_stop.to_numpy()[:-1]).all()
    assert int((ranges.record_stop-ranges.record_start).sum()) == int(selected.n_genes_by_counts.sum())
    ranges.to_parquet(ROOT/'logical_record_ranges.parquet',index=False)
    artifacts = ['selection.parquet','targets.tsv','target_sample_counts.tsv','controls.tsv','logical_record_ranges.parquet']
    save_new(ROOT/'preflight.json',dict(status='COUNT/RANGE PREFLIGHT COMPLETE; independent audit pending; expression cost UNVERIFIED',
        utc=datetime.now(timezone.utc).isoformat(),plan_sha256=rt.sha(PLAN),
        full_metadata_cells=len(cells),expression_version=682,expression_identity_freshly_verified=False,
        panel_targets=79,retained_supported_targets=int((table.action=='retain_existing').sum()),
        append_targets=additions,abstaining_targets=['HMGXB4'],destination_targets_changed=False,
        target_cells=len(pert),control_cells=len(ntc),control_samples=len(controls),
        target_expression_records=int(pert.n_genes_by_counts.sum()),
        control_expression_records=int(ntc.n_genes_by_counts.sum()),
        total_expression_records=int(selected.n_genes_by_counts.sum()),
        touched_historical_fragments=int(ranges.fragment_id.nunique()),logical_spans=len(ranges),
        gene_axis=len(genes),all_original_bank_metadata_counts_replayed=True,
        all_original_cell_metadata_columns_replayed=True,control_selection_identical=True,
        omitted_control_cells_without_expression=len(ntc_before)-len(ntc),
        exact_http_payload_bytes=None,affordability='UNVERIFIED pending pinned expression manifest and encoded range plan',
        why_no_byte_claim='Historical fragment rows and cell n_genes_by_counts give logical record ranges only; encoded page/miniblock offsets and HTTP coalescing/retry costs are absent. No bytes-per-record extrapolation or full-file-download impossibility claim.',
        future_union_targets=len(old)+len(additions),minimum_source_cells=20,
        artifact_sha256={p:rt.sha(ROOT/p) for p in artifacts},
        expression_payloads_read=False,network_io=False,source_activated=False,
        model_scores_read=False,holdout_scored=False,phase2_allowed=False))
    for p,h in plan['sha256'].items():
        assert rt.sha(p)==h,p
    print('PREFLIGHT', len(additions),'new targets',len(pert),'target cells',len(ntc),'controls',
        int(selected.n_genes_by_counts.sum()),'records; transfer bytes UNVERIFIED',flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true')
    args=parser.parse_args()
    freeze() if args.freeze else main()

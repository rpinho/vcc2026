"""Offline scope certificate: cached metadata only; never import original probes."""
import argparse
import ast
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import runtime as rt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

PLAN = Path('artifacts/HCT116_census_scope_plan_round48.json')
OUT = Path('artifacts/HCT116_census_scope_audit_round48.json')
TABLE = Path('artifacts/HCT116_census_scope_targets_round48.tsv')
PARENT = Path('../model_v4')
CELLS = PARENT/'artifacts/HCT116_cells.parquet'
RAW = Path('artifacts/GWPS_raw_obs.parquet')
CD4 = Path('artifacts/CD4_cell_census.parquet')


def read(p):
    return json.loads(Path(p).read_text())


def save_new(p, data):
    with p.open('x') as f:
        f.write(json.dumps(data, indent=2, allow_nan=False) + '\n')


def check_bindings(plan):
    for p, digest in plan['sha256'].items():
        assert rt.sha(p) == digest, p


def freeze():
    assert not PLAN.exists() and not OUT.exists() and not TABLE.exists()
    paths = [Path(__file__), CELLS, RAW, CD4, PARENT/'probe.py',
        PARENT/'prepare_orion.py', PARENT/'orion_offsets.py',
        PARENT/'artifacts/targets.json', PARENT/'artifacts/h1_targets.csv',
        Path('../data/pert_counts.csv'), PARENT/'artifacts/hct_tree.json',
        PARENT/'artifacts/hct_top.json', PARENT/'artifacts/orion_readme.md',
        PARENT/'artifacts/HCT116_offset_mapping.npz',
        PARENT/'artifacts/HCT116_offset_verification.json',
        Path('artifacts/HCT116_matching_audit_round40.json'),
        Path('artifacts/CD4_source_capacity_audit_round44.json'),
        Path('artifacts/calibration_gate.json'), Path('artifacts/benchmark_roles.json'),
        Path('artifacts/tasks/manifest.json'), Path('artifacts/download_allocations.json'),
        Path('artifacts/frozen_metadata_index.json')]
    paths += sorted(Path('artifacts/scores').glob('*/*_s42.json'))
    save_new(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        scope='Cached HCT116 query provenance and metadata population only. No response values, expression payload, model metrics, task activation or network I/O.',
        columns=dict(hct=['cell_integer_id','gene_target','sample','pass_guide_filter'],
                     gwps=['gene'], cd4=['target_gene']),
        estimands='Replay historical437-target query labels; independent cache row/count certification; excluded target-name universe relative to full cached GWPS>=200 and CD4 labels; expression-block presence only; assess whether full metadata may add target evidence. No target or lane selection.',
        identity_limits='Freeze present bytes. Original probe queried unpinned HF HEAD and saved no cell-dataset version or immutable commit in the parquet. Code consistency cannot prove original query completeness or exact release identity. Expression offset version does not identify cells metadata version.',
        sha256={str(p):rt.sha(p) for p in paths}))
    print('FROZEN', PLAN, flush=True)


def streamed(path, columns):
    counter = Counter()
    for batch in pq.ParquetFile(path).iter_batches(batch_size=65536, columns=columns):
        counter.update(zip(*(batch.column(i).to_pylist() for i in range(len(columns)))))
    return counter


def main():
    assert not OUT.exists() and not TABLE.exists(), 'Preserve existing evidence'
    plan = read(PLAN)
    check_bindings(plan)
    source = (PARENT/'probe.py').read_text()
    tree = ast.parse(source)
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Attribute) and n.func.attr == 'to_table']
    query = next(n for n in calls if any(k.arg == 'filter' for k in n.keywords))
    kwargs = {k.arg:k.value for k in query.keywords}
    assert ast.literal_eval(kwargs['columns']) == ['cell_integer_id','cell_start_index',
        'n_genes_by_counts','total_counts','sample','gene_target','pass_guide_filter']
    assert ast.unparse(kwargs['filter']) == "f'gene_target IN ({labels})'"
    labels = next(n.value for n in ast.walk(tree) if isinstance(n, ast.Assign)
                  and any(isinstance(t, ast.Name) and t.id == 'labels' for t in n.targets))
    assert isinstance(labels, ast.Call) and isinstance(labels.args[0], ast.GeneratorExp)
    assert ast.unparse(labels.args[0].generators[0].iter) == "targets + ['Non-Targeting']"
    assert "version=" not in source and "pass_guide_filter==" not in source
    assert "cells.to_parquet(rt/f'artifacts/{name}_cells.parquet')" in source
    public = pd.read_csv('../data/pert_counts.csv').target_gene.tolist()
    h1 = pd.read_csv(PARENT/'artifacts/h1_targets.csv').target_gene.tolist()
    replay = [x for x in dict.fromkeys(public+h1) if x != 'non-targeting']
    declared = read(PARENT/'artifacts/targets.json')
    assert replay == declared and len(declared) == len(set(declared)) == 437
    queried = set(declared) | {'Non-Targeting'}
    c = pd.read_parquet(CELLS, columns=plan['columns']['hct'])
    assert c.cell_integer_id.is_unique and not c.isna().any().any()
    assert set(c.gene_target) <= queried
    counts = Counter(zip(c.gene_target, c['sample'], c.pass_guide_filter))
    independent = streamed(CELLS, ['gene_target','sample','pass_guide_filter'])
    assert counts == independent and sum(counts.values()) == len(c)
    by_target = c.groupby('gene_target', observed=True).size().to_dict()
    by_sample = c.groupby('sample', observed=True).size().to_dict()
    assert all(n == sum(v for (t,b,g),v in independent.items() if t == target)
               for target,n in by_target.items())
    assert all(n == sum(v for (t,b,g),v in independent.items() if b == sample)
               for sample,n in by_sample.items())

    gwps = pd.read_parquet(RAW, columns=['gene']).gene.value_counts()
    gwps_counter = streamed(RAW, ['gene'])
    assert gwps.to_dict() == {k[0]:v for k,v in gwps_counter.items()}
    eligible = set(gwps[gwps >= 200].index)-{'non-targeting'}
    cd4 = set(pd.read_parquet(CD4, columns=['target_gene']).target_gene)-{'non-targeting'}
    assert cd4 == {k[0] for k in streamed(CD4, ['target_gene'])}-{'non-targeting'}
    observed = set(by_target)-{'Non-Targeting'}
    names = sorted(set(declared)|eligible|cd4)
    rows = pd.DataFrame([dict(target=t, queried=t in queried, cached_cells=by_target.get(t,0),
        full_raw_gwps_ge200=t in eligible, cd4_observed=t in cd4) for t in names])
    assert set(rows.loc[rows.queried, 'target']) == set(declared)
    assert set(rows.loc[rows.cached_cells>0,'target']) == observed

    with np.load(PARENT/'artifacts/HCT116_offset_mapping.npz', allow_pickle=False) as d:
        ids = c.cell_integer_id.to_numpy(np.int64)
        block = ids//5000
        assert block.min() >= 0 and block.max() < len(d['logical_start'])
        absent = d['logical_start'][block] < 0
        missing_blocks = np.flatnonzero(d['logical_start']<0).tolist()
        version = int(d['version'])
    offsets = read(PARENT/'artifacts/HCT116_offset_verification.json')
    assert offsets['version'] == version and offsets['missing_blocks'] == missing_blocks
    listing = read(PARENT/'artifacts/hct_tree.json')
    cell_entries = [e for e in listing if e['type']=='file' and e['path'].startswith('data/HCT116/cells.lance/')]
    data_entries = [e for e in cell_entries if '/data/' in e['path'][len('data/'):]]
    assert len(data_entries) == 4 and all(e['size'] > 0 and len(e['oid'])==40 for e in cell_entries)
    readme = (PARENT/'artifacts/orion_readme.md').read_text()
    for token in ['18,903', 'eight million', 'CC-BY-NC-SA-4.0', 'Xaira-Therapeutics/X-Atlas-Orion']:
        assert token in readme
    index = read('artifacts/frozen_metadata_index.json')
    for p, meta in index['entries'].items():
        assert rt.sha(p) == meta['sha256'], p
    assert rt.sha(index['archive']) == index['archive_sha256']
    assert read('artifacts/calibration_gate.json')['status'] == 'CLOSED'
    check_bindings(plan)
    with TABLE.open('x') as f:
        rows.to_csv(f, sep='\t', index=False)
    saved = pd.read_csv(TABLE, sep='\t')
    pd.testing.assert_frame_equal(rows, saved)
    result = dict(status='PASS cached-metadata scope certificate only',
        utc=datetime.now(timezone.utc).isoformat(), plan_sha256=rt.sha(PLAN),
        cache=dict(cells=len(c), targets_including_control=len(by_target),
            perturbation_targets=len(observed), samples=len(by_sample),
            controls=by_target['Non-Targeting'], all_guide_pass=bool(c.pass_guide_filter.all()),
            unique_cell_ids=True, missing_query_labels=sorted(queried-set(by_target)),
            selected_columns=plan['columns']['hct']),
        query=dict(requested_targets=437, public_target_names=len(set(public)),
            h1_target_names=len(set(h1)), exact_target_union_replayed=True,
            predicate='gene_target IN (437 requested labels plus Non-Targeting)',
            guide_filter_in_original_query=False, sample_filter_in_original_query=False,
            limit_or_sampling_in_original_query=False,
            provenance='Present code and cache are consistent; no immutable original cell-query receipt'),
        excluded_name_counts=dict(full_raw_gwps_ge200=len(eligible),
            full_raw_gwps_ge200_not_queried=len(eligible-set(declared)),
            cd4_targets=len(cd4), cd4_targets_not_queried=len(cd4-set(declared)),
            cd4_and_full_raw_gwps_ge200=len(cd4&eligible),
            cd4_and_full_raw_gwps_ge200_not_queried=len((cd4&eligible)-set(declared))),
        expression_block_metadata=dict(expression_version=version, missing_blocks=missing_blocks,
            cache_rows_without_expression_block=int(absent.sum()),
            cache_perturbed_rows_without_expression_block=int((absent & (c.gene_target!='Non-Targeting')).sum()),
            cache_control_rows_without_expression_block=int((absent & (c.gene_target=='Non-Targeting')).sum()),
            cells_metadata_version='UNVERIFIED; not inferable from expression version'),
        cached_listing=dict(total_entries=len(listing), cell_files=len(cell_entries),
            cell_data_files=data_entries, listed_cell_file_bytes=sum(e['size'] for e in cell_entries),
            immutable_repository_commit='UNVERIFIED for this historical listing',
            completeness='UNVERIFIED:1000-entry listing has no retained pagination proof; no whole-release size claim'),
        table=dict(path=str(TABLE), sha256=rt.sha(TABLE), rows=len(rows), saved_recount='PASS'),
        conclusions=[
            'Cache is target-selected, not a full HCT116 census. Original query did not sample cells or filter samples/guide pass, but original query completeness and source identity remain unverified.',
            'Unqueried target labels overlap independent GWPS>=200 and CD4 metadata. Additional guide-eligible cells and sample strata may exist outside cache; presence, counts, expression support and feasibility are UNVERIFIED.',
            'No evidence here that extra cells exist for already queried targets, or that adding missing expression blocks fixes prior precision obstructions.',
            'A separately budgeted immutable full cells-metadata inventory can test omitted targets. Historical listed metadata data files sum to about135MB, within234MB reserve; current identities/schema/complete manifest and retry cost must be checked before allocation/acquisition.',
            'Keep prior cached-scope infeasibility certificates and failed rank gate unchanged. No lower-depth retry or response-driven target selection.'],
        frozen_metadata_members_verified=len(index['entries']),
        expression_payloads_read=False, biological_scores_read=False, network_io=False,
        task_activated=False, source_activated=False, holdout_scored=False, phase2_allowed=False)
    save_new(OUT,result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true')
    args=parser.parse_args()
    freeze() if args.freeze else main()

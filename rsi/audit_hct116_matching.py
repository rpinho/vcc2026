"""Frozen metadata-only HCT116 feasibility audit; never opens expression members."""
import argparse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt
import numpy as np
import pandas as pd

PLAN = Path('artifacts/HCT116_matching_plan_round40.json')
OUT = Path('artifacts/HCT116_matching_audit_round40.json')
CELLS = Path('../model_v4/artifacts/HCT116_cells.parquet')
BULK = Path('artifacts/sources/K562_GWPS_bulk.npz')
RAW = Path('../public_data/model_v4_signature_banks/K562_GWPS_raw.npz')


def load(p):
    return rt.json.loads(Path(p).read_text())


def freeze():
    assert not PLAN.exists() and not OUT.exists(), 'Preserve prior evidence'
    paths = [str(CELLS), str(BULK), str(RAW), __file__, 'recipes.py',
             'vendor/v4_model.py', 'artifacts/tasks/manifest.json',
             'artifacts/benchmark_roles.json', 'artifacts/calibration_gate.json',
             'artifacts/download_allocations.json', 'artifacts/frozen_metadata_index.json',
             'artifacts/call_control_native_plan_round39.json',
             'artifacts/call_control_native_audit_round39.json',
             '../competitors/COMPETITORS.md',
             'artifacts/tasks/HCT116_rich_search_selection.parquet']
    paths += [str(p) for p in Path('artifacts/scores').glob('*/*_s42.json')]
    rt.save(PLAN, dict(status='FROZEN before new HCT116 census calculations',
        utc=datetime.now(timezone.utc).isoformat(),
        preserve_sha256={p: rt.sha(p) for p in paths},
        columns=['cell_integer_id', 'gene_target', 'sample', 'pass_guide_filter'],
        source_members=['targets', 'n_cells'],
        population='Cached HCT116 guide-filtered cells only; not the full release. Exclude H1 sealed target identities. Require >=200 measured source cells in BOTH prepared legacy GWPS bulk and original v4 raw GWPS banks. No destination HCT116, HEK293T, RPE1 or H1 responses enter source eligibility.',
        precision='400 real cells/target is the desired precision, matching emitted depth and H1. Also report feasibility at the existing CD4192 depth and original64 depth without activating either. Panel sizes fixed at79 and300. No lower-depth retry based on scores.',
        control_supply='One common integer quota per sample for all targets; disjoint input/reference arms each8x response quota (H1 precedent). Thus quota <= floor(control_count/16). Also report unlimited-control bounds separately.',
        estimands='For all observed and eligible targets: per-sample kth-largest count upper bounds at k79/300, with/without control cap; counts of targets reaching400 total cells. For eligible only, choose first k by SHA256(hct116-match-r40:<target>), report exact minimum per-sample quota capped by controls//16. No adaptive target or sample search.',
        verification='Independent Counter of selected row tuples versus pandas crosstab; scalar sorted order statistics versus NumPy partition; exact fixed-subset minimum counts and control inequalities. Validate prior native audit PASS/rc0/current bindings and frozen metadata. Recheck bindings before publication.',
        limits='Feasibility only. A coordinate upper bound is necessary, not jointly achievable. A fixed hash subset below a threshold is not global impossibility. Sample labels may not capture guide/processing confounding. No task registration, source-bank extension, acquisition, expression/anchor/model score or calibration-gate replacement.',
        expression_payloads_read=False, new_task_registered=False,
        holdout_scored=False, phase2_allowed=False))
    print('FROZEN', PLAN, flush=True)


def bindings(plan):
    for p, h in plan['preserve_sha256'].items():
        assert rt.sha(p) == h, p


def main():
    assert not OUT.exists(), 'Preserve prior evidence'
    plan = load(PLAN)
    bindings(plan)
    prior_plan = load('artifacts/call_control_native_plan_round39.json')
    prior = load('artifacts/call_control_native_audit_round39.json')
    assert prior['status'] == 'PASS structural audit only'
    assert prior['plan_sha256'] == rt.sha('artifacts/call_control_native_plan_round39.json')
    bindings(prior_plan)
    for e in prior_plan['tasks'].values():
        assert rt.sha(e['path']) == e['expected_sha256']
    ends = [rt.json.loads(s) for s in Path('ledger/steps.jsonl').read_text().splitlines()]
    assert [e for e in ends if e.get('name') == 'audit_call_control_native_round39' and e['event'] == 'end'][-1]['rc'] == 0
    index = load('artifacts/frozen_metadata_index.json')
    for p, e in index['entries'].items():
        assert rt.sha(p) == e['sha256']
    assert rt.sha(index['archive']) == index['archive_sha256']
    obs = pd.read_parquet(CELLS, columns=plan['columns'])
    assert obs.cell_integer_id.is_unique and not obs.isna().any().any()
    obs = obs.loc[obs.pass_guide_filter].copy()
    obs['gene_target'] = obs.gene_target.astype(str)
    obs['sample'] = obs['sample'].astype(str)
    counter = Counter(zip(obs.gene_target, obs['sample']))
    tab = pd.crosstab(obs.gene_target, obs['sample']).sort_index().sort_index(axis=1)
    for t, row in tab.iterrows():
        assert row.tolist() == [counter[t, b] for b in tab.columns]
    assert sum(counter.values()) == int(tab.to_numpy().sum()) == len(obs)
    controls = tab.loc['Non-Targeting'].to_numpy(np.int64)
    tab = tab.drop(index='Non-Targeting')
    sealed = set(load('artifacts/tasks/manifest.json')['H1_holdout']['targets'])
    coverage = {}
    for label, path in [('legacy_bulk', BULK), ('v4_raw', RAW)]:
        with np.load(path, allow_pickle=False) as d:
            targets = d['targets'].astype(str)
            counts = d['n_cells']
        assert len(targets) == len(counts) and len(set(targets)) == len(targets)
        coverage[label] = {t for t, n in zip(targets, counts) if n >= 200}
    eligible = sorted(set(tab.index) & coverage['legacy_bulk'] & coverage['v4_raw'] - sealed)
    scopes = {'all_observed': list(tab.index), 'symmetric_source_ge200_unsealed': eligible}
    result = dict(status='RUNNING', plan_sha256=rt.sha(PLAN),
                  cached_filtered_cells=len(obs), samples=len(tab.columns),
                  controls_per_sample=dict(zip(tab.columns, map(int, controls))),
                  source_counts={k: len(v) for k, v in coverage.items()}, scopes={})
    for scope, ts in scopes.items():
        c = tab.loc[ts].to_numpy(np.int64)
        rec = dict(targets=len(ts), targets_with400_total=int((c.sum(1) >= 400).sum()), panels={})
        for k in [79, 300]:
            if len(ts) < k:
                rec['panels'][str(k)] = dict(status='INSUFFICIENT ELIGIBLE TARGETS')
                continue
            upper = np.partition(c, len(c)-k, axis=0)[len(c)-k]
            independent = [sorted([counter[t, b] for t in ts], reverse=True)[k-1] for b in tab.columns]
            assert upper.tolist() == independent
            capped = np.minimum(upper, controls//16)
            r = dict(unlimited_control_upper_bound=int(upper.sum()),
                     two_8x_control_arms_upper_bound=int(capped.sum()),
                     impossible_at_depth={str(d): bool(capped.sum() < d) for d in [64, 192, 400]},
                     per_sample_upper=dict(zip(tab.columns, map(int, capped))))
            if scope == 'symmetric_source_ge200_unsealed':
                order = sorted(ts, key=lambda t: (rt.hashlib.sha256(('hct116-match-r40:'+t).encode()).hexdigest(), t))
                chosen = order[:k]
                quota = np.minimum(tab.loc[chosen].to_numpy().min(0), controls//16)
                check = [min([counter[t, b] for t in chosen]+[counter['Non-Targeting', b]//16]) for b in tab.columns]
                assert quota.tolist() == check
                assert np.all(tab.loc[chosen].to_numpy() >= quota)
                assert np.all(16*quota <= controls)
                r['fixed_hash_subset'] = dict(targets=chosen, depth=int(quota.sum()),
                    per_sample_quota=dict(zip(tab.columns, map(int, quota))),
                    feasible_at_depth={str(d): bool(quota.sum() >= d) for d in [64, 192, 400]})
            rec['panels'][str(k)] = r
        result['scopes'][scope] = rec
    old = pd.read_parquet('artifacts/tasks/HCT116_rich_search_selection.parquet',
                          columns=['cell_integer_id', 'target_gene', 'evaluation_role', 'batch'])
    assert old.cell_integer_id.is_unique
    lookup = obs.set_index('cell_integer_id')
    assert set(old.cell_integer_id) <= set(lookup.index)
    selected = lookup.loc[old.cell_integer_id]
    assert selected['sample'].tolist() == old.batch.astype(str).tolist()
    result['cached_rich_selection'] = dict(cells=len(old),
        all_ids_in_cached_metadata=True, expression_file_present=Path('artifacts/tasks/HCT116_rich_search.h5ad').exists())
    bindings(plan)
    result.update(status='PASS metadata counts, exact bounds and deterministic witness quotas',
                  frozen_metadata_members_verified=len(index['entries']),
                  prior_native_audit_verified=True,
                  prior_native_summary={t: {v: {k: e[v][k] for k in ['n_cells', 'n_distinct_depths', 'equal_depth_pairs', 'movable_rows', 'movable_genes']} for v in ['all_input', 'first400']} for t, e in prior['tasks'].items()},
                  expression_payloads_read=False, new_task_registered=False,
                  holdout_scored=False, phase2_allowed=False,
                  utc=datetime.now(timezone.utc).isoformat())
    rt.save(OUT, result)
    print(result['status'], flush=True)
    for scope, rec in result['scopes'].items():
        print(scope, 'targets', rec['targets'], 'panels', {k: {x: y for x, y in v.items() if x in ['status', 'unlimited_control_upper_bound', 'two_8x_control_arms_upper_bound']} for k, v in rec['panels'].items()}, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--freeze', action='store_true')
    freeze() if p.parse_args().freeze else main()

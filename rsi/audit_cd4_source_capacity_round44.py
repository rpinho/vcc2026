"""Frozen count-only source-comparable CD4 supply certificate.

No expression members, target selection, prediction or metric access. Pandas
counts are checked by a separate streamed Arrow/Counter pass; scalar sorted
order statistics check every panel supply bound.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import runtime as rt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

PLAN = Path('artifacts/CD4_source_capacity_plan_round44.json')
OUT = Path('artifacts/CD4_source_capacity_audit_round44.json')
TABLE = Path('artifacts/CD4_source_capacity_targets_round44.tsv')
CD4 = 'artifacts/CD4_cell_census.parquet'
HCT = '../model_v4/artifacts/HCT116_cells.parquet'
SELECTION = 'artifacts/tasks/CD4_depth192_search_selection.parquet'
BANKS = dict(legacy='artifacts/cd4_sources/K562_GWPS_bulk.npz',
             raw='artifacts/cd4_banks/K562_GWPS_raw.npz',
             prepared_hct='artifacts/cd4_banks/HCT116.npz')


def load(path):
    return rt.json.loads(Path(path).read_text())


def freeze():
    assert not any(p.exists() for p in [PLAN, OUT, TABLE]), 'Preserve evidence'
    # Bind existing evidence without recalculating any new counts.
    paths = set(load('artifacts/source_support_plan_round43_retry.json')['preserve_sha256'])
    paths.update([__file__, CD4, HCT, SELECTION, *BANKS.values(),
                  'artifacts/source_support_plan_round43_retry.json',
                  'artifacts/source_support_audit_round43.json',
                  'artifacts/source_support_targets_round43.tsv',
                  'artifacts/live_calibration_round43.json',
                  'artifacts/benchmark_depth192_revision.json',
                  'artifacts/tasks/CD4_depth192_search_blueprint.json',
                  'artifacts/HCT116_matching_audit_round40.json',
                  '../competitors/COMPETITORS.md'])
    rt.save(PLAN, dict(status='FROZEN before new count calculations',
        utc=datetime.now(timezone.utc).isoformat(),
        preserve_sha256={p: rt.sha(p) for p in sorted(paths)},
        cd4_columns=['source_row', 'target_gene', 'batch', 'low_quality'],
        hct_columns=['cell_integer_id', 'gene_target', 'pass_guide_filter'],
        selection_columns=['source_row', 'target_gene', 'batch', 'evaluation_role'],
        npz_members=['targets', 'n_cells'],
        population='Cached low-quality-excluded CD4 census; all24 lanes. Exclude non-targeting and sealed H1 target identities. Require >=200 cells in BOTH existing CD4-wide raw and legacy GWPS banks, retaining the original192 design source threshold. This implies K562>=20.',
        scopes='Base symmetric GWPS200; intersect separately with prepared HCT116>=20, cached guide-filtered HCT116 census>=20, and their union. Cached counts are potential biological supply, never prepared source activation or compatible-expression certification.',
        precision='Fixed79-target panel and existing192-cell precision; desired400 also reported. No replacement panel, target ordering, lower depth, source activation, acquisition or scientific score.',
        bounds='Per scope: count targets with >=8 in every lane; uniform common-depth ceiling=24 times79th-largest target lane minimum. For arbitrary common all-positive lane quotas, each quota <=79th-largest count in that lane; sum is a necessary unlimited-control upper bound, not a feasible depth. Report per-lane bounds and impossibility at192/400. Fewer79 targets or any zero lane upper bound implies infeasibility.',
        controls='Independently validate inherited input/reference arms disjoint, census-bound and exactly80/32 per lane. These fixed uniform arms require uniform target mixtures. Exactly400 is not divisible by24 and is impossible with these exact arms; uniform408 supply is only a descriptive >=400 comparison, never an activated depth. Unlimited-control bounds remain valid obstructions for any expanded matched disjoint arms.',
        summaries='Full CD4 count-only target table, per-scope totals/bounds, fixed current79 coverage, and prepared/cache HCT source availability contingency among base targets and base targets meeting all-lane8. Census totals distinguish absent cached biology from missing bank rows; no panel selected.',
        verification='Pandas crosstab versus independent streamed Arrow/Counter; NPZ allowlist and scalar/NumPy eligibility; scalar sort versus NumPy order statistics; exhaustive tiny-subset bound fixture; current prior43 PASS/rc0 bindings and all frozen metadata/archive preserved; rehash before publication.',
        expression_payloads_decoded=False, network_reads=False,
        new_task_registered=False, holdout_scored=False, phase2_allowed=False))
    print('FROZEN', PLAN, flush=True)


def bindings(plan):
    for p, h in plan['preserve_sha256'].items():
        assert rt.sha(p) == h, p


def streamed_counts(path, target, batch=None, filtered=None):
    cols = [target] + ([batch] if batch else []) + ([filtered] if filtered else [])
    result = Counter()
    for piece in pq.ParquetFile(path).iter_batches(batch_size=65536, columns=cols):
        data = piece.to_pydict()
        for i, t in enumerate(data[target]):
            if filtered and not data[filtered][i]:
                continue
            key = (str(t), str(data[batch][i])) if batch else str(t)
            result[key] += 1
    return result


def bank_counts(path):
    with np.load(path, allow_pickle=False) as d:
        ts = d['targets'].astype(str)
        ns = d['n_cells']
    assert len(ts) == len(ns) == len(set(ts))
    assert np.isfinite(ns).all() and (ns >= 0).all() and (ns == np.floor(ns)).all()
    return dict(zip(ts, map(int, ns)))


def bound_fixture():
    c = np.array([[9, 1, 4], [1, 9, 4], [5, 5, 2], [2, 2, 9]])
    k = 2
    upper = np.sort(c, axis=0)[-k]
    uniform = 3 * sorted(c.min(1), reverse=True)[k-1]
    for subset in combinations(range(len(c)), k):
        quotas = c[list(subset)].min(0)
        assert np.all(quotas <= upper) and quotas.sum() <= upper.sum()
        assert 3 * quotas.min() <= uniform
    # Coordinate bounds may be mutually unattainable.
    assert max(int(c[list(s)].min(0).sum()) for s in combinations(range(4), 2)) < upper.sum()
    return dict(subsets=6, status='PASS exhaustive necessary bounds; not sufficiency')


def main():
    assert not OUT.exists() and not TABLE.exists(), 'Preserve evidence'
    plan = load(PLAN)
    bindings(plan)
    prior = load('artifacts/source_support_audit_round43.json')
    assert prior['status'].startswith('PASS')
    prior_plan = load('artifacts/source_support_plan_round43_retry.json')
    bindings(prior_plan)
    ends = [rt.json.loads(s) for s in Path('ledger/steps.jsonl').read_text().splitlines()]
    assert [e for e in ends if e.get('name') == 'audit_source_support_round43_retry' and e['event'] == 'end'][-1]['rc'] == 0
    index = load('artifacts/frozen_metadata_index.json')
    for p, e in index['entries'].items():
        assert rt.sha(p) == e['sha256'], p
    assert rt.sha(index['archive']) == index['archive_sha256']
    fixture = bound_fixture()

    obs = pd.read_parquet(CD4, columns=plan['cd4_columns'])
    assert obs.source_row.is_unique and not obs.isna().any().any()
    assert not obs.low_quality.any(), 'Frozen census should already exclude low quality'
    tab = pd.crosstab(obs.target_gene, obs.batch).sort_index().sort_index(axis=1)
    counter = streamed_counts(CD4, 'target_gene', 'batch')
    assert sum(counter.values()) == len(obs) == int(tab.to_numpy().sum())
    for t, row in tab.iterrows():
        assert row.tolist() == [counter[str(t), str(b)] for b in tab.columns]
    lanes = list(map(str, tab.columns))
    assert len(lanes) == 24
    controls = tab.loc['non-targeting'].to_numpy(np.int64)
    tab = tab.drop(index='non-targeting')

    sel = pd.read_parquet(SELECTION, columns=plan['selection_columns'])
    assert sel.source_row.is_unique
    lookup = obs.set_index('source_row')
    assert set(sel.source_row) <= set(lookup.index)
    for col in ['target_gene', 'batch']:
        assert lookup.loc[sel.source_row, col].astype(str).tolist() == sel[col].astype(str).tolist()
    arms = {}
    for role, quota in [('input_control', 80), ('reference_control', 32)]:
        arm = sel.loc[sel.evaluation_role == role]
        assert set(arm.target_gene.astype(str)) == {'non-targeting'}
        counts = Counter(arm.batch.astype(str))
        assert [counts[b] for b in lanes] == [quota]*24
        arms[role] = set(arm.source_row)
    assert arms['input_control'].isdisjoint(arms['reference_control'])
    assert np.all(controls >= 112)

    hct = pd.read_parquet(HCT, columns=plan['hct_columns'])
    assert hct.cell_integer_id.is_unique and not hct.isna().any().any()
    hct = hct.loc[hct.pass_guide_filter]
    hct_n = {str(t): int(n) for t, n in hct.gene_target.value_counts().items()}
    assert hct_n == dict(streamed_counts(HCT, 'gene_target', filtered='pass_guide_filter'))
    banks = {k: bank_counts(p) for k, p in BANKS.items()}
    sealed = set(load('artifacts/tasks/manifest.json')['H1_holdout']['targets'])
    rows = pd.DataFrame(index=tab.index.astype(str))
    rows.index.name = 'target'
    for k, counts in banks.items():
        rows[k+'_cells'] = [counts.get(t, 0) for t in rows.index]
    rows['cached_hct_cells'] = [hct_n.get(t, 0) for t in rows.index]
    rows['sealed_H1_target'] = rows.index.isin(sealed)
    rows['base_eligible'] = (rows.legacy_cells >= 200) & (rows.raw_cells >= 200) & ~rows.sealed_H1_target
    rows['prepared_hct_ge20'] = rows.prepared_hct_cells >= 20
    rows['cached_hct_ge20'] = rows.cached_hct_cells >= 20
    rows['minimum_lane_cells'] = tab.min(1).to_numpy()
    rows['total_cd4_cells'] = tab.sum(1).to_numpy()
    for b in lanes:
        rows[b] = tab[b].to_numpy()
    for t, row in rows.iterrows():
        scalar = t not in sealed and banks['legacy'].get(t, 0) >= 200 and banks['raw'].get(t, 0) >= 200
        assert scalar == bool(row.base_eligible)
        assert int(row.minimum_lane_cells) == min(counter[t, b] for b in lanes)

    base = rows.base_eligible
    scopes = dict(base_GWPS200=base,
        prepared_two_family=base & rows.prepared_hct_ge20,
        cached_two_family_potential=base & rows.cached_hct_ge20,
        union_two_family_potential=base & (rows.prepared_hct_ge20 | rows.cached_hct_ge20))
    results = {}
    for name, mask in scopes.items():
        ts = rows.index[mask].tolist()
        c = tab.loc[ts].to_numpy(np.int64)
        rec = dict(targets=len(ts), all_lane8_targets=int((c.min(1) >= 8).sum()) if len(ts) else 0,
                   all_lane17_targets=int((c.min(1) >= 17).sum()) if len(ts) else 0)
        if len(ts) < 79:
            rec.update(status='INSUFFICIENT79 TARGETS', impossible192=True, impossible400=True)
        else:
            upper = np.partition(c, len(c)-79, axis=0)[len(c)-79]
            scalar_upper = [sorted((counter[t, b] for t in ts), reverse=True)[78] for b in lanes]
            assert upper.tolist() == scalar_upper
            q = int(np.sort(c.min(1))[-79])
            assert q == sorted((min(counter[t, b] for b in lanes) for t in ts), reverse=True)[78]
            rec.update(status='NECESSARY SUPPLY BOUNDS ONLY',
                uniform_common_quota_ceiling=q, uniform_depth_ceiling=24*q,
                uniform192_supply_possible=q >= 8,
                uniform_at_least400_supply_possible=q >= 17,
                arbitrary_common_depth_upper_bound=int(upper.sum()),
                all24_positive_quotas_possible=bool((upper > 0).all()),
                per_lane_upper=dict(zip(lanes, map(int, upper))),
                impossible192=bool(upper.sum() < 192 or (upper == 0).any()),
                impossible400=bool(upper.sum() < 400 or (upper == 0).any()))
        results[name] = rec

    current = load('artifacts/tasks/manifest.json')['CD4_depth192_search']['targets']
    current_rows = rows.loc[current]
    assert len(current_rows) == 79 and (current_rows.minimum_lane_cells >= 8).all()
    current_support = dict(targets=79,
        base_eligible=int(current_rows.base_eligible.sum()),
        prepared_hct_ge20=int(current_rows.prepared_hct_ge20.sum()),
        cached_hct_ge20=int(current_rows.cached_hct_ge20.sum()),
        union_hct_ge20=int((current_rows.prepared_hct_ge20 | current_rows.cached_hct_ge20).sum()))
    assert current_support['prepared_hct_ge20'] == 6, 'Recheck prior43 identity'
    contingency = {}
    for label, mask in [('base', base), ('base_all_lane8', base & (rows.minimum_lane_cells >= 8))]:
        r = rows.loc[mask]
        contingency[label] = dict(total=len(r),
            prepared_and_cached=int((r.prepared_hct_ge20 & r.cached_hct_ge20).sum()),
            prepared_only=int((r.prepared_hct_ge20 & ~r.cached_hct_ge20).sum()),
            cached_only=int((~r.prepared_hct_ge20 & r.cached_hct_ge20).sum()),
            neither=int((~r.prepared_hct_ge20 & ~r.cached_hct_ge20).sum()))
        assert sum(v for k, v in contingency[label].items() if k != 'total') == len(r)
    bindings(plan)
    rows.to_csv(TABLE, sep='\t')
    result = dict(status='PASS metadata supply certificate only', utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=rt.sha(PLAN), table_sha256=rt.sha(TABLE),
        cd4_census_cells=len(obs), cd4_targets=len(rows), lanes=24,
        hct_guide_filtered_cells=len(hct), hct_observed_noncontrol_targets=len(set(hct_n)-{'Non-Targeting'}),
        current79_support=current_support, scopes=results, hct_support_contingency=contingency,
        inherited_controls=dict(input=1920, reference=768, disjoint=True, census_verified=True,
            input_per_lane=80, reference_per_lane=32,
            exact400_with_inherited_uniform_arms=False),
        independent_streamed_counter_checks=True, fixture=fixture,
        frozen_metadata_members_verified=len(index['entries']),
        preserved_bindings_verified=len(plan['preserve_sha256']),
        limits=plan['bounds']+' '+plan['controls']+' Cached HCT support is limited to the cached selection, not the entire release.',
        expression_payloads_decoded=False, network_reads=False,
        new_task_registered=False, holdout_scored=False, phase2_allowed=False)
    rt.save(OUT, result)
    print(result['status'], flush=True)
    print('current79', current_support, flush=True)
    print('contingency', contingency, flush=True)
    for k, v in results.items():
        print(k, {x: y for x, y in v.items() if x != 'per_lane_upper'}, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    freeze() if parser.parse_args().freeze else main()

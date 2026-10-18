"""Input-control-only scope audit of equal-depth integer exchanges.

No DE, reference controls, response expression, prediction or model selection.
Freeze the selection/estimands before opening any expression arrays.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
import json
from pathlib import Path

import runtime as rt
import h5py
import numpy as np
from anndata.io import read_elem

PLAN = Path('artifacts/call_control_native_plan_round39.json')
OUT = Path('artifacts/call_control_native_audit_round39.json')


def load(path):
    return json.loads(Path(path).read_text())


def freeze():
    assert not PLAN.exists() and not OUT.exists(), 'Preserve existing evidence'
    roles = load('artifacts/benchmark_roles.json')['primary_search']
    assert roles == ['H1_rich_search', 'CD4_depth192_search']
    paths = ['audit_call_control_native.py', 'artifacts/tasks/manifest.json',
             'artifacts/benchmark_roles.json', 'artifacts/calibration_gate.json',
             'artifacts/download_allocations.json', '../competitors/COMPETITORS.md',
             'artifacts/call_control_heterogeneous_plan_round38.json',
             'artifacts/call_control_heterogeneous_audit_round38.json',
             'artifacts/phase1_round36_audit.json']
    tasks = {}
    for task in roles:
        receipt = Path('artifacts/tasks') / (task + '.json')
        paths.append(str(receipt))
        tasks[task] = dict(path=str(receipt.with_suffix('.h5ad')),
                          expected_sha256=load(receipt)['sha256'])
    rt.save(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        status='FROZEN before reading native control expression', tasks=tasks,
        selection='All input_control rows, and first400 of those rows in existing task order. No resampling, depth selection, rounding or alternative subset.',
        estimands='Exact integer readout-axis row-depth histogram, equal-depth cell pairs, legal one-unit 2x2 donor choices and participating rows/genes; one lexicographically first witness per view if available.',
        donor_definition='For each unordered equal-depth row pair r<s, count ordered positive donor genes a in r and b in s with a!=b: nnz(r)*nnz(s)-support_overlap. Move decreases (r,a),(s,b) by1 and increases (r,b),(s,a) by1. Counts are moves available at the initial table, not independent directions or achievable call counts.',
        verification='Independent Python integer row sums, Counter depth histogram, direct support-set donor count, and exact Fraction witness CPM invariants versus NumPy production summaries. Recheck all bindings and task hashes after audit.',
        limits='Depths are sums over each task measured readout axis. Native controls are not a pinned predicted pseudobulk. First400 is a fixed diagnostic, not a production emitter or balanced template. Equal-depth exchange is sufficient, not necessary, for dual-moment preservation. No score, call-set feasibility conclusion, target ranking, benchmark repair or general impossibility claim.',
        preserve_sha256={p: rt.sha(p) for p in paths},
        biological_scored=False, holdout_scored=False, phase2_allowed=False))
    print('FROZEN', PLAN, flush=True)


def rows_from_controls(path):
    """Read only input-control CSR data slices, never other expression rows."""
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs'])
        chosen = np.flatnonzero(obs.evaluation_role.astype(str).to_numpy() == 'input_control')
        assert len(chosen) >= 400
        assert set(obs.iloc[chosen].target_gene.astype(str)) == {'non-targeting'}
        assert len(set(obs.iloc[chosen].source_row)) == len(chosen)
        assert set(obs.iloc[chosen].source_row).isdisjoint(set(obs.loc[obs.evaluation_role.astype(str) != 'input_control', 'source_row']))
        x = f['X']
        assert x.attrs['encoding-type'] == 'csr_matrix'
        shape = tuple(map(int, x.attrs['shape']))
        rows = []
        for i in chosen:
            lo, hi = map(int, x['indptr'][int(i):int(i)+2])
            values = x['data'][lo:hi]
            cols = x['indices'][lo:hi]
            assert np.all(np.isfinite(values)) and np.all(values >= 0)
            assert np.array_equal(values, np.rint(values))
            assert len(set(map(int, cols))) == len(cols)
            row = {int(j): int(v) for j, v in zip(cols, values) if v > 0}
            assert sum(row.values()) > 0
            rows.append(row)
        return rows, shape[1], chosen.tolist(), obs.iloc[chosen].source_row.astype(str).tolist()


def summarize(rows, genes):
    depths = np.array([np.fromiter(r.values(), dtype=np.int64).sum() for r in rows])
    unique, multiplicity = np.unique(depths, return_counts=True)
    histogram = {str(int(d)): int(m) for d, m in zip(unique, multiplicity)}
    independent = Counter(sum(map(int, row.values())) for row in rows)
    assert histogram == {str(d): n for d, n in sorted(independent.items())}
    groups = defaultdict(list)
    for i, row in enumerate(rows):
        groups[sum(row.values())].append(i)
    donor_choices = 0
    pairs = 0
    movable_rows, movable_genes = set(), set()
    witness = None
    # Scope is exactly the current support; no claimed global feasibility test.
    for depth, ids in sorted(groups.items()):
        for u, r in enumerate(ids):
            for s in ids[u+1:]:
                pairs += 1
                ar, bs = set(rows[r]), set(rows[s])
                # Independent count: sum available opposite donors for each a.
                expected = sum(len(bs) - int(a in bs) for a in ar)
                count = len(ar)*len(bs)-len(ar.intersection(bs))
                assert expected == count
                donor_choices += count
                if not count:
                    continue
                movable_rows.update((r, s))
                movable_genes.update(a for a in ar if len(bs-{a}))
                movable_genes.update(b for b in bs if len(ar-{b}))
                if witness is None:
                    a = min(a for a in ar if bs-{a})
                    b = min(bs-{a})
                    before = [dict(rows[r]), dict(rows[s])]
                    after = [dict(row) for row in before]
                    for row, gene, delta in [(0,a,-1),(0,b,1),(1,a,1),(1,b,-1)]:
                        after[row][gene] = after[row].get(gene, 0)+delta
                    assert all(v >= 0 for row in after for v in row.values())
                    assert [sum(row.values()) for row in before] == [sum(row.values()) for row in after]
                    for g in ar | bs:
                        assert sum(row.get(g, 0) for row in before) == sum(row.get(g, 0) for row in after)
                        assert sum((Fraction(row.get(g, 0), sum(row.values())) for row in before), Fraction()) == sum((Fraction(row.get(g, 0), sum(row.values())) for row in after), Fraction())
                    witness = dict(rows=[r,s], genes=[a,b], depth=depth, amount=1,
                                   raw_margins_and_exact_mean_cell_cpm='PASS')
    assert pairs == sum(n*(n-1)//2 for n in independent.values())
    supported = sum(n for n in independent.values() if n >= 2)
    return dict(n_cells=len(rows), n_genes=genes, min_depth=int(depths.min()),
                max_depth=int(depths.max()), n_distinct_depths=len(unique),
                depth_histogram=histogram, cells_in_repeated_depths=supported,
                equal_depth_pairs=pairs, legal_unit_donor_choices=donor_choices,
                movable_rows=len(movable_rows), movable_row_fraction=len(movable_rows)/len(rows),
                movable_genes=len(movable_genes), movable_gene_fraction=len(movable_genes)/genes,
                first_witness=witness,
                status='SUPPORTED limited initial-table moves' if donor_choices else 'UNSUPPORTED equal-depth operator at this table')


def verify(plan):
    for path, digest in plan['preserve_sha256'].items():
        assert rt.sha(path) == digest, path
    for entry in plan['tasks'].values():
        assert rt.sha(entry['path']) == entry['expected_sha256'], entry['path']


def main():
    assert not OUT.exists(), 'Preserve existing evidence'
    plan = load(PLAN)
    verify(plan)
    result = dict(status='RUNNING', plan_sha256=rt.sha(PLAN), tasks={},
                  biological_scored=False, holdout_scored=False, phase2_allowed=False)
    rt.save(OUT, result)
    for task, entry in plan['tasks'].items():
        rows, genes, indices, source_rows = rows_from_controls(entry['path'])
        result['tasks'][task] = dict(input_indices=indices, source_rows=source_rows,
                                    all_input=summarize(rows, genes),
                                    first400=summarize(rows[:400], genes))
        rt.save(OUT, result)
        print(task, {view: {k: v for k,v in result['tasks'][task][view].items() if k not in ['depth_histogram','first_witness']} for view in ['all_input','first400']}, flush=True)
    verify(plan)
    result.update(status='PASS structural audit only', utc=datetime.now(timezone.utc).isoformat())
    rt.save(OUT, result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    freeze() if args.freeze else main()

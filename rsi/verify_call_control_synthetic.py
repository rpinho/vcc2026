"""Independent pairwise/Fraction auditor and official synthetic comparison.

No biological task, source response, HOLDOUT, baseline or anchor matrix is read.
Freeze the plan before execution. Never overwrite a frozen plan or receipt.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import runtime as rt
import numpy as np
import call_control_synthetic as impl

PLAN = Path('artifacts/call_control_synthetic_plan_round37_retry.json')
OUT = Path('artifacts/call_control_synthetic_audit_round37_retry.json')
ROUNDOFF_SIGNS = []


def load(path):
    return json.loads(Path(path).read_text())


def freeze():
    import official_scoring as scoring
    import pdex, numba_mwu
    assert not PLAN.exists() and not OUT.exists(), 'Preserve previous evidence'
    paths = ['call_control_synthetic.py', 'verify_call_control_synthetic.py',
             'runtime.py', 'official_scoring.py', 'artifacts/calibration_gate.json',
             'artifacts/phase1_round36_audit.json', 'artifacts/tasks/manifest.json',
             'artifacts/benchmark_roles.json', 'artifacts/download_allocations.json',
             'artifacts/frozen_metadata_index.json', 'artifacts/frozen_metadata.tar.gz',
             '../competitors/COMPETITORS.md',
             'artifacts/protocol_revisions/call_control_synthetic_r37_initial/failure.json',
             'artifacts/call_control_synthetic_plan_round37.json',
             'artifacts/call_control_synthetic_audit_round37.json']
    deps = {}
    for module in [pdex, numba_mwu]:
        root = Path(module.__file__).parent
        deps[module.__name__] = {str(p.relative_to(root)): rt.sha(p)
                               for p in root.rglob('*.py')}
    rt.save(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        status='FROZEN after preserved initial failure, before revised audit execution',
        revision='Auditor-only: official all-tied uncalled LFC3.203426503814917e-16 versus exact0; allow only uncalled sign disagreement within unchanged pre-existing1e-12 LFC absolute tolerance, retain strict called signs. No input/operator/threshold/arithmetic/margin change.',
        purpose='Competitor priority2 mechanics only; no biological recipe/selection',
        input_domain='Deterministic synthetic integer counts only; main depth2560, auxiliary power-of-two depths; exact representable CPM values',
        synthetic_main='400pred x640 genes /64controls;150 paired genes initially(2,6)/(8,0), desired(5,3);170 neutral pairs(4,4); strict3-source synthetic consensus; self G639',
        main_expected={'before_calls': 0, 'after_calls': 300, 'selected_genes': list(range(300)),
                       'preserve': ['every row depth', 'every gene total', 'integer nonnegative counts']},
        other_cases=['unequal n and depth, sparse zeros, control-only CPM gate', 'all tied p1',
                     'all two-gene tables with4rows/depth2/totals6,2; requests +1,-1 /+1,0 /0,0 /-1,+1'],
        arithmetic='Independent pairwise twice-U using exact count/depth Fractions; pooled tie frequencies; exact rational variance; continuity0.5; two-sided erfc; independent all-upper-ranks BH; arithmetic CPM LFC epsilon1e-9',
        tolerance={'relative': 1e-10, 'absolute_p_q': 1e-300, 'absolute_lfc': 1e-12},
        acceptance='All evaluator/direct/official retained features, p,q,LFC and calls match within frozen tolerance; exact U/ties/margins; replay each2x2 move; exhaustive independent tiny enumeration; reject invalid and unsupported inputs. Report scope and infeasible requests.',
        limits='Synthetic feasibility only. Uniform controls/uniform output are artificial, not a realism or dispersion check, no biological gain or general feasibility solver. No general hidden GPU numerical equivalence claim.',
        official=scoring.provenance(), dependency_sha256=deps,
        preserve_sha256={p: rt.sha(p) for p in paths},
        biological_scored=False, holdout_scored=False, phase2_allowed=False))
    print('FROZEN', PLAN, flush=True)


def direct(pred, controls):
    """Pairwise comparisons on exact rational CPM; no production rank helpers."""
    n, m = len(pred), len(controls)
    out = []
    for g in range(pred.shape[1]):
        x = [F(int(row[g])*1000000, int(sum(row))) for row in pred]
        y = [F(int(row[g])*1000000, int(sum(row))) for row in controls]
        my = sum(y, F(0))/m
        if my <= 5:
            continue
        # Vectorized comparison is still direct all-pairs, independently of search ranks.
        # Fractions are exactly representable in these frozen power-of-two fixtures.
        a, b = np.array(x, dtype=float), np.array(y, dtype=float)
        twice_u = int((2*(a[:, None] > b[None, :])+(a[:, None] == b[None, :])).sum())
        frequencies = Counter(x+y)
        ties = sum(t**3-t for t in frequencies.values())
        variance = F(n*m, 12)*(n+m+1-F(ties, (n+m)*(n+m-1)))
        deviation = abs(F(twice_u, 2)-F(n*m, 2))-F(1, 2)
        p = 1.0 if variance <= 0 else min(1.0, math.erfc(float(deviation)/math.sqrt(2*float(variance))))
        mx = sum(x, F(0))/n
        lfc = math.log2((float(mx)+1e-9)/(float(my)+1e-9))
        out.append(dict(gene=g, u=twice_u/2, ties=ties, variance=float(variance), p=p, lfc=lfc))
    # BH by counting the full upper set, rather than production's ordered suffix scan.
    upper_candidates = [(s['p'], s['p']*len(out)/sum(t['p'] <= s['p'] for t in out))
                        for s in out]
    for r in out:
        r['q'] = min([1.0]+[q for p, q in upper_candidates if p >= r['p']])
        r['call'] = r['q'] < 0.05
        r['sign'] = int(np.sign(r['lfc']))
    return out


def close(x, y, lfc=False):
    assert math.isclose(x, y, rel_tol=1e-10, abs_tol=1e-12 if lfc else 1e-300), (x, y)


def official(pred, controls, sparse_input=True):
    import anndata as ad
    import pandas as pd
    from scipy import sparse
    import official_scoring as scoring
    from dataclasses import replace
    x = np.vstack([controls, pred])
    genes = ['G'+str(i) for i in range(x.shape[1])]
    labels = ['non-targeting']*len(controls)+[genes[-1]]*len(pred)
    data = ad.AnnData(sparse.csr_matrix(x) if sparse_input else x,
        obs=pd.DataFrame({'target_gene': labels}, index=[f'cell{i}' for i in range(len(x))]),
        var=pd.DataFrame(index=genes))
    cfg = replace(scoring.config(), num_threads=1)
    table = scoring.ORIGINAL_DE(data, cfg=cfg, fp=None, store=None, side='pred')
    return {int(row['feature'][1:]): row for row in table.to_dicts()}


def certify(pred, controls, dense=False):
    fast, independent = impl.summarize(pred, controls), direct(pred, controls)
    expected = official(pred, controls)
    if dense:
        dense_rows = official(pred, controls, sparse_input=False)
        assert expected.keys() == dense_rows.keys()
        for g in expected:
            for key in ['p_value', 'p_adj', 'log2_fold_change']:
                close(expected[g][key], dense_rows[g][key], key == 'log2_fold_change')
    assert [r['gene'] for r in fast] == [r['gene'] for r in independent] == sorted(expected)
    for a, b in zip(fast, independent):
        assert a['u'] == b['u'] and a['ties'] == b['ties']
        for key in ['variance', 'p', 'q', 'lfc']:
            close(a[key], b[key], key == 'lfc')
        row = expected[a['gene']]
        close(a['p'], row['p_value'])
        close(a['q'], row['p_adj'])
        close(a['lfc'], row['log2_fold_change'], True)
        assert a['call'] == b['call'] == (row['p_adj'] < 0.05)
        if not (a['sign'] == b['sign'] == int(np.sign(row['log2_fold_change']))):
            assert not a['call'] and not b['call'] and row['p_adj'] >= 0.05
            assert max(abs(a['lfc']), abs(b['lfc']), abs(row['log2_fold_change'])) <= 1e-12
            ROUNDOFF_SIGNS.append(dict(gene=a['gene'], fast_lfc=a['lfc'],
                exact_lfc=b['lfc'], official_lfc=row['log2_fold_change'],
                q=row['p_adj'], scope='uncalled, within frozen original LFC tolerance'))
    return fast


def compositions(total, caps, prefix=()):
    """Independent recursive bounded composition census, no itertools product."""
    if not caps:
        if total == 0:
            yield prefix
        return
    for value in range(max(0, total-sum(caps[1:])), min(caps[0], total)+1):
        yield from compositions(total-value, caps[1:], prefix+(value,))


def digest(x):
    return hashlib.sha256(np.ascontiguousarray(x, dtype='<i8').tobytes()).hexdigest()


def main():
    import official_scoring as scoring
    import pdex, numba_mwu
    plan = load(PLAN)
    assert not OUT.exists(), 'Preserve completed/failed output before revising'
    for path, sha in plan['preserve_sha256'].items():
        assert rt.sha(path) == sha, path
    assert scoring.provenance() == plan['official']
    for module in [pdex, numba_mwu]:
        root = Path(module.__file__).parent
        assert {str(p.relative_to(root)): rt.sha(p) for p in root.rglob('*.py')} == plan['dependency_sha256'][module.__name__]
    receipt = dict(status='RUNNING', utc=datetime.now(timezone.utc).isoformat(),
                   plan_sha256=rt.sha(PLAN), cases={}, biological_scored=False,
                   holdout_scored=False, phase2_allowed=False)
    rt.save(OUT, receipt)
    fixtures = impl.fixtures()
    pred, controls = fixtures.pop('main')
    before = certify(pred, controls)
    assert sum(r['call'] for r in before) == 0
    votes = np.zeros((3, 640), dtype=np.int64)
    votes[:, :300:2] = 1; votes[:, 1:300:2] = -1
    selected = impl.consensus(votes, controls.mean(0) > 0, 300, 639)
    assert selected['status'] == 'SELECTED' and selected['genes'] == list(range(300))
    after = pred.copy()
    moves = []
    for a in range(0, 300, 2):
        after, batch = impl.balance_pair(after, a, a+1)
        moves.extend(batch)
    replay = pred.copy()
    for r, s, a, b, d in moves:
        assert d > 0
        replay[r, a] += d; replay[r, b] -= d
        replay[s, a] -= d; replay[s, b] += d
        assert np.all(replay[[r, s]][:, [a, b]] >= 0)
    assert np.array_equal(replay, after)
    assert np.array_equal(pred.sum(0), after.sum(0))
    assert np.array_equal(pred.sum(1), after.sum(1))
    assert np.all(after >= 0) and np.issubdtype(after.dtype, np.integer)
    final = certify(after, controls)
    called = [r for r in final if r['call'] and r['gene'] != 639]
    assert [r['gene'] for r in called] == selected['genes']
    assert [r['sign'] for r in called] == selected['signs']
    receipt['cases']['main'] = dict(status='PASS', before_calls=0, after_calls=len(called),
        uncalled=640-len(called), moves=len(moves), row_depth=int(after.sum(1)[0]),
        pred_before_sha256=digest(pred), pred_after_sha256=digest(after), controls_sha256=digest(controls),
        exact_row_and_column_margins=True, self_excluded=639,
        worst_called_q=max(r['q'] for r in called), uncalled_min_q=min(r['q'] for r in final if not r['call']),
        limitation='Artificial equal-depth/equal-control fixture, not a prediction/dispersion validation')
    rt.save(OUT, receipt)
    for name, (p, c) in fixtures.items():
        rows = certify(p, c, dense=True)
        if name == 'all_tied':
            assert all(r['p'] == r['q'] == 1 for r in rows)
        if name == 'unequal_depth_filter':
            assert 0 not in [r['gene'] for r in rows]
        receipt['cases'][name] = dict(status='PASS', retained_genes=[r['gene'] for r in rows],
                                     called_genes=[r['gene'] for r in rows if r['call']])
    tiny_c = fixtures['tiny'][1]
    direct_states = {}
    for column in compositions(6, [2]*4):
        x = np.array([[v, 2-v] for v in column], dtype=np.int64)
        rows = certify(x, tiny_c)
        direct_states[tuple(column)] = [r['sign'] if r['call'] else 0 for r in rows]
    queries = {}
    for request in [(1, -1), (1, 0), (0, 0), (-1, 1)]:
        result = impl.enumerate_two_gene(tiny_c, [2]*4, [6, 2], request)
        assert {tuple(r['column']): r['state'] for r in result['tables']} == direct_states
        witnesses = [list(col) for col, state in direct_states.items() if state == list(request)]
        assert (result['status'] == 'FEASIBLE') == bool(witnesses)
        assert result['checked'] == len(direct_states)
        if witnesses:
            assert result['witness'] in witnesses
        queries[str(request)] = {k: v for k, v in result.items() if k != 'tables'}
    receipt['cases']['exhaustive'] = dict(status='PASS', table_count=len(direct_states), requests=queries,
        scope='Every integer nonnegative4x2 table with depths2 and column totals6,2; exact finite certificate, not general solver')
    # Equal column means and mean ranks do not alone fix tie-corrected p-values.
    tied_a = np.array([[1, 1], [1, 1], [2, 0], [2, 0]], dtype=np.int64)
    tied_b = np.array([[0, 2], [2, 0], [2, 0], [2, 0]], dtype=np.int64)
    stats_a, stats_b = direct(tied_a, tiny_c), direct(tied_b, tiny_c)
    assert np.array_equal(tied_a.sum(0), tied_b.sum(0))
    assert [r['u'] for r in stats_a] == [r['u'] for r in stats_b]
    assert all(a['p'] != b['p'] and a['ties'] != b['ties'] for a, b in zip(stats_a, stats_b))
    receipt['cases']['tie_sufficiency'] = dict(status='PASS',
        equal_means_and_u=True, a=stats_a, b=stats_b,
        implication='Pooled tie frequencies are also needed for pinned CPU Wilcoxon p-values; not just predicted means and mean ranks')
    rejects = []
    for name, action in [
        ('fractional', lambda: impl.normalized(np.ones((2, 2), dtype=float))),
        ('negative', lambda: impl.normalized(np.array([[-1, 3]]))),
        ('zero_depth', lambda: impl.normalized(np.zeros((1, 2), dtype=int))),
        ('unsupported_pair', lambda: impl.balance_pair(np.array([[1, 1], [1, 2]]), 0, 1)),
        ('enumeration_limit', lambda: impl.enumerate_two_gene(tiny_c, [8]*8, [48, 16], [1, -1]))]:
        try:
            action()
        except ValueError:
            rejects.append(name)
        else:
            raise AssertionError(name)
    assert impl.consensus(votes, np.ones(640, dtype=bool), 301, 639)['status'].startswith('INFEASIBLE')
    assert impl.enumerate_two_gene(tiny_c, [2]*4, [6, 3], [1, -1])['status'] == 'INFEASIBLE margins'
    receipt['cases']['input_rejection'] = dict(status='PASS', rejected=rejects,
        selection_shortfall_reported=True, inconsistent_margins_reported=True)
    # Frozen scientific inputs and scorer environment remain unchanged after all tests.
    for path, sha in plan['preserve_sha256'].items():
        assert rt.sha(path) == sha, path
    assert scoring.provenance() == plan['official']
    receipt.update(status='PASS synthetic mechanics only', completed_utc=datetime.now(timezone.utc).isoformat(),
                   limits=plan['limits'], uncalled_roundoff_signs=ROUNDOFF_SIGNS,
                   gate_unchanged_sha256=rt.sha('artifacts/calibration_gate.json'))
    rt.save(OUT, receipt)
    print('PASS synthetic call-control mechanics; no biological score', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    freeze() if args.freeze else main()

"""Frozen synthetic stress plan and independent official-DE certification."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import json
from pathlib import Path

import runtime as rt
import numpy as np
import call_control_heterogeneous as impl
import call_control_synthetic as rank_impl
import verify_call_control_synthetic as audit

PLAN = Path('artifacts/call_control_heterogeneous_plan_round38.json')
OUT = Path('artifacts/call_control_heterogeneous_audit_round38.json')
REQUESTS = [(1, -1), (1, 0), (0, 0), (-1, 1)]


def load(p):
    return json.loads(Path(p).read_text())


def freeze():
    import official_scoring as scoring
    import pdex, numba_mwu
    assert not PLAN.exists() and not OUT.exists(), 'Preserve previous evidence'
    previous = load(audit.PLAN)
    paths = list(previous['preserve_sha256']) + [str(audit.PLAN), str(audit.OUT),
        'call_control_heterogeneous.py', 'verify_call_control_heterogeneous.py',
        'artifacts/live_calibration_round37.json']
    rt.save(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        status='FROZEN before any heterogeneous result or official call',
        purpose='Competitor priority2 synthetic mechanics, with variable depth and nonzero dispersion',
        justification='Raw column totals do not fix mean-cell CPM under variable depths. Equal-depth exchanges are a sufficient exact dual-moment domain; test their scope and enumerate a tiny wider domain independently. This uses no biological label or benchmark repair.',
        construction='400pred/64controls/640genes, factors[1,2,4,2], per-depth balanced parity. Controls and neutral pairs(3,5)/(5,3);150 selected pairs initially(2,6)/(8,0), goal(4,4)/(6,2); all scaled by depth factor. Three fixed unanimous source votes; selfG639 excluded. Request300 calls; no seeds.',
        expected_main={'before_calls': 0, 'after_calls': 300, 'unchanged_other_genes': 340,
                       'raw_depths': [2560, 5120, 10240],
                       'preserve': ['integer nonnegative counts', 'every row depth',
                           'raw gene totals', 'per-depth gene totals', 'exact mean-cell CPM'],
                       'dispersion': 'Every retained predicted and control gene has positive CPM variance; no claim of biological realism or matched variance'},
        tiny='4x2 depths[2,4,8,4], totals[11,7], initial a[1,2,5,3], fraction sum19/8;32controls depth8/16/32/16 and fractions3/8,1/2,1/2,5/8. Exhaust all raw-margin tables and exact fraction19/8 subset; independently certify EVERY raw table with official DE. Fixed four requests in each domain.',
        counterexamples={'raw_only': 'Exchange a -1 at row0,+1 at row2; raw margins preserved, mean-cell CPM delta -93750 for gene0,+93750 gene1.',
                         'sufficient_not_necessary': 'a[0,4,3,4] has same raw/CPM margins as initial; operator must report UNSUPPORTED due to different per-depth totals; exact census retains it.'},
        requests=REQUESTS, tolerance=previous['tolerance'],
        acceptance='Exact row, column, per-depth and Fraction CPM invariants; independently replay every equal-depth move; before/after calls and consensus signs match fixed request; sparse+dense official/direct pairwise/tie/BH comparisons. Independent recursive raw census and exact CPM filtering agree with production product census, every table agrees with official DE, requests report scoped feasible/infeasible and witnesses. No requirement that all tiny requests be feasible. Recheck protected inputs/environment at end.',
        limits='Artificial synthetic inputs; sufficient depth-stratified operator only. Finite infeasibility applies only to enumerated tiny margins. No general emitter, realism, competitive score or hidden-GPU claim. No k/seed/amplitude/margin adjustment after results.',
        official=scoring.provenance(), dependency_sha256=previous['dependency_sha256'],
        preserve_sha256={p: rt.sha(p) for p in paths},
        biological_scored=False, holdout_scored=False, phase2_allowed=False))
    print('FROZEN', PLAN, flush=True)


def independent_moments(x):
    """Direct per-row rational accumulation; not the production column helper."""
    result = [F(0)]*x.shape[1]
    for row in x.tolist():
        depth = sum(row)
        for g, count in enumerate(row):
            result[g] += F(count, depth)
    return result


def state(rows):
    result = [0, 0]
    for r in rows:
        result[r['gene']] = r['sign'] if r['call'] else 0
    return result


def verify_bindings(plan):
    import official_scoring as scoring
    import pdex, numba_mwu
    for p, digest in plan['preserve_sha256'].items():
        assert rt.sha(p) == digest, p
    assert scoring.provenance() == plan['official']
    for module in [pdex, numba_mwu]:
        root = Path(module.__file__).parent
        assert {str(p.relative_to(root)): rt.sha(p) for p in root.rglob('*.py')} == plan['dependency_sha256'][module.__name__]


def main():
    assert not OUT.exists(), 'Never overwrite failed/completed evidence'
    plan = load(PLAN)
    verify_bindings(plan)
    receipt = dict(status='RUNNING', utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=rt.sha(PLAN), cases={}, biological_scored=False,
        holdout_scored=False, phase2_allowed=False)
    rt.save(OUT, receipt)
    f = impl.fixtures()
    pred, controls, goal = f['pred'], f['controls'], f['desired']
    before = audit.certify(pred, controls, dense=True)
    assert sum(r['call'] for r in before) == 0
    votes = np.zeros((3, 640), dtype=np.int64)
    votes[:, :300:2] = 1; votes[:, 1:300:2] = -1
    selected = rank_impl.consensus(votes, controls.mean(0) > 0, 300, 639)
    assert selected['status'] == 'SELECTED' and selected['genes'] == list(range(300))
    after = pred.copy()
    moves = []
    for a in range(0, 300, 2):
        result = impl.exchange_within_depth(after, a, a+1, goal[:, a])
        assert result['status'] == 'CONSTRUCTED'
        after = result['counts']; moves.extend(result['moves'])
    assert np.array_equal(after, goal)
    depths = pred.sum(1)
    replay = pred.copy()
    for r, s, a, b, d in moves:
        assert depths[r] == depths[s] and d > 0 and r != s and a != b
        replay[r, a] += d; replay[r, b] -= d
        replay[s, a] -= d; replay[s, b] += d
        assert np.all(replay[[r, s]][:, [a, b]] >= 0)
    assert np.array_equal(replay, after)
    assert np.array_equal(pred.sum(0), after.sum(0))
    assert np.array_equal(depths, after.sum(1))
    for d in np.unique(depths):
        assert np.array_equal(pred[depths == d].sum(0), after[depths == d].sum(0))
    mp, ma = independent_moments(pred), independent_moments(after)
    assert mp == ma == impl.fraction_sums(pred) == impl.fraction_sums(after)
    assert np.all(rank_impl.normalized(after).var(0) > 0)
    assert np.all(rank_impl.normalized(controls).var(0) > 0)
    final = audit.certify(after, controls, dense=True)
    called = [r for r in final if r['call']]
    assert [r['gene'] for r in called] == selected['genes']
    assert [r['sign'] for r in called] == selected['signs']
    assert np.array_equal(pred[:, 300:], after[:, 300:])
    receipt['cases']['main'] = dict(status='PASS', before_calls=0, after_calls=len(called),
        uncalled=len(final)-len(called), moves=len(moves), raw_depths=np.unique(depths).tolist(),
        exact_raw_and_mean_cell_cpm_margins=True, positive_control_and_pred_cpm_variance=True,
        before_sha256=audit.digest(pred), after_sha256=audit.digest(after), controls_sha256=audit.digest(controls),
        worst_called_q=max(r['q'] for r in called), uncalled_min_q=min(r['q'] for r in final if not r['call']))
    rt.save(OUT, receipt)
    print('PASS main variable-depth construction', flush=True)
    tiny, tc = f['tiny'], f['tiny_controls']
    altered = tiny.copy()
    altered[0] += [-1, 1]; altered[2] += [1, -1]
    assert np.array_equal(tiny.sum(0), altered.sum(0))
    assert np.array_equal(tiny.sum(1), altered.sum(1))
    old, new = independent_moments(tiny), independent_moments(altered)
    deltas = [(b-a)*1000000/len(tiny) for a, b in zip(old, new)]
    assert deltas == [F(-93750), F(93750)]
    assert old == impl.fraction_sums(tiny) and new == impl.fraction_sums(altered)
    tiny_rows = audit.certify(tiny, tc, dense=True)
    altered_rows = audit.certify(altered, tc, dense=True)
    receipt['cases']['raw_totals_do_not_pin_cpm'] = dict(status='PASS',
        raw_row_column_margins_equal=True, before_mean_cell_cpm=[str(v*1000000/4) for v in old],
        after_mean_cell_cpm=[str(v*1000000/4) for v in new],
        delta_mean_cell_cpm=list(map(str, deltas)),
        before_lfc=[r['lfc'] for r in tiny_rows], after_lfc=[r['lfc'] for r in altered_rows])
    rt.save(OUT, receipt)
    # Independent bounded-composition recursion and independently calculated moments.
    raw = {}
    dual = {}
    for col in audit.compositions(11, [2, 4, 8, 4]):
        x = np.array([[v, d-v] for v, d in zip(col, [2, 4, 8, 4])], dtype=np.int64)
        rows = audit.certify(x, tc)
        moment = independent_moments(x)[0]
        entry = dict(column=list(col), state=state(rows), fraction_sum=str(moment))
        raw[tuple(col)] = entry
        if moment == F(19, 8):
            dual[tuple(col)] = entry
        print('Certified tiny table', len(raw), list(col), flush=True)
    domains = {}
    for name, target, expected in [('raw_only', None, raw), ('raw_plus_cpm', F(19, 8), dual)]:
        tables = impl.census(tc, [2, 4, 8, 4], [11, 7], target)
        assert {tuple(r['column']): r for r in tables} == expected
        queries = {}
        for request in REQUESTS:
            witnesses = [r['column'] for r in tables if r['state'] == list(request)]
            queries[str(request)] = dict(status='FEASIBLE' if witnesses else 'INFEASIBLE exhaustive',
                witness=witnesses[0] if witnesses else None, witness_count=len(witnesses))
        domains[name] = dict(table_count=len(tables), requests=queries, tables=tables)
    assert dual and len(raw) > len(dual)
    unsupported_goal = np.array([0, 4, 3, 4], dtype=np.int64)
    witness = np.column_stack([unsupported_goal, tiny.sum(1)-unsupported_goal])
    assert independent_moments(witness) == independent_moments(tiny)
    assert np.array_equal(witness.sum(0), tiny.sum(0))
    assert tuple(unsupported_goal) in dual
    result = impl.exchange_within_depth(tiny, 0, 1, unsupported_goal)
    assert result['status'] == 'UNSUPPORTED depth-stratum totals' and not result['moves']
    receipt['cases']['exhaustive'] = dict(status='PASS', domains=domains,
        unsupported_feasible_witness=unsupported_goal.tolist(),
        implication='Equal-depth exchange is sufficient, not necessary; unsupported is not infeasible')
    rejects = []
    for name, args in [('fractional_goal', (tiny, 0, 1, np.ones(4))),
                       ('negative_goal', (tiny, 0, 1, np.array([-1, 0, 0, 0]))),
                       ('same_gene', (tiny, 0, 0, tiny[:, 0])),
                       ('over_capacity', (tiny, 0, 1, np.array([3, 0, 0, 0])))]:
        try:
            impl.exchange_within_depth(*args)
        except ValueError:
            rejects.append(name)
        else:
            raise AssertionError(name)
    receipt['cases']['rejections'] = dict(status='PASS', rejected=rejects)
    verify_bindings(plan)
    receipt.update(status='PASS heterogeneous synthetic mechanics only',
        completed_utc=datetime.now(timezone.utc).isoformat(), limits=plan['limits'],
        uncalled_roundoff_signs=audit.ROUNDOFF_SIGNS,
        gate_unchanged_sha256=rt.sha('artifacts/calibration_gate.json'))
    rt.save(OUT, receipt)
    print('PASS heterogeneous synthetic mechanics; no biological score', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--freeze', action='store_true')
    a = p.parse_args()
    freeze() if a.freeze else main()

"""Synthetic, depth-stratified count exchanges; no biological response inputs.

The operator is sufficient on its stated domain, not a general feasibility test.
All fixed fixtures are deterministic; no seed, amplitude or call-budget search.
"""
from fractions import Fraction
import itertools
import math

import numpy as np
from call_control_synthetic import normalized, summarize


def fraction_sums(counts):
    """Exact sums of cell fractions, distinct from raw column totals."""
    depths = counts.sum(1)
    return [sum((Fraction(int(v), int(d)) for v, d in zip(counts[:, g], depths)),
                Fraction(0)) for g in range(counts.shape[1])]


def exchange_within_depth(counts, a, b, goals):
    """Reach requested a-column using equal-depth 2x2 integer exchanges only.

    Equal-depth exchange preserves both raw totals and mean-cell CPM. A goal
    outside this sufficient domain is UNSUPPORTED, never globally infeasible.
    """
    x = np.asarray(counts).copy()
    normalized(x)
    goal = np.asarray(goals)
    if (a == b or not 0 <= a < x.shape[1] or not 0 <= b < x.shape[1]
            or goal.shape != (len(x),) or not np.issubdtype(goal.dtype, np.integer)):
        raise ValueError('Invalid pair or integer goal')
    capacity = x[:, a] + x[:, b]
    if np.any(goal < 0) or np.any(goal > capacity):
        raise ValueError('Goal exceeds pair capacity')
    depths = x.sum(1)
    groups = [np.flatnonzero(depths == d) for d in np.unique(depths)]
    if any(int(goal[rows].sum()) != int(x[rows, a].sum()) for rows in groups):
        return dict(status='UNSUPPORTED depth-stratum totals', moves=[])
    moves = []
    for rows in groups:
        low = [int(r) for r in rows if x[r, a] < goal[r]]
        high = [int(r) for r in rows if x[r, a] > goal[r]]
        i = j = 0
        while i < len(low) and j < len(high):
            r, s = low[i], high[j]
            delta = min(int(goal[r]-x[r, a]), int(x[s, a]-goal[s]))
            x[r, a] += delta; x[r, b] -= delta
            x[s, a] -= delta; x[s, b] += delta
            moves.append([r, s, a, b, delta])
            i += int(x[r, a] == goal[r])
            j += int(x[s, a] == goal[s])
    assert np.array_equal(x[:, a], goal)
    return dict(status='CONSTRUCTED', counts=x, moves=moves)


def fixtures():
    """400 predicted cells, 64 controls, 640 genes; depths vary fourfold."""
    def base(n):
        factors = np.tile([1, 2, 4, 2], n//4)
        parity = np.empty(n, dtype=np.int64)
        for f in [1, 2, 4]:
            rows = np.flatnonzero(factors == f)
            parity[rows] = np.arange(len(rows)) % 2
        x = np.empty((n, 640), dtype=np.int64)
        x[:, 0::2] = ((3+2*parity)*factors)[:, None]
        x[:, 1::2] = ((5-2*parity)*factors)[:, None]
        return x, factors, parity
    pred, factors, parity = base(400)
    controls, _, _ = base(64)
    desired = pred.copy()
    for a in range(0, 300, 2):
        pred[:, a] = (2+6*parity)*factors
        pred[:, a+1] = (6-6*parity)*factors
        desired[:, a] = (4+2*parity)*factors
        desired[:, a+1] = (4-2*parity)*factors
    tiny = np.array([[1, 1], [2, 2], [5, 3], [3, 1]], dtype=np.int64)
    # Control fractions 3/8, 1/2, 1/2, 5/8; depth8,16,32,16.
    tiny_controls = np.tile([[3, 5], [8, 8], [16, 16], [10, 6]], (8, 1))
    return dict(pred=pred, controls=controls, desired=desired,
                tiny=tiny, tiny_controls=tiny_controls)


def census(controls, depths, totals, fraction_target=None):
    """Exact bounded finite domain; optionally preserve mean-cell CPM too."""
    depths = list(map(int, depths))
    if (len(totals) != 2 or any(d <= 0 for d in depths)
            or sum(depths) != sum(totals) or any(t < 0 for t in totals)):
        raise ValueError('Invalid finite-domain margins')
    if math.prod(d+1 for d in depths) > 100000:
        raise ValueError('Domain exceeds finite enumeration limit')
    tables = []
    for column in itertools.product(*(range(d+1) for d in depths)):
        if sum(column) != totals[0]:
            continue
        moment = sum((Fraction(v, d) for v, d in zip(column, depths)), Fraction(0))
        if fraction_target is not None and moment != fraction_target:
            continue
        x = np.array([[v, d-v] for v, d in zip(column, depths)], dtype=np.int64)
        state = [0, 0]
        for r in summarize(x, controls):
            state[r['gene']] = r['sign'] if r['call'] else 0
        tables.append(dict(column=list(column), fraction_sum=str(moment), state=state))
    return tables

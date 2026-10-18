"""Clean-room synthetic DE mechanics; NOT a biological prediction recipe.

Derived from the installed official scorer and mathematical definitions only.
No stage1-lab code is imported or copied. All inputs are synthetic counts.
"""
import itertools
import math

import numpy as np


def normalized(counts):
    counts = np.asarray(counts)
    if counts.ndim != 2 or not np.issubdtype(counts.dtype, np.integer):
        raise ValueError('Expected integer count matrix')
    if np.any(counts < 0) or np.any(counts.sum(1) <= 0):
        raise ValueError('Expected nonnegative counts and positive row depths')
    return counts.astype(float) / (counts.sum(1)[:, None] / 1e6)


def summarize(pred, controls, threshold=5.0, alpha=0.05):
    """Control search-ranks + pooled ties + per-target BH, without official calls."""
    x, y = normalized(pred), normalized(controls)
    n, m = len(x), len(y)
    mean_x, mean_y = x.mean(0), y.mean(0)
    keep = mean_y > threshold
    indices = np.flatnonzero(keep)
    records = []
    for gene in indices:
        a, b = np.log1p(x[:, gene]), np.log1p(y[:, gene])
        ordered = np.sort(b)
        left = np.searchsorted(ordered, a, side='left')
        right = np.searchsorted(ordered, a, side='right')
        u = float((left + (right-left)/2).sum())
        _, frequencies = np.unique(np.concatenate([a, b]), return_counts=True)
        ties = sum(int(t)**3-int(t) for t in frequencies)
        variance = n*m/12*((n+m+1)-ties/((n+m)*(n+m-1)))
        p = 1.0 if variance <= 0 else min(1.0, math.erfc(
            (abs(u-n*m/2)-0.5)/math.sqrt(2*variance)))
        records.append(dict(gene=int(gene), u=u, ties=ties, variance=variance,
            p=p, lfc=float(np.log2((mean_x[gene]+1e-9)/(mean_y[gene]+1e-9)))))
    if records:
        order = np.argsort([r['p'] for r in records], kind='stable')
        q = 1.0
        for k in range(len(order)-1, -1, -1):
            i = int(order[k])
            q = min(q, records[i]['p']*len(order)/(k+1))
            records[i]['q'] = q
        for r in records:
            r['call'] = bool(r['q'] < alpha)
            r['sign'] = int(np.sign(r['lfc']))
    return records


def consensus(votes, eligible, k, self_gene):
    """Strict synthetic vote consensus; fixed index tie-break, no response labels."""
    votes = np.asarray(votes)
    if not np.isin(votes, [-1, 0, 1]).all():
        raise ValueError('Votes must be -1, 0, 1')
    agreement = np.abs(votes.sum(0))
    candidates = [g for g in range(votes.shape[1]) if eligible[g]
                  and g != self_gene and agreement[g] == votes.shape[0]]
    candidates.sort(key=lambda g: (-int(agreement[g]), g))
    if len(candidates) < k:
        return dict(status='INFEASIBLE eligibility/consensus', available=len(candidates))
    selected = candidates[:k]
    return dict(status='SELECTED', genes=selected,
                signs=[int(np.sign(votes[:, g].sum())) for g in selected])


def balance_pair(counts, a, b):
    """2x2 count exchanges to uniform columns, preserving both margins exactly.

    Supported domain: constant pair row depth and divisible column totals.
    This constructive toy operator does not claim general DE feasibility.
    """
    x = np.asarray(counts).copy()
    normalized(x)  # validate first
    pair_depth = x[:, a]+x[:, b]
    total = int(x[:, a].sum())
    if np.any(pair_depth != pair_depth[0]) or total % len(x):
        raise ValueError('Unsupported margins; no general infeasibility claim')
    goal = total//len(x)
    low = list(np.flatnonzero(x[:, a] < goal))
    high = list(np.flatnonzero(x[:, a] > goal))
    moves = []
    i = j = 0
    while i < len(low) and j < len(high):
        r, s = int(low[i]), int(high[j])
        delta = min(goal-int(x[r, a]), int(x[s, a])-goal)
        x[r, a] += delta; x[r, b] -= delta
        x[s, a] -= delta; x[s, b] += delta
        moves.append([r, s, a, b, delta])
        if x[r, a] == goal:
            i += 1
        if x[s, a] == goal:
            j += 1
    assert np.all(x[:, a] == goal) and np.all(x >= 0)
    return x, moves


def enumerate_two_gene(controls, row_depths, totals, desired):
    """Exact finite enumeration ONLY for a tiny two-column count matrix.

    A 0 request means uncalled; +/-1 means called with that LFC sign. Every
    feasible count table is checked, so INFEASIBLE is scoped and certified.
    """
    row_depths = [int(d) for d in row_depths]
    if len(totals) != 2 or len(desired) != 2 or sum(row_depths) != sum(totals):
        return dict(status='INFEASIBLE margins', checked=0)
    if math.prod(d+1 for d in row_depths) > 100000:
        raise ValueError('Enumeration limit; no feasibility conclusion')
    feasible = []
    for column in itertools.product(*(range(d+1) for d in row_depths)):
        if sum(column) != totals[0]:
            continue
        x = np.array([[v, d-v] for v, d in zip(column, row_depths)], dtype=np.int64)
        rows = summarize(x, controls)
        state = [0, 0]
        for r in rows:
            state[r['gene']] = r['sign'] if r['call'] else 0
        feasible.append(dict(column=list(column), state=state))
    witnesses = [r['column'] for r in feasible if r['state'] == list(desired)]
    return dict(status='FEASIBLE' if witnesses else 'INFEASIBLE exhaustive',
                checked=len(feasible), witness=witnesses[0] if witnesses else None,
                tables=feasible)


def fixtures():
    """Frozen deterministic fixtures, no seeds or parameter search."""
    # Each pair has constant depth8; the first150 pairs have means5/3.
    controls = np.full((64, 640), 4, dtype=np.int64)
    pred = np.full((400, 640), 4, dtype=np.int64)
    for pair in range(150):
        a, b = 2*pair, 2*pair+1
        pred[:200, a] = 2; pred[:200, b] = 6
        pred[200:, a] = 8; pred[200:, b] = 0
    unequal_p = np.array([[0, 1, 1, 6], [1, 1, 2, 12], [2, 4, 2, 24],
                         [0, 0, 4, 4], [4, 0, 0, 12]], dtype=np.int64)
    unequal_c = np.array([[0, 0, 2, 6], [0, 1, 3, 12], [0, 2, 6, 24],
                         [0, 0, 2, 6], [0, 1, 3, 12], [0, 2, 6, 24],
                         [0, 0, 2, 6]], dtype=np.int64)
    # The zero-control gene is filtered even if predicted expression is positive.
    return dict(main=(pred, controls), unequal_depth_filter=(unequal_p, unequal_c),
                all_tied=(np.full((9, 4), 2, dtype=np.int64),
                          np.full((7, 4), 2, dtype=np.int64)),
                tiny=(np.array([[1, 1], [1, 1], [2, 0], [2, 0]], dtype=np.int64),
                      np.full((16, 2), 1, dtype=np.int64)))

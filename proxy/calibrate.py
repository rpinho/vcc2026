"""Calibrate the LOCO proxy against official validation scores (PIE_AND_EVIDENCE_PLAN fix #3).

Pass rule (fixed in the task brief BEFORE any proxy score was looked at):
  Spearman(proxy macro overall, official overall) >= 0.8 over the scored recipes  AND  proxy ranks v4a above v4b.
Macro = mean over destinations scored for every recipe. Also reports Kendall, per-destination and per-component
Spearman, and the leave-one-recipe-out error of official ~ a + b * proxy.
Official v5a/v5b are read from model_v5/submit_v5_tonight.log when present (arrive after 2026-10-18 00:03 UTC).
"""
import json, re, sys, itertools
from pathlib import Path
import numpy as np
from scipy import stats
P = Path(__file__).resolve().parent; V = P.parent
SIX = ['pds_cosine', 'expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae', 'de_wilcoxon_direction_fidelity_yield_raw', 'de_wilcoxon_direction_reach_raw', 'de_wilcoxon_sig_jaccard']
SHORT = dict(zip(SIX, ['pds', 'mse', 'nmae', 'fid', 'reach', 'jac']))
OFFICIAL = {'v1': dict(overall=.0877, pds=.517, mse=0, nmae=-.006, fid=-.067, reach=.099, jac=-.017),
            'v2b': dict(overall=.0617, pds=.446, mse=0, nmae=.042, fid=-.188, reach=.092, jac=-.020),
            'v3b': dict(overall=-.2482, pds=.050, mse=0, nmae=.014, fid=-1.539, reach=.028, jac=-.042),
            'v4a': dict(overall=.1438, pds=.678, mse=0, nmae=.129, fid=-.019, reach=.076, jac=-.001),
            'v4b': dict(overall=.1243, pds=.658, mse=0, nmae=.070, fid=-.025, reach=.042, jac=-.000),
            # 2026-10-18 official (OFFICIAL_SCORES.md); the log parser below never matched the CLI's output format.
            'v5a': dict(overall=.1499, pds=.678, mse=.078, nmae=.094, fid=-.020, reach=.071, jac=-.001),
            'v5b': dict(overall=.1528, pds=.680, mse=.082, nmae=.109, fid=-.031, reach=.077, jac=-.000)}


def official_v5():
    log = V / 'model_v5/submit_v5_tonight.log'
    out = {}
    if not log.exists():
        return out
    text = log.read_text(errors='replace')
    for arm in ['v5a', 'v5b']:
        # accept JSON-ish or table-ish status lines naming the arm and an overall score
        m = re.findall(rf'{arm}[^\n]*?(?:overall|score)["\s:=]+(-?\d+\.\d+)', text, flags=re.I)
        if m:
            out[arm] = dict(overall=float(m[-1]))
    return out


def main():
    off = dict(OFFICIAL); off.update(official_v5())
    rows = {}
    for f in sorted((P / 'artifacts/scores').glob('*/*.json')):
        if f.name == 'anchors.json' or '_fill' in f.name:
            continue
        r = json.loads(f.read_text()); rows.setdefault(r['recipe'], {})[r['dest']] = r
    recipes = [r for r in rows if r in off]
    dests = sorted(set.intersection(*[set(rows[r]) for r in recipes])) if recipes else []
    res = dict(official_used={r: off[r] for r in recipes}, destinations=dests, recipes=recipes)
    if len(recipes) < 3 or not dests:
        res['verdict'] = 'INSUFFICIENT: need >= 3 recipes with official scores on a common destination set'
        print(json.dumps(res, indent=1)); return res
    macro = {r: float(np.mean([rows[r][d]['overall'] for d in dests])) for r in recipes}
    x = np.array([macro[r] for r in recipes]); y = np.array([off[r]['overall'] for r in recipes])
    rho = stats.spearmanr(x, y).statistic; tau = stats.kendalltau(x, y).statistic
    res.update(macro_proxy=macro, spearman=float(rho), kendall=float(tau), pearson=float(np.corrcoef(x, y)[0, 1]))
    res['per_destination_spearman'] = {d: float(stats.spearmanr([rows[r][d]['overall'] for r in recipes], y).statistic) for d in dests}
    comp = {}
    for m in SIX:
        k = SHORT[m]
        if all(k in off[r] for r in recipes):
            px = [np.mean([rows[r][d]['scaled'][m] for d in dests]) for r in recipes]; oy = [off[r][k] for r in recipes]
            comp[k] = dict(proxy=dict(zip(recipes, map(float, px))), official=dict(zip(recipes, oy)),
                           spearman=float(stats.spearmanr(px, oy).statistic) if len(set(oy)) > 1 else None)
    res['per_component'] = comp
    if len(recipes) >= 3:
        err = []
        for i in range(len(recipes)):
            m = np.ones(len(recipes), bool); m[i] = False
            b, a = np.polyfit(x[m], y[m], 1); err.append(abs(a + b * x[i] - y[i]))
        res['loo_abs_error'] = dict(zip(recipes, map(float, err))); res['loo_mae'] = float(np.mean(err))
    v4 = ('v4a' in macro and 'v4b' in macro and macro['v4a'] > macro['v4b'])
    res['v4a_above_v4b'] = bool(v4)
    res['verdict'] = 'PASS' if (rho >= .8 and v4) else 'FAIL'
    res['rule'] = 'Spearman >= 0.8 AND proxy v4a > v4b (fixed before scoring)'
    (P / 'artifacts/calibration.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1)); return res


if __name__ == '__main__':
    main()

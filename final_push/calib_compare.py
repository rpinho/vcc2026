"""Does the H1_rich_search harness order our official DE history correctly? Compares local scaled deltas with the
official validation deltas for fidelity, reach, nmae and jaccard. Usage: python calib_compare.py [task]"""
import sys, json
from pathlib import Path
HERE = Path(__file__).resolve().parent
task = sys.argv[1] if len(sys.argv) > 1 else 'H1_rich_search'
S = json.loads((HERE / 'artifacts/de_check' / task / 'calibration_summary.json').read_text())
K = {'fid': 'de_wilcoxon_direction_fidelity_yield_raw', 'reach': 'de_wilcoxon_direction_reach_raw',
     'nmae': 'de_wilcoxon_lfc_nmae', 'jac': 'de_wilcoxon_sig_jaccard', 'pds': 'pds_cosine'}
# official validation components (OFFICIAL_SCORES.md / PUSH_1021.md)
OFF = {'v5b': dict(fid=-0.031, reach=0.077, nmae=0.109, jac=-0.000, pds=0.680),
       'v8a': dict(fid=-0.064, reach=0.162, nmae=0.102, jac=-0.015, pds=0.680),
       'v9a': dict(fid=-0.068, reach=0.181, nmae=0.108, jac=-0.016, pds=0.680),
       'v8b': dict(fid=-0.060, reach=0.178, nmae=-0.059, jac=-0.015, pds=0.708)}
loc = {arm: {k: S[arm]['scaled'][v] for k, v in K.items()} for arm in S}
print('local scaled:'); [print(f"  {a:10s}", {k: round(v, 4) for k, v in d.items()}) for a, d in loc.items()]
contrasts = [('v9a_885', 'v5b', 'v9a', 'v5b'), ('v9a_694', 'v5b', 'v8a', 'v5b'), ('v8b_like', 'v9a_885', 'v8b', 'v8a')]
noise = 0.01   # treat |official delta| < 0.01 as 'flat'
total = 0; hit = 0
for la, lb, oa, ob in contrasts:
    if la not in loc or lb not in loc:
        continue
    print(f"\n{la} - {lb}  vs official {oa} - {ob}")
    for k in ['fid', 'reach', 'nmae', 'jac', 'pds']:
        dl = loc[la][k] - loc[lb][k]; do = OFF[oa][k] - OFF[ob][k]
        so = 'flat' if abs(do) < noise else ('up' if do > 0 else 'down')
        sl = 'flat' if abs(dl) < noise else ('up' if dl > 0 else 'down')
        ok = (so == sl) or (so == 'flat' and abs(dl) < 0.03)
        if k != 'pds':
            total += 1; hit += ok
        print(f"  {k:5s} local {dl:+.4f} ({sl:4s})   official {do:+.4f} ({so:4s})   {'OK' if ok else 'MISS'}")
print(f"\nDE-metric sign agreement: {hit}/{total}")

"""Supplement: remove only the conservative sealed-target-name exclusion.

This admits more CD4 targets without using any H1 response source. A low bound
in this larger universe implies impossibility in the original narrower scope.
The supplemental estimand is frozen after the primary count-only certificate;
no expression or model outcomes enter either calculation.
"""
import argparse
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt
import numpy as np
import pandas as pd

PLAN = Path('artifacts/CD4_capacity_superset_plan_round44.json')
OUT = Path('artifacts/CD4_capacity_superset_audit_round44.json')
TABLE = 'artifacts/CD4_source_capacity_targets_round44.tsv'
PRIMARY = 'artifacts/CD4_source_capacity_audit_round44.json'


def freeze():
    assert not PLAN.exists() and not OUT.exists()
    rt.save(PLAN, dict(status='FROZEN supplemental count-only estimand after primary counts',
        utc=datetime.now(timezone.utc).isoformat(),
        preserve_sha256={p: rt.sha(p) for p in [__file__, TABLE, PRIMARY,
            'artifacts/CD4_source_capacity_plan_round44.json',
            'artifacts/CD4_source_capacity_saved_check_round44.json']},
        scope='All CD4 targets with both prepared legacy and raw GWPS>=200 and either prepared or cached HCT116>=20, including H1 sealed target names. This is an eligibility relaxation only: no H1 response is accessed or source activated.',
        estimand='Fixed79 target arbitrary common positive24-lane quota necessary upper bound, unlimited controls; uniform common quota ceiling; count all-lane8 targets. Scalar and NumPy order statistics must agree. Compare192 and400 only; no target selection or lower-depth design.',
        primary_scope_limitation='Primary audit excludes H1 sealed target identities conservatively; those targets need not be excluded when H1 response sources are entirely unused. This supplement checks robustness to removing that extra restriction.',
        expression_payloads_decoded=False, holdout_scored=False, phase2_allowed=False))
    print('FROZEN supplemental superset bound', flush=True)


def main():
    assert not OUT.exists()
    plan = rt.json.loads(PLAN.read_text())
    for p, h in plan['preserve_sha256'].items():
        assert rt.sha(p) == h, p
    primary = rt.json.loads(Path(PRIMARY).read_text())
    assert primary['status'].startswith('PASS') and rt.sha(TABLE) == primary['table_sha256']
    t = pd.read_csv(TABLE, sep='\t', index_col='target')
    lanes = [c for c in t.columns if c.startswith('CD4i_R2L')]
    assert len(lanes) == 24
    s = t.loc[(t.legacy_cells >= 200) & (t.raw_cells >= 200) &
              ((t.prepared_hct_cells >= 20) | (t.cached_hct_cells >= 20))]
    assert len(s) >= 79
    upper = {b: int(sorted(s[b], reverse=True)[78]) for b in lanes}
    assert list(upper.values()) == np.sort(s[lanes].to_numpy(), axis=0)[-79].tolist()
    q = int(sorted(s[lanes].min(axis=1), reverse=True)[78])
    assert q == int(np.sort(s[lanes].to_numpy().min(1))[-79])
    assert all(upper[b] >= primary['scopes']['union_two_family_potential']['per_lane_upper'][b] for b in lanes)
    rt.save(OUT, dict(status='PASS relaxed-superset necessary bound only',
        utc=datetime.now(timezone.utc).isoformat(), plan_sha256=rt.sha(PLAN),
        targets=len(s), additional_sealed_name_targets=int(s.sealed_H1_target.sum()),
        all_lane8_targets=int((s[lanes].min(axis=1) >= 8).sum()),
        uniform_depth_ceiling=24*q, arbitrary_common_depth_upper_bound=sum(upper.values()),
        per_lane_upper=upper,
        impossible192=bool(sum(upper.values()) < 192 or min(upper.values()) == 0),
        impossible400=bool(sum(upper.values()) < 400 or min(upper.values()) == 0),
        limits='Upper bound is necessary, not attainable depth. Cached selection only, not full public releases. This relaxes the primary certificate and activates no source, task or target panel.',
        expression_payloads_decoded=False, new_task_registered=False,
        network_reads=False, holdout_scored=False, phase2_allowed=False))
    print('PASS relaxed superset', len(s), 'targets;', sum(upper.values()),
          'cell arbitrary-mixture ceiling; uniform', 24*q, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--freeze', action='store_true')
    freeze() if p.parse_args().freeze else main()

"""Serial fixed Phase1 recalibration in a distinct source/cache namespace."""
from datetime import datetime, timezone
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import gc
import time
import runtime as rt
import hct116_revision_round54 as revision


def score(task, recipe, plan):
    import benchmark
    import bounded_de
    import bounded_score
    import stat_cache
    import official_scoring as scoring
    from audit_phase1_calibration import exact_scaled
    output = revision.ROOT/'scores'/task
    output.mkdir(parents=True, exist_ok=True)
    path = output/(recipe+'_s42.json')
    anchors_path = Path('artifacts/scores')/task/'anchors.json'
    anchors = revision.load(anchors_path)
    if path.exists():
        record = revision.load(path)
        assert record['plan_sha256'] == rt.sha(revision.PLAN)
        assert record['audit_sha256'] == rt.sha(revision.ROOT/'audit.json')
        exact_scaled(record, anchors)
        return record
    start = time.monotonic()
    entry, ctrl, real = benchmark.setup(task)
    assert entry == plan['tasks'][task] and entry['split'] == 'search'
    assert anchors['n_splits'] == 5 and anchors['real_sha256'] == scoring.matrix_id(real)
    assert anchors['driver_sha256'] == rt.sha('bounded_score.py')
    assert anchors['scorer_sha256'] == rt.sha('official_scoring.py')
    scoring.scale(anchors['baseline'], anchors['baseline'], anchors['replicate'])
    bounded_de.install()
    stat_cache.install()
    stat_cache.namespace += ('|HCT116_source509_r54|'+rt.sha(revision.PLAN)).encode()
    seen = []
    def iterator():
        for target, matrix in revision.predictions(ctrl, task, recipe, plan):
            rt.guard(512*1024**2)
            seen.append(target)
            if len(seen)%25 == 1:
                print('Revision emission',task,recipe,len(seen),'/',len(entry['targets']),flush=True)
            yield target,matrix
    tidy,digest = bounded_score.score_iter(real,iterator(),
        scoring.config(cache_real=str(output/'real_cache')),release=benchmark.release_models)
    assert seen == entry['targets']
    raw = benchmark.aggregate(tidy)
    scaled, overall = scoring.scale(raw, anchors['baseline'],anchors['replicate'])
    record = dict(revision_id=plan['revision_id'],task=task,recipe=recipe,seed=42,
        effective_seed=20260910,raw=raw,scaled=scaled,overall=overall,
        wall_s=time.monotonic()-start,prediction_sha256=digest,
        plan_sha256=rt.sha(revision.PLAN),audit_sha256=rt.sha(revision.ROOT/'audit.json'),
        anchors_sha256=rt.sha(anchors_path),reference_sha256=anchors['real_sha256'],
        selection_sha256=entry['selection_sha256'],n_targets=79,cells_pred=400,
        status='Revised source-population recipe port; no exact official score pairing or promotion',
        utc=datetime.now(timezone.utc).isoformat())
    exact_scaled(record,anchors)
    tidy.write_csv(output/(recipe+'_tidy.csv'))
    rt.save(path,record)
    print('Revision SCORE',task,recipe,overall,scaled,flush=True)
    del real,ctrl
    gc.collect()
    return record


def main():
    plan = revision.check_plan()
    audit = revision.load(revision.ROOT/'audit.json')
    assert audit['status'].startswith('PASS')
    assert audit['plan_sha256'] == rt.sha(revision.PLAN)
    assert audit['verifier_sha256'] == rt.sha('verify_hct116_revision_round54.py')
    assert __import__('official_scoring').provenance() == plan['official_environment']
    assert not (revision.ROOT/'results.json').exists()
    records = {}
    from audit_phase1_calibration import exact_scaled, correlations, OFFICIAL
    for task in revision.TASKS:
        records[task] = {}
        for recipe in ['v1','v2b','v3b']:
            path = Path('artifacts/scores')/task/(recipe+'_s42.json')
            record = revision.load(path)
            assert record['implementation_sha256'] == plan['original_implementation_sha256'][task]
            exact_scaled(record,revision.load(Path('artifacts/scores')/task/'anchors.json'))
            records[task][recipe] = dict(overall=record['overall'],path=str(path),sha256=rt.sha(path),
                reused_unchanged_legacy_result=True)
        for recipe in revision.RECIPES:
            record = score(task,recipe,plan)
            path = revision.ROOT/'scores'/task/(recipe+'_s42.json')
            records[task][recipe] = dict(overall=record['overall'],path=str(path),sha256=rt.sha(path),
                reused_unchanged_legacy_result=False)
    macro = {r:float(sum(Fraction(records[t][r]['overall']) for t in revision.TASKS)/2)
             for r in ['v1','v2b','v3b','v4a','v4b']}
    pairs = [dict(a=a,b=b,matches=(macro[a]-macro[b])*(OFFICIAL[a]-OFFICIAL[b]) > 0)
             for a,b in combinations(macro,2)]
    assert all(__import__('provenance').identity(t) == h for t,h in plan['original_implementation_sha256'].items())
    revision.check_plan()
    rt.save(revision.ROOT/'results.json',dict(status='Completed fixed source-adequacy recalibration; independent final result audit pending',
        utc=datetime.now(timezone.utc).isoformat(),plan_sha256=rt.sha(revision.PLAN),
        runner_sha256=rt.sha(__file__),tasks=records,macro=macro,correlation=correlations(macro),
        all_pairs=pairs,all_pairs_rank_consistent=all(p['matches'] for p in pairs),
        without_v3b_diagnostic_only=correlations({r:v for r,v in macro.items() if r != 'v3b'}),
        historical_gate_preserved=True,phase2_allowed=False,holdout_scored=False,
        limitations='Dependent n5 recipe ports, two79-target panels; CD4 small fidelity scale. Source revision does not reproduce hidden submitted source bytes. Exact external baselines still unavailable. No official forecast.'))
    print('Completed revised calibration',macro,correlations(macro),flush=True)


if __name__ == '__main__':
    main()

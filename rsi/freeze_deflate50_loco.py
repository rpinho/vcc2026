"""Freeze one description-based, available-input priority4 baseline."""
from pathlib import Path
from datetime import datetime, timezone
import runtime as rt
import atlas_loco_baseline as atlas
import deflate50_loco_baseline as baseline


def main():
    assert not baseline.SPEC.exists()
    parent = atlas.read_spec()
    assert rt.json.loads(Path('artifacts/calibration_gate.json').read_text())['status'] == 'CLOSED'
    descriptions = Path('../competitors/raw/leaderboard-descriptions-2026-10-17.txt')
    selected = [line for line in descriptions.read_text().splitlines()
                if line.startswith('376 | Hao Tian |') or line.startswith('423 | Becker Lab |')]
    assert len(selected) == 2
    rt.save('artifacts/deflate50_primary_descriptions_round45.json', {
        'source': str(descriptions), 'source_sha256': rt.sha(descriptions), 'author_descriptions': selected,
        'scope': 'Cached primary author statements; Becker commit c6d81cc code absent from local snapshots. No new public refresh or exact numerical release claim.'})
    spec = dict(parent)
    spec.update({'utc': datetime.now(timezone.utc).isoformat(), 'baseline_id': baseline.BASELINE_ID,
        'axis_seed': 20260910, 'parent_spec_sha256': rt.sha(atlas.SPEC),
        'recipe': dict(parent['recipe'], deflation={'k':50, 'eta':1,
            'space':'column-centered input-control cell log1p CP50K',
            'location':'after add_cd4_family / before desired_mean amplitude and clip',
            'channel':'bulk_delta only', 'max_controls':4000, 'PCA_dtype':'float32',
            'randomized_svd_iterations':5, 'PCA_seed':20260910}),
        'adaptations': parent['adaptations'] + [
            'Description-based reconstruction, not exact Becker c6d81cc code or official0.1669 reproduction.',
            'Retain released AtlasShift amplitudes0.6/0.3; Hao V11 amplitudes0.6/0.6 not adopted.',
            'Bounded float32 randomized PCA convention matches archived local v5b helper; seed20260910 and max4000 fixed without outcomes. All current input controls fit under cap.',
            'Project fused bulk_delta before released desired_mean; released mean-CPM branch, promoter, template and integer emitter unchanged. No B/claim calibration or DE call gate.'],
        'external_official_label':None, 'exact_full_release_reproduction':False,
        'official_environment': parent['official_environment']})
    tracked = ['deflate50_loco_baseline.py','freeze_deflate50_loco.py','verify_deflate50_loco.py','run_deflate50_loco.py']
    spec['implementation_sha256'] = dict(parent['implementation_sha256'])
    spec['implementation_sha256'].update({p: rt.sha(p) for p in tracked})
    extra = [str(atlas.SPEC), 'artifacts/atlas_loco_audit.json', 'artifacts/calibration_gate.json',
        'artifacts/deflate50_primary_descriptions_round45.json',
        'artifacts/reproduction_inputs/v5b_round33/v5.py']
    for task in baseline.TASKS:
        extra += [f'artifacts/tasks/{task}.h5ad']
        extra += [f'artifacts/scores/{task}/{r}_s42.json' for r in ['v1','v2b','v3b','v4a','v4b']]
    spec['protocol_sha256'] = dict(parent['protocol_sha256'])
    spec['protocol_sha256'].update({p:rt.sha(p) for p in extra})
    rt.save(baseline.SPEC,spec)
    print('FROZEN', baseline.SPEC, rt.sha(baseline.SPEC), flush=True)


if __name__ == '__main__': main()

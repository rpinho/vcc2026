"""Activate the round-13 count-only protocol after its independent count audit.

No response scores are read. Originals and expected replacements are archived
before mutation, so a partial activation can be resumed without changing intent.
"""
import copy
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt

NAME = 'CD4_depth120_search'
OLD = 'CD4_matched_search'
REVISION = Path('artifacts/protocol_revisions/CD4_depth120')
FILES = ['artifacts/tasks/manifest.json', 'artifacts/benchmark_roles.json',
         'artifacts/anchor_policy.json']


def main():
    plan_path = Path('artifacts/benchmark_depth120_revision.json')
    plan = rt.json.loads(plan_path.read_text())
    audit_path = Path('artifacts/CD4_depth120_materialization_audit.json')
    audit = rt.json.loads(audit_path.read_text())
    assert audit['status'].startswith('PASS independent count/')
    assert audit['task'] == NAME and audit['holdout_scored'] is False
    extension_path = Path(f'artifacts/tasks/{NAME}_manifest_extension.json')
    assert rt.sha(extension_path) == audit['manifest_extension_sha256']
    extension = rt.json.loads(extension_path.read_text())
    assert rt.sha(f'artifacts/tasks/{NAME}.h5ad') == audit['sha256']
    assert rt.sha(f'artifacts/tasks/{NAME}_blueprint.json') == audit['blueprint_sha256']
    assert plan['primary_search'] == ['H1_rich_search', NAME]
    assert plan['anchors'] == {'H1_rich_search': 5, NAME: 5}
    # The historical archive/index are not rewritten by an active revision.
    index = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    assert rt.sha(index['archive']) == index['archive_sha256']
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    REVISION.mkdir(parents=True, exist_ok=True)
    intent_path = REVISION/'intent.json'
    if not intent_path.exists():
        originals = {}
        for path in FILES:
            source = Path(path)
            archive = REVISION/(source.name + '.before')
            data = source.read_bytes()
            if archive.exists():
                assert archive.read_bytes() == data
            else:
                archive.write_bytes(data)
            originals[path] = dict(archive=str(archive), sha256=rt.sha(archive))
        manifest, roles, policy = [rt.json.loads(Path(originals[p]['archive']).read_text()) for p in FILES]
        assert NAME not in manifest and OLD in roles['primary_search']
        manifest[OLD]['eligible_for_primary_calibration'] = False
        manifest[OLD]['superseded_by'] = NAME
        manifest[OLD]['qualification'] = 'INVALID six-metric scale; diagnostic only, never retried or included in calibration macro'
        manifest[NAME] = copy.deepcopy(extension[NAME])
        # The immutable blueprint inherits the old three-split field. The later
        # frozen round-13 protocol explicitly supersedes it before any new score.
        manifest[NAME]['replicate_splits'] = plan['anchors'][NAME]
        manifest[NAME]['anchor_policy_revision'] = str(plan_path)
        roles['primary_search'] = plan['primary_search']
        roles['diagnostic_search'] = list(dict.fromkeys(roles['diagnostic_search'] + [OLD]))
        roles['cd4_depth120_revision'] = dict(plan=str(plan_path),
            materialization_audit=str(audit_path), policy=plan['policy'])
        roles['cd4_final_count_design'] = '300 unchanged targets x 120 actual cells (5/lane); unchanged 1920 input and 768 reference controls; five anchor splits'
        roles['unresolved'] = 'New CD4 six-metric scales and two-context all-pair ranks UNVERIFIED; exact public AtlasShift/Sidechain reproductions unresolved'
        policy[NAME] = plan['anchors'][NAME]
        policy['depth120_revision'] = str(plan_path)
        replacements = {}
        for path, value in zip(FILES, [manifest, roles, policy]):
            expected = REVISION/(Path(path).name + '.after')
            rt.save(expected, value)
            replacements[path] = dict(archive=str(expected), sha256=rt.sha(expected))
        rt.save(intent_path, dict(utc=datetime.now(timezone.utc).isoformat(),
            plan_sha256=rt.sha(plan_path), audit_sha256=rt.sha(audit_path),
            before=originals, after=replacements))
    intent = rt.json.loads(intent_path.read_text())
    assert intent['plan_sha256'] == rt.sha(plan_path)
    assert intent['audit_sha256'] == rt.sha(audit_path)
    for path in FILES:
        before, after = intent['before'][path], intent['after'][path]
        assert rt.sha(before['archive']) == before['sha256']
        assert rt.sha(after['archive']) == after['sha256']
        assert rt.sha(path) in [before['sha256'], after['sha256']], path
    for path in FILES:
        expected = Path(intent['after'][path]['archive'])
        temporary = Path(path + '.activation.tmp')
        temporary.write_bytes(expected.read_bytes())
        temporary.replace(path)
        assert rt.sha(path) == intent['after'][path]['sha256']
    rt.save(REVISION/'activation.json', dict(status='ACTIVE count-audited protocol; scales/ranks UNVERIFIED',
        utc=datetime.now(timezone.utc).isoformat(), intent_sha256=rt.sha(intent_path),
        script_sha256=rt.sha(__file__), before=intent['before'], after=intent['after'],
        materialization_audit_sha256=rt.sha(audit_path), holdout_scored=False,
        frozen_metadata_members_unchanged=len(index['entries'])))
    print('PASS activated frozen CD4 depth120 protocol; historical originals archived; scales/ranks UNVERIFIED', flush=True)


if __name__ == '__main__':
    main()
